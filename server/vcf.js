import { gunzipSync } from "node:zlib";
import { badRequest } from "./errors.js";
import { normalizeChromosome } from "./database.js";
import { config } from "./config.js";

const MAX_WARNINGS = 100;
const integerPattern = /^\d+$/;
const sequencePattern = /^[ACGTN]+$/;
const buildAliases = [
  { build: "GRCh38", pattern: /(?:GRCh38|hg38|GCF_000001405\.2[6-9])/i },
  { build: "GRCh37", pattern: /(?:GRCh37|hg19|GCF_000001405\.25)/i }
];

const addWarning = (warnings, line, code, message) => {
  if (warnings.length < MAX_WARNINGS) warnings.push({ line, code, message });
};

function decodeUpload(file) {
  const looksGzip =
    file.originalname?.toLowerCase().endsWith(".gz") ||
    (file.buffer[0] === 0x1f && file.buffer[1] === 0x8b);
  try {
    return (looksGzip
      ? gunzipSync(file.buffer, { maxOutputLength: config.maxUploadBytes })
      : file.buffer).toString("utf8");
  } catch (error) {
    if (error?.code === "ERR_BUFFER_TOO_LARGE") {
      throw badRequest("DECOMPRESSED_FILE_TOO_LARGE", "The decompressed VCF exceeds the configured upload limit.");
    }
    throw badRequest("INVALID_GZIP", "The uploaded gzip file could not be decompressed.");
  }
}

function parseInfo(value) {
  const info = new Map();
  if (!value || value === ".") return info;
  for (const part of value.split(";")) {
    const separator = part.indexOf("=");
    if (separator === -1) info.set(part, true);
    else info.set(part.slice(0, separator), part.slice(separator + 1));
  }
  return info;
}

function parseNullableNumber(value) {
  if (value === undefined || value === "." || value === "") return null;
  const number = Number(value);
  return Number.isFinite(number) ? number : null;
}

export function inspectVcfSamples(text) {
  const header = text.replace(/^\uFEFF/, "").split(/\r?\n/).find((line) => line.startsWith("#CHROM"));
  if (!header) return [];
  return header.split("\t").slice(9).filter(Boolean);
}

function detectGenomeBuild(referenceValues, contigAssemblies) {
  const candidates = [...referenceValues, ...contigAssemblies];
  const builds = new Set();
  for (const value of candidates) {
    const match = buildAliases.find((alias) => alias.pattern.test(value));
    if (match) builds.add(match.build);
  }
  return builds.size === 1 ? [...builds][0] : null;
}

export function parseVcfUpload(file, requestedSample) {
  if (!file?.buffer?.length) {
    throw badRequest("VCF_REQUIRED", "Choose a non-empty VCF or VCF.GZ file.");
  }
  const text = decodeUpload(file).replace(/^\uFEFF/, "");
  if (text.includes("\0")) {
    throw badRequest("INVALID_ENCODING", "VCF files must be UTF-8 text, optionally gzip-compressed.");
  }

  const lines = text.split(/\r?\n/);
  const warnings = [];
  let fileFormat = null;
  const referenceValues = [];
  const contigAssemblies = [];
  let header = null;
  let sampleNames = [];
  let selectedSample = null;
  let sampleColumn = -1;
  const records = [];
  let dataLineCount = 0;

  for (let index = 0; index < lines.length; index += 1) {
    const lineNumber = index + 1;
    const line = lines[index];
    if (!line.trim()) continue;
    if (line.startsWith("##fileformat=")) fileFormat = line.slice(13).trim();
    if (line.startsWith("##reference=")) referenceValues.push(line.slice(12).trim());
    if (line.startsWith("##contig=<")) {
      const assembly = line.match(/(?:^|,)assembly=([^,>]+)/i)?.[1];
      if (assembly) contigAssemblies.push(assembly.replace(/^"|"$/g, ""));
    }
    if (line.startsWith("##")) continue;

    if (line.startsWith("#CHROM")) {
      header = line.split("\t");
      if (header.length < 8) {
        throw badRequest("INVALID_VCF_HEADER", "The #CHROM header must contain at least 8 tab-separated columns.");
      }
      sampleNames = header.slice(9);
      if (requestedSample) {
        const sampleIndex = sampleNames.indexOf(requestedSample);
        if (sampleIndex === -1) {
          throw badRequest("SAMPLE_NOT_FOUND", `Sample "${requestedSample}" is not present in this VCF.`, { sampleNames });
        }
        selectedSample = requestedSample;
        sampleColumn = sampleIndex + 9;
      } else if (sampleNames.length) {
        selectedSample = sampleNames[0];
        sampleColumn = 9;
      }
      continue;
    }
    if (line.startsWith("#")) continue;
    if (!header) {
      throw badRequest("MISSING_VCF_HEADER", "VCF data appeared before a #CHROM header.");
    }

    dataLineCount += 1;
    const columns = line.split("\t");
    if (columns.length < 8) {
      addWarning(warnings, lineNumber, "MALFORMED_RECORD", "Record has fewer than 8 columns and was skipped.");
      continue;
    }

    const [chromosomeRaw, positionRaw, id, refRaw, altRaw, qualRaw, filterRaw, infoRaw, formatRaw] = columns;
    if (!integerPattern.test(positionRaw) || Number(positionRaw) < 1) {
      addWarning(warnings, lineNumber, "INVALID_POSITION", "Record position is invalid and was skipped.");
      continue;
    }
    const ref = refRaw.toUpperCase();
    const alts = altRaw === "." ? [] : altRaw.split(",").map((alt) => alt.toUpperCase());
    if (!ref || ref === "." || !sequencePattern.test(ref) || alts.some((alt) => !alt)) {
      addWarning(warnings, lineNumber, "INVALID_ALLELE", "Record alleles are invalid and were skipped.");
      continue;
    }
    for (const alt of alts) {
      if (!sequencePattern.test(alt)) {
        addWarning(warnings, lineNumber, "UNSUPPORTED_ALLELE", `ALT "${alt}" is symbolic or complex and cannot be matched by this service.`);
      }
    }

    const formatKeys = formatRaw && formatRaw !== "." ? formatRaw.split(":") : [];
    const sampleValues = sampleColumn >= 0 && columns[sampleColumn]
      ? columns[sampleColumn].split(":")
      : [];
    const sample = new Map(formatKeys.map((key, fieldIndex) => [key, sampleValues[fieldIndex]]));
    const genotypeRaw = sample.get("GT") ?? null;
    let alleleIndexes = [];

    if (genotypeRaw && genotypeRaw !== "." && genotypeRaw !== "./." && genotypeRaw !== ".|.") {
      const tokens = genotypeRaw.split(/[\/|]/);
      if (tokens.every((token) => token === "." || integerPattern.test(token))) {
        alleleIndexes = tokens.filter((token) => token !== ".").map(Number);
        if (alleleIndexes.some((allele) => allele > alts.length)) {
          addWarning(warnings, lineNumber, "INVALID_GENOTYPE", "Genotype references an ALT allele that does not exist.");
          alleleIndexes = [];
        }
      } else {
        addWarning(warnings, lineNumber, "INVALID_GENOTYPE", "Genotype syntax is unsupported.");
      }
    }

    const info = parseInfo(infoRaw);
    const depth = parseNullableNumber(sample.get("DP") ?? info.get("DP"));
    const genotypeQuality = parseNullableNumber(sample.get("GQ"));
    const quality = parseNullableNumber(qualRaw);
    if (qualRaw !== "." && quality === null) {
      addWarning(warnings, lineNumber, "INVALID_QUALITY", "QUAL is not numeric; it was treated as missing.");
    }

    records.push({
      line: lineNumber,
      chromosome: normalizeChromosome(chromosomeRaw),
      position: Number(positionRaw),
      id: id === "." ? null : id,
      ref,
      alts,
      filter: filterRaw,
      sampleFilter: sample.get("FT") ?? null,
      quality,
      genotypeQuality,
      depth,
      genotype: genotypeRaw,
      phased: genotypeRaw?.includes("|") ?? false,
      alleleIndexes
    });
  }

  if (!header) throw badRequest("MISSING_VCF_HEADER", "No #CHROM header was found.");
  if (dataLineCount === 0) throw badRequest("EMPTY_VCF", "The VCF contains no variant records.");
  if (!sampleNames.length) {
    addWarning(warnings, null, "NO_SAMPLES", "This sites-only VCF has no genotype sample columns.");
  }
  if (!fileFormat) {
    addWarning(warnings, null, "MISSING_FILEFORMAT", "The ##fileformat declaration is missing.");
  }
  const genomeBuild = detectGenomeBuild(referenceValues, contigAssemblies);
  if (!genomeBuild) {
    addWarning(warnings, null, "GENOME_BUILD_NOT_DETECTED", "Genome build could not be detected from ##reference or contig assembly metadata.");
  }
  if (warnings.length === MAX_WARNINGS) {
    warnings.push({ line: null, code: "WARNINGS_TRUNCATED", message: "Additional warnings were omitted." });
  }

  return {
    metadata: {
      fileFormat,
      reference: referenceValues[0] ?? null,
      genomeBuild,
      sampleNames,
      selectedSample,
      totalRecords: dataLineCount,
      parsedRecords: records.length
    },
    records,
    warnings
  };
}
