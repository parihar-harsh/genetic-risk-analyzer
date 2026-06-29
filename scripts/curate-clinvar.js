import { createReadStream } from "node:fs";
import { mkdir, readFile, writeFile } from "node:fs/promises";
import { createGunzip } from "node:zlib";
import { createHash } from "node:crypto";
import path from "node:path";
import readline from "node:readline";
import { fileURLToPath } from "node:url";
import { locusKey, normalizeSmallVariant } from "../server/database.js";

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const argumentsList = process.argv.slice(2);
const option = (name) => {
  const index = argumentsList.indexOf(name);
  return index === -1 ? null : argumentsList[index + 1];
};
const inputPath = option("--input");
const apply = argumentsList.includes("--apply");
const trustedPanel = argumentsList.includes("--trusted-panel");
const expectedMd5 = option("--expected-md5");
const reportPath = option("--report") ?? path.join(root, "data", "clinvar-curation-report.json");

if (!inputPath) {
  console.error("Usage: node scripts/curate-clinvar.js --input /path/to/clinvar.vcf.gz [--apply]");
  process.exit(1);
}

const digest = createHash("md5");
for await (const chunk of createReadStream(inputPath)) digest.update(chunk);
const sourceMd5 = digest.digest("hex");
if (expectedMd5 && sourceMd5.toLowerCase() !== expectedMd5.toLowerCase()) {
  throw new Error(`ClinVar checksum mismatch: expected ${expectedMd5}, received ${sourceMd5}`);
}

const legacyVariants = JSON.parse(await readFile(path.join(root, "disease_data.json"), "utf8"));
const targetIndex = new Map();
for (const variant of legacyVariants) {
  const key = locusKey(variant.chromosome, variant.position, variant.ref, variant.alt);
  const targets = targetIndex.get(key) ?? [];
  targets.push(variant);
  targetIndex.set(key, targets);
}

const decodeClinVar = (value) => {
  if (!value) return null;
  try {
    return decodeURIComponent(value).replaceAll("_", " ");
  } catch {
    return value.replaceAll("_", " ");
  }
};

const parseInfo = (raw) => {
  const info = new Map();
  for (const field of raw.split(";")) {
    const separator = field.indexOf("=");
    if (separator === -1) info.set(field, true);
    else info.set(field.slice(0, separator), field.slice(separator + 1));
  }
  return info;
};

const reviewStars = (status = "") => {
  if (status.includes("no_assertion_criteria") || status.includes("no_classification")) return 0;
  if (status.includes("practice_guideline")) return 4;
  if (status.includes("reviewed_by_expert_panel")) return 3;
  if (status.includes("multiple_submitters") && status.includes("no_conflicts")) return 2;
  if (status.includes("criteria_provided")) return 1;
  return 0;
};

const trustedRules = [
  { label: "Hereditary Breast and Ovarian Cancer", condition: /Hereditary_breast_ovarian_cancer/i, genes: new Set(["BRCA1", "BRCA2", "PALB2"]) },
  { label: "Cystic Fibrosis", condition: /Cystic_fibrosis/i, genes: new Set(["CFTR"]) },
  { label: "Lynch Syndrome", condition: /Lynch_syndrome/i, genes: new Set(["MLH1", "MSH2", "MSH6", "PMS2", "EPCAM"]) },
  { label: "Familial Hypercholesterolemia", condition: /Familial_hypercholesterolemia/i, genes: new Set(["LDLR", "APOB", "PCSK9", "LDLRAP1"]) },
  { label: "Marfan Syndrome", condition: /Marfan_syndrome/i, genes: new Set(["FBN1"]) },
  { label: "Neurofibromatosis Type 1", condition: /Neurofibromatosis(?:,_|_)type_1/i, genes: new Set(["NF1"]) },
  { label: "Hereditary Hemochromatosis", condition: /Hemochromatosis/i, genes: new Set(["HFE", "HJV", "HAMP", "TFR2", "SLC40A1"]) },
  { label: "Alpha-1 Antitrypsin Deficiency", condition: /Alpha-1.antitrypsin_deficiency/i, genes: new Set(["SERPINA1"]) }
];
const acceptedClassifications = new Set([
  "Pathogenic",
  "Likely_pathogenic",
  "Pathogenic/Likely_pathogenic"
]);

const source = createReadStream(inputPath);
const lines = readline.createInterface({
  input: inputPath.endsWith(".gz") ? source.pipe(createGunzip()) : source,
  crlfDelay: Infinity
});

let sourceRelease = null;
let sourceBuild = null;
let sourceRecords = 0;
const matchedRecords = new Map();
const trustedRecords = new Map();

for await (const line of lines) {
  if (line.startsWith("##fileDate=")) sourceRelease = line.slice(11);
  if (line.startsWith("##reference=")) sourceBuild = line.slice(12);
  if (!line || line.startsWith("#")) continue;
  sourceRecords += 1;
  const [chromosome, positionRaw, variationId, ref, altRaw, , , infoRaw] = line.split("\t");
  const position = Number(positionRaw);
  const info = parseInfo(infoRaw);
  for (const alt of altRaw.split(",")) {
    const key = locusKey(chromosome, position, ref, alt);
    const normalized = normalizeSmallVariant(chromosome, position, ref, alt);
    const record = {
      key,
      chromosome: normalized.chromosome,
      position: normalized.position,
      ref: normalized.ref,
      alt: normalized.alt,
      variationId,
      alleleId: info.get("ALLELEID") ?? null,
      rsId: info.get("RS") ? `rs${info.get("RS")}` : null,
      classification: decodeClinVar(info.get("CLNSIG")),
      reviewStatus: decodeClinVar(info.get("CLNREVSTAT")),
      reviewStars: reviewStars(info.get("CLNREVSTAT")),
      conditions: decodeClinVar(info.get("CLNDN")),
      conditionIdentifiers: decodeClinVar(info.get("CLNDISDB")),
      hgvs: decodeClinVar(info.get("CLNHGVS")),
      genes: decodeClinVar(info.get("GENEINFO")),
      molecularConsequences: decodeClinVar(info.get("MC")),
      scvAccessions: decodeClinVar(info.get("CLNSIGSCV")),
      variantType: decodeClinVar(info.get("CLNVC"))
    };
    if (targetIndex.has(key)) {
      const existing = matchedRecords.get(key) ?? [];
      existing.push(record);
      matchedRecords.set(key, existing);
    }

    if (trustedPanel) {
      const geneSymbols = new Set((info.get("GENEINFO") ?? "")
        .split("|")
        .map((gene) => gene.split(":")[0]));
      const conditionsRaw = info.get("CLNDN") ?? "";
      const matchingRule = trustedRules.find((rule) =>
        rule.condition.test(conditionsRaw) &&
        [...geneSymbols].some((gene) => rule.genes.has(gene))
      );
      const isSequenceAllele = /^[ACGTN]+$/i.test(ref) && /^[ACGTN]+$/i.test(alt);
      if (
        matchingRule &&
        isSequenceAllele &&
        acceptedClassifications.has(info.get("CLNSIG")) &&
        record.reviewStars >= 2
      ) {
        record.matchedConditions = conditionsRaw
          .split("|")
          .filter((condition) => matchingRule.condition.test(condition))
          .map(decodeClinVar)
          .join("; ");
        record.conditionLabel = matchingRule.label;
        record.matchedGenes = [...geneSymbols].filter((gene) => matchingRule.genes.has(gene)).join("|");
        const existing = trustedRecords.get(key) ?? [];
        existing.push(record);
        trustedRecords.set(key, existing);
      }
    }
  }
}

const matches = [];
const unmatched = [];
for (const [key, variants] of targetIndex) {
  const records = matchedRecords.get(key);
  if (records?.length) matches.push({ key, legacyEntries: variants, clinvarRecords: records });
  else unmatched.push(...variants.map((variant) => ({
    key,
    disease: variant.disease,
    classification: variant.classification,
    description: variant.description
  })));
}

const recordsToCurate = trustedPanel
  ? [...trustedRecords.values()].flat()
  : matches.flatMap(({ clinvarRecords }) => clinvarRecords);
const curated = recordsToCurate
  .filter((record) => record.classification && record.conditions && record.reviewStatus)
  .map((record) => ({
    chromosome: record.chromosome,
    position: record.position,
    ref: record.ref,
    alt: record.alt,
    disease: record.conditionLabel ?? record.matchedConditions ?? record.conditions.replaceAll("|", "; "),
    description: `ClinVar aggregate germline classification for ${record.matchedConditions ?? record.conditions.replaceAll("|", "; ")}.`,
    isAgeSpecific: false,
    agePenetrance: [],
    alleleFrequency: null,
    gene: record.matchedGenes ?? record.genes,
    geneImpact: record.molecularConsequences ?? record.variantType ?? "Not provided",
    classification: record.classification,
    populationFrequencies: {},
    environmentalFactors: {},
    riskFactors: {},
    confidenceIntervals: [],
    oddsRatio: 1,
    evidence: {
      status: "source-matched",
      sourceName: "ClinVar",
      accession: `VariationID:${record.variationId}`,
      alleleId: record.alleleId,
      rsId: record.rsId,
      reviewStatus: record.reviewStatus,
      reviewStars: record.reviewStars,
      lastEvaluated: null,
      sourceRelease,
      sourceUrl: `https://www.ncbi.nlm.nih.gov/clinvar/variation/${record.variationId}/`,
      hgvs: record.hgvs,
      conditionIdentifiers: record.conditionIdentifiers,
      sourceConditions: record.conditions,
      matchedConditions: record.matchedConditions ?? null,
      sourceGenes: record.genes,
      scvAccessions: record.scvAccessions,
      effectEstimateValidated: false
    }
  }));
const countBy = (values) => Object.fromEntries([...values.reduce((counts, value) => {
  counts.set(value, (counts.get(value) ?? 0) + 1);
  return counts;
}, new Map()).entries()].sort(([left], [right]) => left.localeCompare(right)));

const report = {
  generatedAt: new Date().toISOString(),
  source: {
    name: "ClinVar",
    url: "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz",
    release: sourceRelease,
    genomeBuild: sourceBuild,
    md5: sourceMd5,
    recordsScanned: sourceRecords
  },
  inputVariants: legacyVariants.length,
  exactMatchedLoci: matches.length,
  trustedPanel,
  trustedPanelLoci: trustedRecords.size,
  curatedEntries: curated.length,
  unmatchedEntries: unmatched.length,
  panelSummary: {
    diseases: countBy(curated.map((entry) => entry.disease)),
    genes: countBy(curated.map((entry) => entry.gene)),
    classifications: countBy(curated.map((entry) => entry.classification)),
    reviewStars: countBy(curated.map((entry) => String(entry.evidence.reviewStars)))
  },
  legacyComparison: trustedPanel ? {
    exactMatchedLoci: matches.length,
    unmatchedEntries: unmatched.length
  } : {
    matches,
    unmatched
  }
};

await mkdir(path.dirname(reportPath), { recursive: true });
await writeFile(reportPath, `${JSON.stringify(report, null, 2)}\n`);

if (apply) {
  if (sourceBuild !== "GRCh38") throw new Error(`Expected GRCh38 ClinVar input, received ${sourceBuild}`);
  if (!sourceRelease) throw new Error("ClinVar source release date is missing");
  if (!curated.length) throw new Error("No exact ClinVar matches were found; refusing to replace the database");
  await writeFile(path.join(root, "disease_data.json"), `${JSON.stringify(curated, null, 2)}\n`);
  const metadataPath = path.join(root, "disease_data.metadata.json");
  const metadata = JSON.parse(await readFile(metadataPath, "utf8"));
  await writeFile(metadataPath, `${JSON.stringify({
    ...metadata,
    datasetVersion: `clinvar-${sourceRelease}`,
    releasedAt: sourceRelease,
    provenanceStatus: "source-matched",
    clinicalUseApproved: false,
    description: "Exact GRCh38 allele matches imported from the named ClinVar release. Classifications remain submitted assertions and are not independent clinical validation.",
    source: {
      name: "ClinVar",
      url: report.source.url,
      release: sourceRelease,
      md5: sourceMd5,
      recordsScanned: sourceRecords
    }
  }, null, 2)}\n`);
}

console.log(JSON.stringify({
  sourceRelease,
  sourceBuild,
  sourceRecords,
  sourceMd5,
  inputVariants: legacyVariants.length,
  exactMatchedLoci: matches.length,
  trustedPanel,
  trustedPanelLoci: trustedRecords.size,
  curatedEntries: curated.length,
  unmatchedEntries: unmatched.length,
  reportPath,
  applied: apply
}, null, 2));
