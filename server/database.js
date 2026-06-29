import { readFile } from "node:fs/promises";
import path from "node:path";
import { fileURLToPath } from "node:url";

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), "..");
const databasePath = path.join(root, "disease_data.json");
const metadataPath = path.join(root, "disease_data.metadata.json");
const supportedBuilds = new Set(["GRCh37", "GRCh38"]);

const requiredStrings = [
  "chromosome", "ref", "alt", "disease", "description", "geneImpact", "classification"
];

const assertObject = (value, label) => {
  if (!value || typeof value !== "object" || Array.isArray(value)) {
    throw new Error(`${label} must be an object`);
  }
};

function validateEvidence(evidence, index) {
  if (evidence === undefined) {
    return Object.freeze({
      status: "unverified",
      sourceName: null,
      accession: null,
      reviewStatus: "not provided",
      reviewStars: 0,
      lastEvaluated: null,
      effectEstimateValidated: false
    });
  }
  assertObject(evidence, `Variant ${index}.evidence`);
  const reviewStars = Number(evidence.reviewStars ?? 0);
  if (!Number.isInteger(reviewStars) || reviewStars < 0 || reviewStars > 4) {
    throw new Error(`Variant ${index}.evidence.reviewStars must be from 0 through 4`);
  }
  const allowedStatuses = new Set(["unverified", "source-matched", "verified"]);
  return Object.freeze({
    status: allowedStatuses.has(evidence.status) ? evidence.status : "unverified",
    sourceName: typeof evidence.sourceName === "string" ? evidence.sourceName : null,
    accession: typeof evidence.accession === "string" ? evidence.accession : null,
    alleleId: typeof evidence.alleleId === "string" ? evidence.alleleId : null,
    rsId: typeof evidence.rsId === "string" ? evidence.rsId : null,
    reviewStatus: typeof evidence.reviewStatus === "string" ? evidence.reviewStatus : "not provided",
    reviewStars,
    lastEvaluated: typeof evidence.lastEvaluated === "string" ? evidence.lastEvaluated : null,
    sourceRelease: typeof evidence.sourceRelease === "string" ? evidence.sourceRelease : null,
    sourceUrl: typeof evidence.sourceUrl === "string" ? evidence.sourceUrl : null,
    hgvs: typeof evidence.hgvs === "string" ? evidence.hgvs : null,
    sourceConditions: typeof evidence.sourceConditions === "string" ? evidence.sourceConditions : null,
    matchedConditions: typeof evidence.matchedConditions === "string" ? evidence.matchedConditions : null,
    sourceGenes: typeof evidence.sourceGenes === "string" ? evidence.sourceGenes : null,
    conditionIdentifiers: typeof evidence.conditionIdentifiers === "string"
      ? evidence.conditionIdentifiers
      : null,
    scvAccessions: typeof evidence.scvAccessions === "string" ? evidence.scvAccessions : null,
    effectEstimateValidated: evidence.effectEstimateValidated === true
  });
}

function validateVariant(entry, index, metadata) {
  assertObject(entry, `Variant ${index}`);
  for (const key of requiredStrings) {
    if (typeof entry[key] !== "string" || !entry[key].trim()) {
      throw new Error(`Variant ${index}.${key} must be a non-empty string`);
    }
  }
  if (!Number.isInteger(entry.position) || entry.position < 1) {
    throw new Error(`Variant ${index}.position must be a positive integer`);
  }
  if (typeof entry.isAgeSpecific !== "boolean") {
    throw new Error(`Variant ${index}.isAgeSpecific must be boolean`);
  }
  if (!Array.isArray(entry.agePenetrance)) {
    throw new Error(`Variant ${index}.agePenetrance must be an array`);
  }
  for (const [ageIndex, point] of entry.agePenetrance.entries()) {
    assertObject(point, `Variant ${index}.agePenetrance[${ageIndex}]`);
    if (!Number.isInteger(point.age) || point.age < 0 || point.age > 120) {
      throw new Error(`Variant ${index} contains an invalid penetrance age`);
    }
    if (!Number.isFinite(point.penetrance) || point.penetrance < 0 || point.penetrance > 1) {
      throw new Error(`Variant ${index} contains an invalid penetrance value`);
    }
  }
  if (!Number.isFinite(entry.oddsRatio) || entry.oddsRatio <= 0) {
    throw new Error(`Variant ${index}.oddsRatio must be positive`);
  }
  for (const key of ["populationFrequencies", "environmentalFactors", "riskFactors"]) {
    assertObject(entry[key] ?? {}, `Variant ${index}.${key}`);
  }
  return Object.freeze({
    ...entry,
    chromosome: normalizeChromosome(entry.chromosome),
    ref: entry.ref.toUpperCase(),
    alt: entry.alt.toUpperCase(),
    gene: typeof entry.gene === "string" ? entry.gene : null,
    genomeBuild: entry.genomeBuild ?? metadata.genomeBuild,
    evidence: validateEvidence(entry.evidence, index),
    agePenetrance: Object.freeze([...entry.agePenetrance].sort((a, b) => a.age - b.age))
  });
}

export function normalizeChromosome(value) {
  const normalized = String(value).trim().replace(/^chr/i, "").toUpperCase();
  return normalized === "M" ? "MT" : normalized;
}

export function normalizeSmallVariant(chromosome, position, ref, alt) {
  let normalizedPosition = Number(position);
  let normalizedRef = ref.toUpperCase();
  let normalizedAlt = alt.toUpperCase();

  while (
    normalizedRef.length > 1 &&
    normalizedAlt.length > 1 &&
    normalizedRef.at(-1) === normalizedAlt.at(-1)
  ) {
    normalizedRef = normalizedRef.slice(0, -1);
    normalizedAlt = normalizedAlt.slice(0, -1);
  }
  while (
    normalizedRef.length > 1 &&
    normalizedAlt.length > 1 &&
    normalizedRef[0] === normalizedAlt[0]
  ) {
    normalizedRef = normalizedRef.slice(1);
    normalizedAlt = normalizedAlt.slice(1);
    normalizedPosition += 1;
  }
  return {
    chromosome: normalizeChromosome(chromosome),
    position: normalizedPosition,
    ref: normalizedRef,
    alt: normalizedAlt
  };
}

export function locusKey(chromosome, position, ref, alt) {
  const normalized = normalizeSmallVariant(chromosome, position, ref, alt);
  return `${normalized.chromosome}:${normalized.position}:${normalized.ref}:${normalized.alt}`;
}

function validateMetadata(metadata) {
  assertObject(metadata, "Database metadata");
  if (!Number.isInteger(metadata.schemaVersion) || metadata.schemaVersion < 1) {
    throw new Error("Database metadata schemaVersion must be a positive integer");
  }
  if (typeof metadata.datasetVersion !== "string" || !metadata.datasetVersion) {
    throw new Error("Database metadata datasetVersion is required");
  }
  if (!supportedBuilds.has(metadata.genomeBuild)) {
    throw new Error("Database metadata genomeBuild must be GRCh37 or GRCh38");
  }
  if (typeof metadata.clinicalUseApproved !== "boolean") {
    throw new Error("Database metadata clinicalUseApproved must be boolean");
  }
  return Object.freeze(metadata);
}

export async function loadDatabase() {
  const [databaseText, metadataText] = await Promise.all([
    readFile(databasePath, "utf8"),
    readFile(metadataPath, "utf8")
  ]);
  const parsed = JSON.parse(databaseText);
  const metadata = validateMetadata(JSON.parse(metadataText));
  if (!Array.isArray(parsed) || parsed.length === 0) {
    throw new Error("Disease database must be a non-empty array");
  }

  const variants = parsed.map((entry, index) => validateVariant(entry, index, metadata));
  const index = new Map();
  for (const variant of variants) {
    const key = locusKey(variant.chromosome, variant.position, variant.ref, variant.alt);
    const existing = index.get(key) ?? [];
    existing.push(variant);
    index.set(key, existing);
  }

  return Object.freeze({
    metadata,
    variants: Object.freeze(variants),
    index,
    diseases: Object.freeze([...new Set(variants.map((item) => item.disease))].sort())
  });
}
