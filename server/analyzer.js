import { locusKey } from "./database.js";
import { badRequest } from "./errors.js";

const pathogenic = new Set([
  "Pathogenic",
  "Likely Pathogenic",
  "Likely pathogenic",
  "Pathogenic/Likely pathogenic"
]);
const riskClassifications = new Set(["Risk Allele"]);
const informational = new Set(["Benign", "Likely Benign"]);

function zygosity(record, alleleNumber, dosage) {
  const calledAlleles = record.alleleIndexes.length;
  if (calledAlleles === 1) return dosage ? "hemizygous/haploid alternate" : "haploid reference";
  if (dosage === calledAlleles) return "homozygous alternate";
  return "heterozygous";
}

function callQuality(record) {
  const issues = [];
  if (record.filter && record.filter !== "PASS" && record.filter !== ".") issues.push(`FILTER=${record.filter}`);
  if (record.sampleFilter && record.sampleFilter !== "PASS" && record.sampleFilter !== ".") {
    issues.push(`sample FT=${record.sampleFilter}`);
  }
  if (record.quality !== null && record.quality < 20) issues.push(`low site quality (${record.quality})`);
  if (record.genotypeQuality !== null && record.genotypeQuality < 20) issues.push(`low genotype quality (${record.genotypeQuality})`);
  if (record.depth !== null && record.depth < 10) issues.push(`low read depth (${record.depth})`);
  if (record.depth === null) issues.push("read depth unavailable");
  return { status: issues.length ? "review" : "pass", issues };
}

function ageAnnotation(variant, age) {
  if (!variant.isAgeSpecific || !variant.agePenetrance.length) return null;
  const eligible = variant.agePenetrance.filter((point) => point.age <= age);
  if (!eligible.length) {
    return { status: "below-first-threshold", firstReferenceAge: variant.agePenetrance[0].age };
  }
  const point = eligible.at(-1);
  return {
    status: "reference-only",
    age: point.age,
    databasePenetrance: point.penetrance,
    note: "Unvalidated database annotation; not a personalized probability."
  };
}

function priorityFor(classification, quality, evidence) {
  if (quality.status !== "pass") return "review";
  if (evidence.status !== "verified" && evidence.status !== "source-matched") return "review";
  const conflicting = /conflict/i.test(evidence.reviewStatus) && !/no conflicts?/i.test(evidence.reviewStatus);
  if (evidence.reviewStars < 2 || conflicting) return "review";
  if (pathogenic.has(classification)) return "high";
  if (riskClassifications.has(classification)) return "moderate";
  return "informational";
}

export function validateAnalysisInput(body = {}) {
  const ageRaw = body.age;
  const age = ageRaw === undefined || ageRaw === "" ? null : Number(ageRaw);
  if (age !== null && (!Number.isInteger(age) || age < 1 || age > 120)) {
    throw badRequest("INVALID_AGE", "Age must be a whole number between 1 and 120.");
  }

  let factors = [];
  if (body.factors !== undefined) {
    try {
      factors = typeof body.factors === "string" ? JSON.parse(body.factors) : body.factors;
    } catch {
      throw badRequest("INVALID_FACTORS", "Factors must be a JSON array of strings.");
    }
    if (!Array.isArray(factors) || factors.some((item) => typeof item !== "string")) {
      throw badRequest("INVALID_FACTORS", "Factors must be a JSON array of strings.");
    }
    if (factors.some((item) => item.length > 100)) {
      throw badRequest("INVALID_FACTORS", "Factor names cannot exceed 100 characters.");
    }
  }
  if (typeof body.sampleName === "string" && body.sampleName.length > 200) {
    throw badRequest("INVALID_SAMPLE", "Sample names cannot exceed 200 characters.");
  }
  return {
    age,
    sampleName: typeof body.sampleName === "string" && body.sampleName.trim() ? body.sampleName.trim() : null,
    genomeBuild: typeof body.genomeBuild === "string" && body.genomeBuild.trim()
      ? body.genomeBuild.trim()
      : null,
    factors: [...new Set(factors.map((item) => item.trim()).filter(Boolean))].slice(0, 20)
  };
}

export function resolveGenomeBuild(parsed, database, requestedBuild) {
  const supported = new Set(["GRCh37", "GRCh38"]);
  if (requestedBuild && !supported.has(requestedBuild)) {
    throw badRequest("INVALID_GENOME_BUILD", "Genome build must be GRCh37 or GRCh38.");
  }
  if (requestedBuild && parsed.metadata.genomeBuild && requestedBuild !== parsed.metadata.genomeBuild) {
    throw badRequest(
      "GENOME_BUILD_CONFLICT",
      `Selected build ${requestedBuild} conflicts with VCF metadata (${parsed.metadata.genomeBuild}).`
    );
  }
  const effectiveBuild = parsed.metadata.genomeBuild ?? requestedBuild;
  if (!effectiveBuild) {
    throw badRequest(
      "GENOME_BUILD_REQUIRED",
      "Genome build could not be detected. Select GRCh37 or GRCh38 explicitly."
    );
  }
  if (effectiveBuild !== database.metadata.genomeBuild) {
    throw badRequest(
      "GENOME_BUILD_MISMATCH",
      `This database uses ${database.metadata.genomeBuild}; the uploaded VCF uses ${effectiveBuild}. Liftover is required before analysis.`
    );
  }
  return effectiveBuild;
}

export function analyzeVcf(parsed, database, input) {
  const findings = [];
  const scoreByDisease = new Map();
  const seenFindings = new Set();
  const analysisWarnings = [];

  for (const record of parsed.records) {
    if (record.filter && record.filter !== "PASS" && record.filter !== ".") {
      analysisWarnings.push({
        line: record.line,
        code: "FILTERED_CALL_EXCLUDED",
        message: `Record was excluded because FILTER=${record.filter}.`
      });
      continue;
    }
    if (record.sampleFilter && record.sampleFilter !== "PASS" && record.sampleFilter !== ".") {
      analysisWarnings.push({
        line: record.line,
        code: "SAMPLE_FILTERED_CALL_EXCLUDED",
        message: `Record was excluded because sample FT=${record.sampleFilter}.`
      });
      continue;
    }
    for (let altIndex = 0; altIndex < record.alts.length; altIndex += 1) {
      if (!/^[ACGTN]+$/.test(record.alts[altIndex])) continue;
      const alleleNumber = altIndex + 1;
      const dosage = record.alleleIndexes.filter((allele) => allele === alleleNumber).length;
      if (!dosage) continue;
      const matches = database.index.get(
        locusKey(record.chromosome, record.position, record.ref, record.alts[altIndex])
      ) ?? [];

      for (const variant of matches) {
        const findingKey = `${locusKey(record.chromosome, record.position, record.ref, record.alts[altIndex])}:${variant.disease}`;
        if (seenFindings.has(findingKey)) {
          analysisWarnings.push({
            line: record.line,
            code: "DUPLICATE_MATCH",
            message: "A duplicate database match was skipped to avoid double-counting."
          });
          continue;
        }
        seenFindings.add(findingKey);
        const quality = callQuality(record);
        const activeModifiers = input.factors.flatMap((factor) => {
          const multiplier = variant.environmentalFactors?.[factor] ?? variant.riskFactors?.[factor];
          return multiplier === undefined ? [] : [{ factor, databaseMultiplier: multiplier }];
        });
        const logOddsContribution = Math.log(variant.oddsRatio) * dosage;
        if (
          variant.evidence.status === "verified" &&
          variant.evidence.effectEstimateValidated &&
          Number.isFinite(logOddsContribution) &&
          !informational.has(variant.classification)
        ) {
          const score = scoreByDisease.get(variant.disease) ?? { logOddsScore: 0, contributors: 0 };
          score.logOddsScore += logOddsContribution;
          score.contributors += 1;
          scoreByDisease.set(variant.disease, score);
        }

        findings.push({
          locus: `${record.chromosome}:${record.position}`,
          id: record.id,
          ref: record.ref,
          alt: record.alts[altIndex],
          genotype: record.genotype,
          dosage,
          zygosity: zygosity(record, alleleNumber, dosage),
          disease: variant.disease,
          description: variant.description,
          gene: variant.gene,
          classification: variant.classification,
          geneImpact: variant.geneImpact,
          alleleFrequency: variant.alleleFrequency,
          populationFrequencies: variant.populationFrequencies,
          oddsRatio: variant.oddsRatio,
          callQuality: quality,
          priority: priorityFor(variant.classification, quality, variant.evidence),
          evidence: variant.evidence,
          ageAnnotation: input.age ? ageAnnotation(variant, input.age) : null,
          activeModifiers,
          sourceLine: record.line
        });
      }
    }
  }

  const order = { high: 0, review: 1, moderate: 2, informational: 3 };
  findings.sort((a, b) => order[a.priority] - order[b.priority] || a.disease.localeCompare(b.disease));
  const relativeScores = [...scoreByDisease.entries()]
    .map(([disease, score]) => {
      const multiplier = Math.exp(score.logOddsScore);
      return {
        disease,
        ...score,
        relativeOddsMultiplier: Number.isFinite(multiplier) ? Number(multiplier.toPrecision(6)) : null,
        interpretation: "Relative index only; not an absolute probability or validated PRS."
      };
    })
    .sort((a, b) => b.logOddsScore - a.logOddsScore);

  return {
    generatedAt: new Date().toISOString(),
    disclaimer: database.metadata.clinicalUseApproved
      ? "This software does not provide a diagnosis; findings require professional review."
      : "For research and education only. The bundled evidence is unverified and must not be used for clinical decisions.",
    database: database.metadata,
    input: {
      ...parsed.metadata,
      genomeBuild: input.genomeBuild ?? parsed.metadata.genomeBuild,
      age: input.age,
      activeFactors: input.factors
    },
    summary: {
      databaseVariants: database.variants.length,
      recordsProcessed: parsed.metadata.parsedRecords,
      exactAlternateMatches: findings.length,
      highPriority: findings.filter((item) => item.priority === "high").length,
      reviewRequired: findings.filter((item) => item.priority === "review").length
    },
    findings,
    relativeScores,
    warnings: [...parsed.warnings, ...analysisWarnings].slice(0, 101)
  };
}
