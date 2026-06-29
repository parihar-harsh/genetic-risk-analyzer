import { loadDatabase } from "../server/database.js";

const database = await loadDatabase();
const sourceMatched = database.variants.filter((variant) =>
  variant.evidence.status === "source-matched" || variant.evidence.status === "verified"
);
const independentlyVerified = database.variants.filter((variant) => variant.evidence.status === "verified");
const consensusOrExpertReviewed = sourceMatched.filter((variant) => variant.evidence.reviewStars >= 2);
const effectValidated = sourceMatched.filter((variant) => variant.evidence.effectEstimateValidated);
const withAccessions = database.variants.filter((variant) => variant.evidence.accession);
const withEvaluationDates = database.variants.filter((variant) => variant.evidence.lastEvaluated);
const withSourceRelease = database.variants.filter((variant) => variant.evidence.sourceRelease);
const counts = {
  datasetVersion: database.metadata.datasetVersion,
  genomeBuild: database.metadata.genomeBuild,
  totalVariants: database.variants.length,
  sourceMatchedEvidence: sourceMatched.length,
  independentlyVerifiedEvidence: independentlyVerified.length,
  consensusOrExpertReviewed: consensusOrExpertReviewed.length,
  sourceAccessions: withAccessions.length,
  evaluationDates: withEvaluationDates.length,
  sourceReleaseDates: withSourceRelease.length,
  validatedEffectEstimates: effectValidated.length,
  clinicalUseApproved: database.metadata.clinicalUseApproved
};

console.log(JSON.stringify(counts, null, 2));

const incomplete =
  sourceMatched.length !== database.variants.length ||
  withAccessions.length !== database.variants.length ||
  withSourceRelease.length !== database.variants.length;

if (process.argv.includes("--strict") && incomplete) {
  console.error("Dataset evidence audit failed. See MANUAL_ACTIONS.md.");
  process.exitCode = 1;
}
