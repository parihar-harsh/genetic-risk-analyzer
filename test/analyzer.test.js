import assert from "node:assert/strict";
import { readFile } from "node:fs/promises";
import { gzipSync } from "node:zlib";
import test from "node:test";
import { analyzeVcf, resolveGenomeBuild, validateAnalysisInput } from "../server/analyzer.js";
import { loadDatabase, locusKey } from "../server/database.js";
import { parseVcfUpload } from "../server/vcf.js";

const database = await loadDatabase();
const knownVariant = database.variants.find((variant) => variant.ref.length === 1 && variant.alt.length === 1);
const knownRecord = ({ genotype = "0/1", quality = 90, filter = "PASS", format = "GT:GQ:DP", sample = null } = {}) =>
  `${knownVariant.chromosome}\t${knownVariant.position}\t.\t${knownVariant.ref}\t${knownVariant.alt}\t${quality}\t${filter}\t.\t${format}\t${sample ?? `${genotype}:90:30`}`;
const upload = (text, name = "test.vcf") => ({
  originalname: name,
  buffer: Buffer.from(text)
});
const vcf = (records, samples = ["SAMPLE"]) => [
  "##fileformat=VCFv4.2",
  "##reference=GRCh38",
  `#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${samples.join("\t")}`,
  ...records
].join("\n");

test("database loads and indexes the expected variants", () => {
  assert.equal(database.variants.length, 11511);
  assert.equal(database.diseases.length, 8);
  assert.ok(database.index.size > 11_000);
  assert.equal(database.metadata.genomeBuild, "GRCh38");
  assert.equal(database.metadata.provenanceStatus, "source-matched");
  assert.ok(database.variants.every((variant) => variant.evidence.status === "source-matched"));
  assert.ok(database.variants.every((variant) => variant.evidence.reviewStars >= 2));
  assert.ok(database.variants.every((variant) => variant.evidence.effectEstimateValidated === false));
  assert.ok(database.variants.every((variant) =>
    ["Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"].includes(variant.classification)
  ));
});

test("bundled synthetic fixtures resolve only sourced passing alternate calls", async () => {
  const sampleBuffer = await readFile(new URL("../sample_input.vcf", import.meta.url));
  const sampleParsed = parseVcfUpload({ originalname: "sample_input.vcf", buffer: sampleBuffer });
  const sampleReport = analyzeVcf(sampleParsed, database, { age: null, factors: [], sampleName: null });
  assert.equal(sampleReport.findings.length, 3);
  assert.ok(sampleReport.findings.every((finding) => finding.priority === "high"));

  const patientBuffer = await readFile(new URL("../patient001.vcf", import.meta.url));
  const patientParsed = parseVcfUpload({ originalname: "patient001.vcf", buffer: patientBuffer });
  const patientReport = analyzeVcf(patientParsed, database, { age: null, factors: [], sampleName: null });
  assert.equal(patientReport.findings.length, 5);
  assert.equal(patientReport.warnings.at(-1).code, "FILTERED_CALL_EXCLUDED");
});

test("VCF FORMAT fields use field indexes rather than character offsets", () => {
  const parsed = parseVcfUpload(upload(vcf([
    "7\t117559594\trs113993960\tCTT\tC\t1500\tPASS\tDP=110\tGT:GQ:DP\t1/1:99:110"
  ])));
  assert.equal(parsed.records[0].genotype, "1/1");
  assert.equal(parsed.records[0].genotypeQuality, 99);
  assert.equal(parsed.records[0].depth, 110);
});

test("exact allele analysis returns annotations without fabricated probabilities", () => {
  const parsed = parseVcfUpload(upload(vcf([knownRecord()])));
  const report = analyzeVcf(parsed, database, { age: 50, factors: [], sampleName: null });
  assert.equal(report.findings.length, 1);
  assert.equal(report.findings[0].disease, knownVariant.disease);
  assert.equal(report.findings[0].callQuality.status, "pass");
  assert.equal(report.findings[0].priority, "high");
  assert.equal(report.findings[0].evidence.status, "source-matched");
  assert.equal(report.findings[0].ageAnnotation, null);
  assert.ok(!("risk" in report.findings[0]));
});

test("chromosome aliases, phased genotypes, multi-allelic ALT, and sample selection work", () => {
  const unusedAlt = ["A", "C", "G", "T"].find((allele) =>
    allele !== knownVariant.ref && allele !== knownVariant.alt
  );
  const text = vcf([
    `chr${knownVariant.chromosome}\t${knownVariant.position}\t.\t${knownVariant.ref}\t${unusedAlt},${knownVariant.alt}\t95\tPASS\t.\tGQ:GT:DP\t99:0/0:40\t91:0|2:35`
  ], ["CONTROL", "CASE"]);
  const parsed = parseVcfUpload(upload(text), "CASE");
  const report = analyzeVcf(parsed, database, { age: null, factors: [], sampleName: "CASE" });
  assert.equal(parsed.metadata.selectedSample, "CASE");
  assert.equal(report.findings.length, 1);
  assert.equal(report.findings[0].alt, knownVariant.alt);
  assert.equal(report.findings[0].genotype, "0|2");
  assert.equal(report.findings[0].zygosity, "heterozygous");
});

test("a reference genotype creates no finding or baseline score", () => {
  const parsed = parseVcfUpload(upload(vcf([
    knownRecord({ format: "GT:DP", sample: "0/0:35" })
  ])));
  const report = analyzeVcf(parsed, database, { age: null, factors: [], sampleName: null });
  assert.deepEqual(report.findings, []);
  assert.deepEqual(report.relativeScores, []);
});

test("unmodeled user factors cannot change sourced classifications", () => {
  const parsed = parseVcfUpload(upload(vcf([knownRecord()])));
  const inactive = analyzeVcf(parsed, database, { age: null, factors: [], sampleName: null });
  const active = analyzeVcf(parsed, database, { age: null, factors: ["Smoking"], sampleName: null });
  assert.deepEqual(inactive.findings[0].activeModifiers, []);
  assert.deepEqual(active.findings[0].activeModifiers, []);
  assert.equal(active.findings[0].oddsRatio, inactive.findings[0].oddsRatio);
  assert.ok(!("risk" in active.findings[0]));
});

test("failed filters and low call metrics are flagged for review", () => {
  const parsed = parseVcfUpload(upload(vcf([
    knownRecord({ quality: 18, sample: "0/1:15:9" })
  ])));
  const report = analyzeVcf(parsed, database, { age: null, factors: [], sampleName: null });
  assert.equal(report.findings[0].priority, "review");
  assert.equal(report.findings[0].callQuality.issues.length, 3);
});

test("failed site and sample filters are excluded rather than interpreted", () => {
  const parsed = parseVcfUpload(upload(vcf([
    knownRecord({ filter: "q10", format: "GT:DP", sample: "0/1:30" }),
    knownRecord({ format: "GT:DP:FT", sample: "0/1:30:failed" })
  ])));
  const report = analyzeVcf(parsed, database, { age: null, factors: [], sampleName: null });
  assert.equal(report.findings.length, 0);
  assert.deepEqual(
    report.warnings.map((warning) => warning.code),
    ["FILTERED_CALL_EXCLUDED", "SAMPLE_FILTERED_CALL_EXCLUDED"]
  );
});

test("gzip uploads are decoded", () => {
  const text = vcf(["7\t117559594\t.\tCTT\tC\t.\tPASS\t.\tGT\t0/1"]);
  const parsed = parseVcfUpload({ originalname: "test.vcf.gz", buffer: gzipSync(text) });
  assert.equal(parsed.records.length, 1);
});

test("invalid records generate bounded warnings instead of crashing", () => {
  const parsed = parseVcfUpload(upload(vcf([
    "10\tbad\t.\tC\tT\t.\tPASS\t.\tGT\t0/1",
    "10\t114758349\t.\tC\tT\tnot-a-number\tPASS\t.\tGT\t0/1"
  ])));
  assert.equal(parsed.records.length, 1);
  assert.deepEqual(parsed.warnings.map((item) => item.code), ["INVALID_POSITION", "INVALID_QUALITY"]);
});

test("duplicate VCF records are not double-counted", () => {
  const record = knownRecord({ format: "GT:DP", sample: "0/1:30" });
  const parsed = parseVcfUpload(upload(vcf([record, record])));
  const report = analyzeVcf(parsed, database, { age: null, factors: [], sampleName: null });
  assert.equal(report.findings.length, 1);
  assert.deepEqual(report.relativeScores, []);
  assert.equal(report.warnings.at(-1).code, "DUPLICATE_MATCH");
});

test("minimal allele normalization matches equivalent padded representations", () => {
  assert.equal(locusKey("chr1", 100, "AAC", "AGC"), locusKey("1", 101, "A", "G"));
});

test("genome build is required and must match the database", () => {
  const parsed = parseVcfUpload(upload(vcf([knownRecord({ format: "GT", sample: "0/1" })])));
  assert.equal(resolveGenomeBuild(parsed, database, null), "GRCh38");
  assert.throws(
    () => resolveGenomeBuild({ metadata: { genomeBuild: "GRCh37" } }, database, null),
    /Liftover is required/
  );
  assert.throws(
    () => resolveGenomeBuild({ metadata: { genomeBuild: null } }, database, null),
    /Genome build could not be detected/
  );
});

test("analysis input validation rejects invalid age and factors", () => {
  assert.throws(() => validateAnalysisInput({ age: "0" }), /Age must/);
  assert.throws(() => validateAnalysisInput({ factors: "not-json" }), /Factors must/);
  assert.deepEqual(validateAnalysisInput({ age: "42", factors: '["Smoking","Smoking"]' }), {
    age: 42,
    sampleName: null,
    genomeBuild: null,
    factors: ["Smoking"]
  });
});
