import assert from "node:assert/strict";
import { once } from "node:events";
import test, { after, before } from "node:test";
import { createApp } from "../server/app.js";
import { loadDatabase } from "../server/database.js";

let server;
let baseUrl;
let database;

before(async () => {
  database = await loadDatabase();
  const app = createApp(database);
  server = app.listen(0, "127.0.0.1");
  await once(server, "listening");
  baseUrl = `http://127.0.0.1:${server.address().port}`;
});

after(async () => {
  await new Promise((resolve, reject) => server.close((error) => error ? reject(error) : resolve()));
});

test("health endpoint reports readiness", async () => {
  const response = await fetch(`${baseUrl}/api/v1/health`);
  const body = await response.json();
  assert.equal(response.status, 200);
  assert.equal(body.status, "ok");
  assert.equal(body.databaseVariants, 11511);
  assert.ok(response.headers.get("x-request-id"));
  assert.equal(response.headers.get("x-content-type-options"), "nosniff");
});

test("analysis endpoint processes multipart VCF and disables caching", async () => {
  const variant = database.variants.find((item) => item.ref.length === 1 && item.alt.length === 1);
  const content = [
    "##fileformat=VCFv4.2",
    "##reference=GRCh38",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    `${variant.chromosome}\t${variant.position}\t.\t${variant.ref}\t${variant.alt}\t95\tPASS\t.\tGT:GQ:DP\t0/1:90:40`
  ].join("\n");
  const form = new FormData();
  form.append("vcf", new Blob([content], { type: "text/plain" }), "sample.vcf");
  form.append("age", "50");
  form.append("factors", "[]");
  const response = await fetch(`${baseUrl}/api/v1/analyze`, { method: "POST", body: form });
  const body = await response.json();
  assert.equal(response.status, 200);
  assert.equal(response.headers.get("cache-control"), "no-store");
  assert.equal(body.findings.length, 1);
});

test("analysis endpoint returns structured errors", async () => {
  const response = await fetch(`${baseUrl}/api/v1/analyze`, { method: "POST", body: new FormData() });
  const body = await response.json();
  assert.equal(response.status, 400);
  assert.equal(body.error.code, "VCF_REQUIRED");
  assert.ok(body.error.requestId);
});

test("analysis rejects unknown or incompatible genome builds", async () => {
  const content = [
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    "10\t114758349\t.\tC\tT\t95\tPASS\t.\tGT\t0/1"
  ].join("\n");
  const missingBuild = new FormData();
  missingBuild.append("vcf", new Blob([content]), "sample.vcf");
  const missingResponse = await fetch(`${baseUrl}/api/v1/analyze`, { method: "POST", body: missingBuild });
  assert.equal(missingResponse.status, 400);
  assert.equal((await missingResponse.json()).error.code, "GENOME_BUILD_REQUIRED");

  const mismatch = new FormData();
  mismatch.append("vcf", new Blob([content]), "sample.vcf");
  mismatch.append("genomeBuild", "GRCh37");
  const mismatchResponse = await fetch(`${baseUrl}/api/v1/analyze`, { method: "POST", body: mismatch });
  assert.equal(mismatchResponse.status, 400);
  assert.equal((await mismatchResponse.json()).error.code, "GENOME_BUILD_MISMATCH");
});

test("unknown API routes return JSON 404 rather than the SPA", async () => {
  const response = await fetch(`${baseUrl}/api/v1/unknown`);
  assert.equal(response.status, 404);
  assert.equal((await response.json()).error.code, "API_NOT_FOUND");
});

test("SPA fallback serves the application", async () => {
  const response = await fetch(`${baseUrl}/any/client/route`);
  assert.equal(response.status, 200);
  assert.match(await response.text(), /GeneCheck/);
});
