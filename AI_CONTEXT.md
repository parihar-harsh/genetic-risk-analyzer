# AI Project Context and Progress

Last updated: 2026-06-29  
Project: `genetic-risk-analyzer`  
Current application version: 2.2.0  
Branch observed during migration: `master`

## Objective

Turn the original educational C++ genetic-risk CLI into a deployment-ready web application with a usable frontend, robust backend edge-case handling, explicit medical-safety boundaries, automated tests, and reproducible deployment.

## Current status

Implementation is complete for the version-2 scope:

- [x] Node.js backend and versioned API
- [x] Responsive frontend
- [x] Exact allele matching
- [x] Defensive VCF/VCF.GZ parsing
- [x] Multi-sample and multi-allelic handling
- [x] Quality annotations
- [x] Removal of fabricated absolute-risk output
- [x] Security and privacy middleware
- [x] Unit and integration tests
- [x] Docker production definition
- [x] README and environment documentation
- [x] Dependency audit
- [x] Reproducible ClinVar data curation

Latest verification:

- JavaScript syntax checks: passed
- Automated tests: 21 passed, 0 failed
- Dependency audit: 0 known vulnerabilities
- Curated dataset audit: 11,511/11,511 records source-matched; strict audit passed
- Docker build: not executed because no Docker daemon was running in the development environment
- Visual browser automation: not executed because the browser tool was unavailable in the development environment

Run `npm run check` and `npm audit` before trusting this status after future changes.

## Architecture

The application is a single deployable Node.js service:

```text
Browser UI (public/)
       |
       | multipart POST /api/v1/analyze
       v
Express app (server/app.js)
       |
       +-- VCF parser (server/vcf.js)
       +-- validated/indexed JSON database (server/database.js)
       +-- exact-match analysis (server/analyzer.js)
```

There is no database server and no persistent upload storage. The service holds one bounded upload in memory, produces JSON, and releases it after the request.

## Important design decisions

1. **The C++ code is legacy, not the production entrypoint.** It is retained for reference only.
2. **FASTA/manual sequence input was removed.** Valid variant calling requires alignment and genome-build context.
3. **No absolute disease probability is calculated.** The old algorithm returned 50% with no variants and could exceed 100%.
4. **Risk scores are disabled.** ClinVar classifications are not a validated patient-level risk model.
5. **Matching is exact:** normalized chromosome + position + REF + specific ALT.
6. **Lifestyle-factor calculations were removed.** The curated source does not support them.
7. **Uploads are not persisted or logged.** Responses are `no-store`.
8. **The JSON database is generated.** It is pinned to ClinVar GRCh38 2026-06-27 and every record carries source/review metadata.

Version 2.1 adds dataset-level metadata and accuracy gates:

- The annotation database is explicitly GRCh38.
- VCF build is detected from `##reference`/contig assembly metadata or must be selected.
- GRCh37 input is rejected rather than silently matched to GRCh38 coordinates.
- Small sequence alleles are reduced to minimal representation before matching.
- Failed site/sample filters are excluded.
- Every legacy assertion defaults to unverified evidence.
- Unverified assertions cannot be high-priority or contribute effect estimates.
- `npm run audit:data` reports evidence-curation completeness.

Version 2.2 replaces the legacy 94-record dataset:

- Scanned 4,439,382 official ClinVar GRCh38 records.
- Legacy audit found only two exact alleles, both contradicting the project claim.
- Generated 11,511 exact sequence alleles across eight inherited disease groups.
- Restricted classifications to pathogenic/likely pathogenic.
- Required two-star consensus, expert panel, or practice guideline status.
- Restricted each condition to configured established genes.
- Stored Variation ID, Allele ID, SCVs, review state, release, HGVS, identifiers, URL, and source condition.
- Kept effect estimates disabled and `clinicalUseApproved=false`.

`MANUAL_ACTIONS.md` describes remaining independent scientific and deployment review.

Do not weaken these safety boundaries without a validated scientific model and explicit product requirements.

## Backend behavior

`server/vcf.js` handles:

- UTF-8 BOM and CRLF
- Plain and gzip VCF
- Missing `##fileformat`
- Missing/invalid records as bounded warnings
- Sites-only VCFs
- Multiple samples with exact sample selection
- FORMAT fields in arbitrary order
- `GT`, `GQ`, `DP`, and INFO/DP fallback
- `/` and `|` genotypes
- Haploid, diploid, partial-missing, and multi-allelic calls
- Missing or non-numeric QUAL
- Invalid allele indexes
- `chr` aliases and mitochondrial `M`/`MT`

Hard request errors use:

```json
{
  "error": {
    "code": "STABLE_MACHINE_CODE",
    "message": "Human-readable message",
    "requestId": "..."
  }
}
```

## Frontend behavior

`public/index.html`, `styles.css`, and `app.js` provide:

- Drag/drop and file picker
- Client-side extension and size checks
- Sample discovery for uncompressed VCF
- Optional age and explicit genome-build selection
- Loading and structured error states
- Summary, findings, quality issues, relative scores, warnings
- Safe DOM rendering through `textContent`
- Downloadable JSON report
- Responsive and reduced-motion styles

No third-party browser assets are loaded.

## Security and deployment

- Helmet/CSP and common security headers
- API rate limiting
- Strict multipart file/field/count limits
- Fixed accepted extensions
- Request IDs
- Production errors do not expose stack traces
- Graceful SIGTERM/SIGINT shutdown
- Docker image runs as a non-root user
- Docker health check calls `/api/v1/health`
- Exact dependency versions and lockfile

Deployment still requires TLS, appropriately configured proxy limits/timeouts, log retention controls, monitoring, and an organizational privacy policy.

## Known limitations and recommended next work

Priority 0 — scientific/data governance:

1. Independently review the configured gene–condition scopes.
2. Use ClinVar VCV/RCV XML if per-condition submission dates and complete evidence are required.
3. Left-align variants upstream against the exact GRCh38 reference; backend normalization is minimal representation only.
4. Establish expert review and release approval for database updates.

Priority 1 — scaling:

1. Replace in-memory multipart buffering with a streaming parser if uploads above 25 MB are required.
2. Index only target loci during streaming.
3. Add request cancellation and reverse-proxy timeout tests.

Priority 2 — product:

1. Add a separate upload/inspect step for sample selection in gzip files.
2. Add database filtering and finding search.
3. Generate accessible PDF reports only after report language is professionally reviewed.
4. Add localization if required.

Priority 3 — engineering:

1. Add CI for `npm ci`, `npm run check`, `npm audit`, and Docker build.
2. Add browser end-to-end tests.
3. Add structured metrics without genomic labels or file identifiers.
4. Decide whether to remove legacy binaries/generated artifacts from Git history.

## Commands

```bash
npm ci
npm run dev
npm run check
npm audit
docker build -t genetic-risk-analyzer .
docker run --read-only --tmpfs /tmp -p 3000:3000 genetic-risk-analyzer
```

## Files to read first

1. `README.md`
2. `AI_CONTEXT.md`
3. `server/vcf.js`
4. `server/analyzer.js`
5. `server/app.js`
6. `test/analyzer.test.js`

## Definition of done for future changes

- No genomic contents or filenames logged/persisted.
- Errors remain structured and safe.
- New parsing behavior has edge-case tests.
- No output is described as clinical risk without a validated model.
- `npm run check` passes.
- `npm audit` reports zero known vulnerabilities.
- Docker image builds and health check passes.
- This file is updated with decisions, status, and remaining work.
