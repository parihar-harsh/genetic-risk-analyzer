# Helix Genetic Variant Explorer

A privacy-first educational web application that matches VCF alternate alleles against a reproducibly curated ClinVar panel and reports sourced classifications and genotype quality.

> **Not for clinical use.** Results are not diagnoses, validated polygenic risk scores, or personalized probabilities. Confirm findings through an accredited laboratory and a qualified genetics professional.

## What changed in version 2

The original C++ CLI remains in `genetic_risk_analyzer.cpp` as historical reference. The production application is now a Node.js web service with:

- Responsive, accessible browser interface
- VCF and gzip-compressed VCF uploads
- Exact chromosome, position, REF, and ALT matching
- Required genome-build compatibility with the GRCh38 database
- Minimal representation normalization for small variants
- Multi-sample, phased, haploid, and multi-allelic genotype support
- FORMAT-aware `GT`, `GQ`, and `DP` parsing with INFO/DP fallback
- Structured parsing warnings and quality review flags
- No fake baseline probabilities or unvalidated absolute PRS values
- 11,511 source-matched ClinVar GRCh38 alleles across eight inherited disease groups
- Only pathogenic/likely pathogenic records with at least two-star consensus or expert review
- No penetrance, odds ratio, or PRS values are invented
- In-memory request processing with no genomic-file persistence
- Security headers, request IDs, rate limits, upload limits, and no-store API responses
- Automated unit/integration tests and a non-root Docker image

FASTA and manual DNA entry were intentionally removed. Raw sequence cannot be interpreted as genomic variants without alignment, a reference genome/build, chromosome coordinates, and a proper variant-calling pipeline.

## Local development

Requirements: Node.js 22 or newer.

```bash
npm ci
npm run dev
```

Open `http://localhost:3000`.

Run verification:

```bash
npm run check
npm audit
```

## Production

### Node

```bash
NODE_ENV=production PORT=3000 TRUST_PROXY=1 npm start
```

Use a TLS-terminating reverse proxy or managed container service. Set `TRUST_PROXY` to the exact proxy hop count; do not blindly trust arbitrary forwarded headers.

### Docker

```bash
docker build -t genetic-risk-analyzer .
docker run --read-only --tmpfs /tmp -p 3000:3000 genetic-risk-analyzer
```

Health endpoint:

```text
GET /api/v1/health
```

## Configuration

| Variable | Default | Purpose |
|---|---:|---|
| `NODE_ENV` | `development` | Enables production cache behavior |
| `HOST` | `0.0.0.0` | Listen address |
| `PORT` | `3000` | Listen port |
| `TRUST_PROXY` | `0` | Number of trusted reverse-proxy hops |
| `MAX_UPLOAD_MB` | `25` | Upload limit, constrained to 1–100 MB |
| `RATE_LIMIT_WINDOW_MS` | `900000` | API rate-limit window |
| `RATE_LIMIT_MAX` | `60` | Requests per client per window |

See `.env.example`.

## API

### `POST /api/v1/analyze`

`multipart/form-data` fields:

- `vcf`: required `.vcf` or `.vcf.gz`
- `sampleName`: optional exact sample column name; defaults to the first sample
- `age`: optional integer from 1 through 120
- `genomeBuild`: required when the build cannot be detected; `GRCh37` or `GRCh38`

The response contains input metadata, summary counts, exact findings, ClinVar evidence metadata, and parser warnings. Scores remain disabled because the curated source does not provide a validated patient-level risk model.

## Privacy and operational limits

- Uploaded bytes live only in process memory for the request.
- Filenames and genomic contents are not logged.
- API responses use `Cache-Control: no-store`.
- The default maximum upload is 25 MB. Whole-genome VCFs should be normalized, filtered, and reduced upstream for this small reference database.
- Browser/client and reverse-proxy logs remain the deployer’s responsibility.
- The bundled panel is source-matched to ClinVar GRCh38 release 2026-06-27. ClinVar classifications are submitted assertions, not independent diagnoses.
- `disease_data.metadata.json` records the dataset build, version, provenance state, and clinical-use approval.

## Repository map

```text
public/                  Browser UI
server/app.js            HTTP middleware and routes
server/vcf.js            Defensive VCF parsing
server/analyzer.js       Exact matching and annotations
server/database.js       Database validation and indexing
server/index.js          Process startup and shutdown
test/                    Unit and API integration tests
disease_data.json        Generated ClinVar panel
disease_data.metadata.json  Pinned source/build/release metadata
scripts/curate-clinvar.js   Reproducible curation importer
data/                    Curation and legacy-audit reports
sample_input.vcf         Synthetic passing-call fixture
patient001.vcf           Synthetic mixed-call/filter fixture
genetic_risk_analyzer.*  Legacy C++ implementation
AI_CONTEXT.md            Current state and AI handoff notes
MANUAL_ACTIONS.md         Human scientific and deployment work
```

## License

MIT. See `LICENSE`.
