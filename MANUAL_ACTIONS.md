# Manual Actions Required for Accuracy

The backend and bundled panel now have reproducible source provenance, but source matching is not independent clinical validation. Complete the following work before describing this project as clinically or scientifically validated.

## 1. Establish the intended use

Choose and document one:

- Educational variant explorer
- Research-only annotation tool under an approved protocol
- Clinical decision-support software

The current system is approved only for the first option. Clinical use requires qualified medical, regulatory, privacy, quality-system, and security review. Do not change `clinicalUseApproved` to `true` yourself unless the applicable approval process is complete and documented.

## 2. Dataset curation completed

The original 94 project-authored records were audited and replaced. Only two matched the ClinVar GRCh38 allele coordinates, and both contradicted the claimed gene or condition. See `data/legacy-dataset-audit.json`.

The generated database now contains 11,511 records from ClinVar GRCh38 release 2026-06-27. The importer requires:

1. A sequence REF/ALT allele on GRCh38.
2. Pathogenic, likely pathogenic, or combined pathogenic/likely pathogenic classification.
3. At least two-star multiple-submitter consensus, expert-panel review, or practice guideline review.
4. A matching ClinVar condition term.
5. An established gene within the configured eight disease scopes.
6. No conflicting aggregate classification.

Every entry records its ClinVar Variation ID, Allele ID, SCV accessions, review status, review stars, release date, HGVS expression, condition identifiers, source URL, and matched gene.

Reproduce a release with:

```bash
curl -o /tmp/clinvar-GRCh38.vcf.gz \
  https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar.vcf.gz
node scripts/curate-clinvar.js \
  --input /tmp/clinvar-GRCh38.vcf.gz \
  --expected-md5 27166aac6887cb7ac10829de69587993 \
  --trusted-panel \
  --apply
```

The remaining manual scientific work is independent review of the selected scope. The ClinVar VCF is summary-level and does not provide per-condition evaluation dates or complete submission evidence. Use the VCV/RCV XML releases if that level of audit is required.

Run:

```bash
npm run audit:data
npm run audit:data:strict
```

The strict command now passes for source provenance. It does not mean the application is approved for clinical use.

## 3. Do not validate effect estimates casually

Set `effectEstimateValidated` to `true` only when all of these are documented:

- Exact effect allele and other allele
- Genome build and variant normalization
- Published source and model version
- Target population/ancestry
- Covariates and baseline model
- Units and scale of the weight
- Independent validation cohort
- Intended disease endpoint

An odds ratio from one publication is not a PRS. Do not combine odds ratios from unrelated studies into an absolute probability.

## 4. Prepare input VCFs correctly

Before upload, you need to:

1. Confirm the sample and genome build.
2. Normalize and split multiallelic records against the exact GRCh38 reference.
3. Check REF alleles against that reference.
4. Apply the originating caller’s validated FILTER/GQ/DP criteria.
5. Confirm sample identity and contamination checks.
6. Upload only the intended sample.

Typical upstream commands, after installing `bcftools` and obtaining the matching reference:

```bash
bcftools norm -f GRCh38.fa -m -any input.vcf.gz -Oz -o normalized.vcf.gz
bcftools index -t normalized.vcf.gz
bcftools view -f PASS,. -s SAMPLE_ID normalized.vcf.gz -Oz -o sample.pass.vcf.gz
```

Do not run these commands with an arbitrary reference FASTA. The contig names and assembly must match the VCF.

## 5. Add clinical interpretation rules only with experts

The backend deliberately does not infer that:

- One pathogenic allele causes a recessive disease
- A VUS is disease-causing
- A benign/common risk allele is diagnostic
- Lack of a database match means a negative result
- An age penetrance table is a personal probability

Adding inheritance, phase, sex-chromosome, mosaicism, CNV, repeat-expansion, mitochondrial heteroplasmy, or phenotype logic requires domain-specific validated rules and expert sign-off.

## 6. Operational work before public deployment

- Build and run the Docker image in CI.
- Put the service behind TLS.
- Configure exact proxy trust, body limits, and timeouts.
- Add authentication if reports are retained or shared.
- Verify that proxy, CDN, analytics, and application logs contain no genomic data or filenames.
- Define retention, deletion, incident-response, and consent policies.
- Complete threat modeling and dependency/container scanning.
- Add monitoring that uses only non-genomic metrics.
- Obtain legal/privacy review for the jurisdictions where it will run.

## 7. Release checklist

Before each dataset or application release:

```bash
npm ci
npm run check
npm audit
npm run audit:data:strict
docker build -t genetic-risk-analyzer .
```

Then perform a documented review with known positive, negative, multiallelic, filtered, build-mismatch, and no-call reference samples.

## Authoritative references

- [VCF specification](https://samtools.github.io/hts-specs/VCFv4.5.pdf)
- [ClinVar downloads and release scope](https://www.ncbi.nlm.nih.gov/clinvar/docs/downloads/)
- [ClinVar review status](https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/)
- [Using ClinVar data programmatically](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/)
