# vcf2phylip audit roadmap

This document tracks the first `plant-genome-audit` branch for a fork of `vcf2phylip`.

The goal is not to replace the original tool. The goal is to extend it with audit-oriented features for complex, non-model, repeat-rich, and plant genomes where a SNP matrix can be affected by missing data, ambiguous mapping, paralogy, repeats, and reference bias.

## Scope of the first branch

Branch name:

```text
plant-genome-audit
```

Initial philosophy:

```text
Before making a tree, prove the SNP matrix is not cursed.
```

## Completed starter changes

### 1. Add audit roadmap

This file documents the rationale, implemented features, and future development plan.

### 2. Add small synthetic VCF test cases

Synthetic test files are included under:

```text
tests/data/
```

These test files cover:

- biallelic SNPs
- multiallelic SNPs
- non-PASS FILTER values
- missing genotypes
- BED exclusion
- distance-based thinning
- multi-digit allele IDs such as `10/10`

### 3. Fix genotype allele ID parsing

The original logic removed separators from genotype strings before parsing alleles. That works for simple allele IDs such as `0/1`, but allele IDs greater than 9 are unsafe.

Example problem:

```text
10/10
```

Unsafe character-based parsing sees this as:

```text
1, 0, 1, 0
```

Safe allele-token parsing sees this as:

```text
10, 10
```

The audit branch now parses genotype fields using `/` and `|` as separators while preserving full allele tokens.

### 4. Add reproducible IUPAC resolution

New option:

```bash
--seed 42
```

This only affects runs using:

```bash
--resolve-IUPAC
```

The random allele resolver is now reproducible when a seed is supplied.

### 5. Add basic audit summary report

New option:

```bash
--audit-summary
```

This writes:

```text
<prefix>.audit_summary.json
<prefix>.audit_summary.tsv
```

Tracked counts include records processed, accepted records, and skipped records by reason.

### 6. Add PASS-only filtering

New option:

```bash
--filter-pass-only
```

This keeps only records with:

```text
FILTER=PASS
FILTER=.
```

### 7. Add biallelic-only filtering

New option:

```bash
--biallelic-only
```

This keeps only records with exactly one ALT allele after lightweight `<NON_REF>` normalization.

### 8. Add repeat/blacklist BED exclusion

New option:

```bash
--exclude-bed repeats_or_blacklist.bed
```

BED files are interpreted as standard 0-based, half-open intervals. VCF `POS` is interpreted as 1-based.

This is intended for excluding repeats, low-complexity regions, problematic regions, or user-defined blacklists.

### 9. Add distance-based SNP thinning

New option:

```bash
--thin-distance 1000
```

This keeps retained SNPs at least the specified number of base pairs apart on the same contig.

This is useful when clustered SNPs may represent local mapping artifacts, repeats, or overly linked regions.

## Recommended first test command

```bash
python vcf2phylip.py \
  -i tests/data/synthetic_basic.vcf \
  --output-folder test_output \
  --output-prefix synthetic_basic \
  --filter-pass-only \
  --biallelic-only \
  --exclude-bed tests/data/exclude_regions.bed \
  --thin-distance 25 \
  --audit-summary \
  --write-used-sites \
  -f
```

Expected behavior:

- low-quality FILTER record is skipped
- multiallelic record is skipped by `--biallelic-only`
- BED-overlapping record is skipped
- one close SNP is skipped by `--thin-distance`
- accepted SNPs are written to PHYLIP and FASTA
- audit reports are written as JSON and TSV

## Future features for MANDRAGORA integration

Potential future commands:

```bash
mandragora vcf-audit
mandragora vcf2phylo
mandragora snp-filter-report
```

Future audit layers:

- per-sample missingness report
- per-sample heterozygosity report
- per-locus heterozygosity/paralog warning
- depth-aware filtering using FORMAT/DP
- genotype-quality filtering using FORMAT/GQ
- allele-balance filtering for diploid genotypes
- contig-level SNP density report
- linked-SNP pruning options
- invariant-site/ascertainment-bias warning for IQ-TREE/RAxML-style analysis
- optional HTML report
- optional integration with RepeatMasker/BED outputs

## Recommended Git commit sequence

```text
docs: add audit roadmap

test: add synthetic VCF and BED fixtures

fix: parse genotype allele IDs safely

feat: add seeded IUPAC resolution

feat: add audit summary reports

feat: add PASS-only filtering

feat: add biallelic-only filtering

feat: add BED exclusion for blacklist regions

feat: add distance-based SNP thinning
```

## Development principle

Keep the first fork respectful and small. The original script remains the foundation. This branch adds safeguards for cases where variant-derived phylogenetic matrices require stronger auditability before biological interpretation.
