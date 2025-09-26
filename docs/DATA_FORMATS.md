# Data Formats (Schema Only)

This repo ships **no data**. Provide your own files matching the schemas below.
Paths correspond to `config/config.yaml`.

---

## Inputs you must provide

### Union meta-GWAS (Bellenguez + Wightman)
**Path:** `data/raw/meta-GWAS-Wightman-Bellenguez-raw.tbl`  
**Format:** TSV (tab-delimited), header required  
**Minimum columns (case-sensitive):**
- `MarkerName` (rsID)
- `Allele1` (effect allele in source GWAS)
- `Allele2` (other allele)
- `Effect` (β or log(OR))
- `StdErr`
- `P-value`

**Tiny example**
MarkerName Allele1 Allele2 Effect StdErr P-value
rs429358 C T 0.35 0.05 1e-50


### dbSNP b151 (GRCh37/hg19)
**Path:** `data/raw/common_all_20180423.vcf.gz`  
**Format:** VCF 4.1/4.2 (bgzipped + tabix-indexed recommended)

### Cohort genotypes (PLINK2 pfile set, hg19)
**Prefix:** `data/raw/COHORT/merged_cohort.hg19`  
**Files:** `.pgen`, `.pvar`, `.psam`  

**`.pvar` header example**
##fileformat=VCFv4.2
#CHROM POS ID REF ALT QUAL FILTER INFO
19 45411941 19:45411941:T:C T C . PASS .

---

## Intermediate files produced by the workflow

- `data/interim/clump/clumped.tsv`  
  Columns (example): `SNP`, `P`
- `data/interim/matched/matched.tsv`  
  Columns: `RSID`, `ADNI_VARIANT_ID`, `REF`, `ALT`, `CHR`, `POS`
- `data/interim/scores/global.score` (PLINK2 `--score` 3-col file)  
  Columns: `SNP(ADNI_ID)`, `A1`, `BETA`
- `data/interim/clusters/snp_weighted_mapping_noAPOE_4clust.csv`  
  Columns: `SNP`, `cluster`, `weighted_factor`

---

## Outputs (when running full pipeline)

- `data/processed/plink_scores/global_noAPOE.sscore`  
  Standard PLINK2 `.sscore` with `SCORE1_SUM`, `SCORE1_AVG`, etc.
- `data/processed/figures/prs_histograms.png`  
  QC visualization.

For full commands and thresholds see [`docs/PARAMETERS.md`](PARAMETERS.md).
For the narrative steps see [`docs/METHODS.md`](METHODS.md).

