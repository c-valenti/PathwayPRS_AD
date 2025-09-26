# Methods for Pathway-Specific Polygenic Risk Score (Pathway-PRS) Analysis

This repository provides **all scripts, workflow logic, and parameter settings** needed to reproduce
the analysis on properly licensed datasets such as ADNI and A4.  
It contains **no real genetic data**.

---

## 1. Genotyping and Quality Control

We constructed pathway-specific polygenic risk scores (PRS) using genome-wide significant
variants reported by **Bellenguez et al. (2022)** and **Wightman et al. (2020)**.

**Data acquisition and preprocessing**

- **Input source**: ADNI genotype data (ADNI1, ADNIGO/ADNI2, ADNI3) obtained under ADNI’s data-use agreement.
- **Imputation**: Minimac4 on the Michigan Imputation Server with the HRC reference panel (v1.1, GRCh37).
- **Array-specific handling**: ADNI genotypes come from three Illumina arrays  
  (Human610-Quad, HumanOmniExpress, HumanOmni2.5-8); imputation was performed separately per array.
- **Pre-imputation QC**: Strand orientation, chromosomal positions, and reference/alternate alleles
  were aligned to HRC.  
  SNPs were removed if any of the following applied:
  - Allele mismatch or unresolved strand flips
  - Allele-frequency difference > 0.2 versus reference
  - Absent from the reference panel
  - Palindromic SNP with minor allele frequency (MAF) > 0.4
- **Post-imputation QC**: Variants with imputation r² < 0.5 or MAF < 0.01 were excluded.

---

## 2. Meta-analysis of GWAS Summary Statistics

To combine Bellenguez and Wightman GWAS results we used **METAL**
([Willer et al., 2010](https://genome.sph.umich.edu/wiki/METAL)).

- We combined the **union** of all SNPs present in either study, not merely the overlap.
- METAL was run with standard inverse-variance weighting to generate
  `meta-GWAS-Wightman-Bellenguez-raw.tbl`, a full meta-analysis table.
  
---
## 3. Clumping of Independent Genome-wide Significant Variants

Independent SNPs were extracted with PLINK 2 clumping using the 1000 Genomes Phase 3 EUR reference.

Exact command line and parameter values (p-value threshold, r², kb window, and field names)
are provided in: [docs/PARAMETERS.md](PARAMETERS.md#3-clumping-of-independent-genome-wide-significant-variants).
  
---

## 4. SNP-to-Cohort Variant Matching

To harmonize GWAS summary statistics with ADNI genotypes we:

- Downloaded dbSNP b151 (GRCh37) common SNP VCF  
  [`common_all_20180423.vcf.gz`](https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/common_all_20180423.vcf.gz).
- Used `extract_snps_from_vcf.sh` to filter to the genome-wide significant variants, removing rare SNPs (MAF < 0.01).
- Used `match_dbsnp_with_adni.sh` to match dbSNP rsIDs to ADNI `chr:pos:ref:alt` variant IDs, automatically resolving strand flips and multi-allelic sites.

Exact filtering parameters and full commands are provided in  
[docs/PARAMETERS.md](PARAMETERS.md#4-snp-to-cohort-variant-matching).

---

## 5. Construction of PLINK Score Files

We built per-subject PRS inputs as follows:

- Merged matched SNPs with meta-analysis summary statistics using `get_score.R`.
- Determined the risk-increasing allele for each variant and converted all effect sizes to positive values.
- Created a PLINK `--score` file with three columns  
  `ADNI_VARIANT_ID   RISK_ALLELE   FINAL_EFFECT`.

These score files form the input to PLINK2 for computing global and cluster-specific PRS.

Exact R logic (allele matching, risk-allele determination, weighted cluster scores)  
and PLINK2 scoring commands are listed in: [docs/PARAMETERS.md](PARAMETERS.md#5-construction-of-plink-score-files).


