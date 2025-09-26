# Methods for Pathway-Specific Polygenic Risk Score (Pathway-PRS) Analysis

This repository provides **all scripts, workflow logic, and parameter settings** needed to reproduce
the analysis on properly licensed datasets such as ADNI and A4.  
It contains **no real genetic data**.

See exact command and thresholds in [Parameters §3 PLINK clumping](PARAMETERS.md#3-plink-clumping).

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

*** Scripts: ***
[extract_snps_from_vcf.sh](../scripts/20_match/extract_snps_from_vcf.sh),
[match_dbsnp_with_cohort.sh](../scripts/20_match/match_dbsnp_with_cohort.sh).

The first script produces a compact VCF (`extracted_snps.vcf`) containing only meta-GWAS significant SNPs present in dbSNP.  
The second script produces a tab-separated table (`matched.tsv`) mapping each GWAS rsID to the corresponding ADNI `chr:pos:ref:alt` identifier, including reference and alternate alleles and genomic position.

---

## 5. Construction of PLINK Score Files

We built per-subject PRS inputs as follows:

- Merged matched SNPs with meta-analysis summary statistics using `get_score.R`.
- Determined the risk-increasing allele for each variant and converted all effect sizes to positive values.
- Created a PLINK `--score` file with three columns  
  `ADNI_VARIANT_ID   RISK_ALLELE   FINAL_EFFECT`.

***R script:*** [get_score.R](../scripts/30_scores/get_score.R).

These score files form the input to PLINK2 for computing global and cluster-specific PRS. Exact R logic (allele matching, risk-allele determination, weighted cluster scores) and PLINK2 scoring commands are listed in: [docs/PARAMETERS.md](PARAMETERS.md#5-construction-of-plink-score-files).

### 5.1 PRS Computation (Global and Cluster-specific)

Using the cleaned and weighted score files, we computed polygenic risk scores (PRS)
for each ADNI participant with PLINK2.  
Both a **global** PRS and **four pathway-specific cluster PRS** were calculated by
running `plink2 --score` once per score file.

All commands and parameter settings (e.g., `cols=+scoresums,+scoreavgs`)
are detailed in [docs/PARAMETERS.md](PARAMETERS.md#51-plink-score-calculation-global-and-cluster-specific).

---

## 6. Pathway-PRS Cluster Creation
This module identifies biologically meaningful pathway clusters from genome-wide significant SNPs and computes weighted SNP-to-cluster assignments for pathway-specific PRS.

It consists of three sub-steps: (i) SNP->gene mapping and GO enrichment, (ii) clustering of enriched pathways, and (iii) weighted mapping of SNPs to pathway clusters.

See exact scripts and command-line parameters in
[docs/PARAMETERS.md](PARAMETERS.md#6-pathway-prs-cluster-creation).

### 6.1 Gene-to-Pathway Mapping

1. SNP -> gene annotation
We used the Gene-set Enrichment Analysis function in SNPXplorer with tissue context GTEx – Whole Blood.
- Input: ≈291 independent significant SNPs (after PLINK clumping).
- Output: a SNP–gene annotation table (data/external/RESULTS_SnpXplorer/snp_annotation.txt).

2. GO functional enrichment
Genes mapped by SNPXplorer were tested for Gene Ontology (GO) Biological Process enrichment using g:Profiler (Raudvere et al., 2019).
- Only pathways with FDR-adjusted p < 0.01 were retained.
- Output: gene–GO intersections (data/external/gProfiler_intersections.csv).

3. Preparation for similarity analysis
Significant GO terms were exported as go_terms_BP_all.txt for downstream semantic similarity calculations.

### 6.2 Similarity Matrix and Cluster Creation
To identify biologically coherent pathway groups:

1. REVIGO-style similarity matrix
We used a customised version of the approach in Tesi et al. (2021) implemented in scripts/40_pathways/Alternative_REVIGO.py.
- This script computes Lin semantic similarity (Lin et al., 1998) among enriched GO terms and outputs a distance matrix (revigo_lin_distance.txt).

2. **Hierarchical clustering**
The distance matrix was clustered in R (scripts/40_pathways/ClusterCreation.R) using Ward.D2 linkage and the cutreeDynamic algorithm (2–15 candidate clusters).
- This yielded four stable clusters representing major biological themes: **1)** Immune activation, **2)**  Leukocyte regulation, **3)** Amyloid processing and **4)** Stimulus transduction

3. **Visualization**
Cluster dendrograms and cluster-specific word clouds were generated to highlight key GO terms and their relationships.
- Background terms were filtered (≥5 % frequency) and tokenised with tidytext to enhance interpretability.

### 6.3 Weighted Mapping of SNPs to GO Clusters
To link SNPs quantitatively to pathway clusters:

1. **Filtering and QC** 
- Removed APOE locus variants (rs429358, rs7412) to prevent single-locus dominance.
- Removed SNPs mapping exclusively to the IGH locus to avoid spurious immune enrichment.
- Retained 254 SNPs (212 genes) after filtering.

2. **Merging SNP–gene–pathway relationships**
- The cleaned SNP–gene table was intersected with significant gene–GO associations from Gprofiler.
- Two aggregated views were generated: **1)** ***snp2gene2pathway_01.csv*** (SNP->gene->cluster) and **2)** ***pathway2gene2snp_01_4clust.csv*** (cluster->gene->SNP)

3. **Computing SNP–cluster weights**
For each SNP, a weighted mapping factor was calculated as: (#links to a given cluster) ÷ (total #links across all clusters).

- The resulting file ***snp_weighted_mapping_noAPOE_.csv*** lists all SNPs with weights across the four clusters and serves as input for pathway-specific PRS computation.

---
All scripts (Alternative_REVIGO.py, ClusterCreation.R) and parameter settings for these steps are listed with runnable examples in
docs/PARAMETERS.md