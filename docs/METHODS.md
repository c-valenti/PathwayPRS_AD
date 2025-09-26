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

This module identifies biologically meaningful pathway clusters from genome-wide significant SNPs and computes weighted SNP-to-cluster assignments for pathway-specific PRS. It consists of three sub-steps: (i) SNP→gene mapping and GO enrichment, (ii) clustering of enriched pathways, and (iii) weighted mapping of SNPs to pathway clusters.

See exact scripts and command-line parameters in
[docs/PARAMETERS.md](PARAMETERS.md#6-pathway-prs-cluster-creation).


### 6.1 Gene-to-Pathway Mapping

1. **SNP→gene annotation**  
   We used the *Gene-set Enrichment Analysis* function in **SNPXplorer** (GTEx: **Whole Blood**).  
   - **Input:** ~291 independent significant SNPs (post-PLINK clumping).  
   - **Output:** SNP–gene annotation table: `data/external/RESULTS_SnpXplorer/snp_annotation.txt`.

2. **GO functional enrichment**  
   Genes mapped by SNPXplorer were tested for GO Biological Process enrichment using **g:Profiler** (Raudvere et al., 2019).  
   - **Threshold:** FDR-adjusted *p* < 0.01 (retained terms only).  
   - **Output:** gene–GO intersections: `data/external/gProfiler_intersections.csv`.

3. **Preparation for similarity analysis**  
   Significant GO BP terms were exported for semantic similarity computation as:  
   `data/external/go_terms_BP_all.txt`.


### 6.2 Similarity Matrix and Cluster Creation

1. **REVIGO-style similarity matrix**  
   We implemented a customised version of the Tesi et al. (2021) approach in  
   `scripts/40_pathways/Alternative_REVIGO.py`, which computes **Lin** semantic similarity (Lin, 1998) among enriched GO terms and writes a distance matrix:  
   `data/interim/clusters/revigo_lin_distance.txt`.

2. **Hierarchical clustering**  
   The distance matrix was clustered in R (`scripts/40_pathways/ClusterCreation.R`) using **Ward.D2** linkage and **cutreeDynamic** (searching 2–15 clusters).  
   - **Result:** four stable clusters representing major biological themes:  
     1) Immune activation, 2) Leukocyte regulation, 3) Amyloid processing, 4) Stimulus transduction.

3. **Visualisation**  
   Cluster dendrograms and cluster-specific word clouds were generated to highlight key GO terms and their relationships. Background tokens (≥5% frequency) were filtered; tokenisation used `tidytext`.


### 6.3 Weighted Mapping of SNPs to GO Clusters

1. **Filtering and QC**  
   - Removed **APOE** locus variants (rs429358, rs7412) to prevent single-locus dominance.  
   - Removed SNPs mapping exclusively to the **IGH** locus to avoid spurious immune enrichment.  
   - **Retained:** 254 SNPs (212 genes).

2. **Merging SNP–gene–pathway relationships**  
   The cleaned SNP–gene table was intersected with significant gene–GO associations (g:Profiler *q* < 0.01). Two aggregated views were produced:  
   - `data/interim/clusters/snp2gene2pathway_01.csv` (SNP → gene → cluster)  
   - `data/interim/clusters/pathway2gene2snp_01.csv` (cluster → gene → SNP)

3. **Computing SNP–cluster weights**  
   For each SNP, the **weighted mapping factor** to a cluster was defined as:  
   *(# links to that cluster) ÷ (total # links across all clusters)*.  
   - **Output:** `data/interim/clusters/snp_weighted_mapping_noAPOE.csv`, providing per-SNP weights across all four clusters (used downstream for pathway-specific PRS).

---

## 7. Pathway-PRS Calculation with PLINK
_This section follows the cluster creation described in [Methods §6](METHODS.md#6-pathway-prs-cluster-creation)._

This stage converts the weighted SNP–cluster mappings from §6 into global and cluster-specific polygenic risk scores (PRS) for each ADNI participant.
See exact scripts and command-line parameters in [docs/PARAMETERS.md](PARAMETERS.md#7-pathway-prs-calculation-with-plink).

### 7.1 Rationale and APOE Handling
To capture polygenic effects beyond the strong APOE locus, we:
- **Included** APOE SNPs in gene/pathway mapping and cluster construction (§6),
  ensuring biologically accurate clusters.
- **Excluded** APOE SNPs from PRS scoring, using the filtered
  `snp_weighted_mapping_noAPOE.csv` to prevent single-locus dominance.

### 7.2 Weighted PRS Score Files
We merged the three-column PLINK score file from §5 (`global.score`, containing ADNI SNP ID, effect allele, and Final_Effect) with the cluster weights from §6 (`snp_weighted_mapping_noAPOE.csv`).

- For each SNP and each cluster:
\[\text{FinalEffectWeighted} = \text{FinalEffect} \times \text{WeightedFactor}\]

**Outputs:**
- `global_weighted.score` – global PRS (all SNPs weighted by their cluster factors)
- `weighted_cluster1.score` ... `weighted_cluster4.score` – one file per cluster

### 7.3 PLINK2 Scoring
We ran PLINK2 `--score` to compute per-individual PRS values.
For each cluster (1–4) the same command was repeated, substituting thecorresponding cluster-specific weighted_clusterX.score and output prefix.

**Example (global PRS):**
```bash
plink2 \
  --pfile data/raw/COHORT/merged_cohort.hg19 \
  --score data/interim/scores/global_weighted.score 1 2 3 \
  cols=+scoresums,+scoreavgs list-variants header \
  --out data/processed/plink_scores/global_noAPOE
  ```
  
**PLINK outputs** .sscore files containing:
- SCORE1_SUM – sum of weighted allele dosages per individual 
- SCORE1_AVG – mean weighted score per SNP


### 7.4 Quality Control and Integration
Using `scripts/60_qc_and_viz/ADNI_PLINK_input_output.R` we:
- Removed duplicate individuals and scaled SCORE1_AVG to a 0–1 range
- Combined global and cluster-specific PRS into one tidy dataset
- Generated histograms and boxplots of PRS distributions for QC
- Reported summary statistics (mean/SD per cluster, number of valid scores)

Final cluster SNP counts (after APOE removal):
- Cluster 1: 88 SNPs
- Cluster 2: 55 SNPs
- Cluster 3: 105 SNPs 
- Cluster 4: 157 SNPs