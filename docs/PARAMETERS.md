# Parameters for Pathway-Specific Polygenic Risk Score (Pathway-PRS) Analysis

See the overall narrative in [docs/METHODS.md](METHODS.md#3-clumping-of-independent-genome-wide-significant-variants).

## 3. PLINK Clumping
_See narrative in [Methods §3](METHODS.md#3-clumping-of-independent-genome-wide-significant-variants)._

Independent SNPs were extracted with **PLINK 2** clumping
([documentation](https://zzz.bwh.harvard.edu/plink/clump.shtml)):

```bash
plink2 --bfile 1KG_EUR \
    --clump meta-GWAS-Wightman-Bellenguez-raw.tbl \
    --clump-p1 5e-8 \
    --clump-r2 0.1 \
    --clump-kb 250 \
    --clump-id-field MarkerName \
    --clump-p-field P-value
```
Reference panel: 1000 Genomes Phase 3 EUR
Expected output: ~320 independent genome-wide significant SNPs

## 4. SNP-to-Cohort Variant Matching
_See narrative in [Methods §4]_

Filter dbSNP to genome-wide significant variants and match rsIDs to ADNI chr:pos:ref:alt identifiers:

```bash
bash scripts/20_match/extract_snps_from_vcf.sh \
    data/raw/common_all_20180423.vcf.gz \
    data/raw/meta-GWAS-Wightman-Bellenguez-raw.tbl \
    data/interim/matched/extracted_snps.vcf
```
```bash
bash scripts/20_match/match_dbsnp_with_cohort.sh \
    data/interim/matched/extracted_snps.vcf \
    data/raw/COHORT/merged_cohort.hg19.pvar \
    data/interim/matched/matched.tsv
```
**Outputs**
- `data/interim/matched/extracted_snps.vcf` — dbSNP VCF subset containing only genome-wide significant SNPs (MAF ≥0.01).
- `data/interim/matched/matched.tsv` — TSV mapping each rsID to its ADNI `chr:pos:ref:alt` identifier with REF/ALT alleles and genomic position.


## 5. Construction of PLINK Score Files

_See narrative in [Methods §5](METHODS.md#5-construction-of-plink-score-files)._

Generate the three-column PLINK score file used for PRS calculation:

```bash
Rscript scripts/30_scores/get_score.R \
    --matched data/interim/matched/matched.tsv \
    --sumstats data/raw/meta-GWAS-Wightman-Bellenguez-raw.tbl \
    --out data/interim/scores/global.score \
    --exclude-apoe true
```

## 5.1 PLINK Score Calculation (Global and Cluster-specific)

_See narrative in [Methods §5](METHODS.md#5-construction-of-plink-score-files)._

To compute global and cluster-specific pathway-PRS, run PLINK2 `--score` using the weighted
score files generated in the previous step:

### Global Pathway-PRS
```bash
plink2 \
  --pfile data/raw/COHORT/merged_cohort.hg19 \
  --score data/interim/scores/global_weighted.score 1 2 3 \
  cols=+scoresums,+scoreavgs list-variants header \
  --out data/processed/plink_scores/global_noAPOE
```

Repeat for each cluster (1–4), changing the input and output file names:

### Cluster-specific Pathway-PRS
```bash
plink2 \
  --pfile data/raw/COHORT/merged_cohort.hg19 \
  --score data/interim/scores/weighted_cluster{1..4}.score 1 2 3 \
  cols=+scoresums,+scoreavgs list-variants header \
  --out data/processed/plink_scores/cluster{1..4}_noAPOE
```

---

## 6. Pathway-PRS Cluster Creation
_See narrative in [Methods §6](METHODS.md#6-pathway-prs-cluster-creation)._
The script internally evaluates 2–15 candidate clusters (cutreeDynamic) before fixing the final cluster count to k=4 as reported in Methods.

This stage produces the SNP→cluster weights used to build pathway-specific PRS.  
It consists of: (1) SNP→gene and GO enrichment, (2) GO similarity + clustering, (3) SNP→cluster weighting.

---

### 6.1 Gene-to-Pathway Mapping

**Inputs**
- `data/external/RESULTS_SnpXplorer/snp_annotation.txt` — SNP→gene mapping exported from SNPXplorer (GTEx Whole Blood)  
- `data/external/gProfiler_intersections.csv` — gene–GO intersections from g:Profiler (FDR < 0.01)  
- `data/external/go_terms_BP_all.txt` — list of significant GO BP terms (one per line)

**Command (semantic similarity / REVIGO-style distance)**
```bash
python3 scripts/40_pathways/Alternative_REVIGO.py \
  data/external/go_terms_BP_all.txt \
  data/interim/clusters/revigo_lin_distance.txt
```

**Outputs**
- `data/interim/clusters/revigo_lin_distance.txt` — Lin distance matrix among enriched GO terms


### 6.2 Cluster Creation

**Inputs**
- `data/external/RESULTS_SnpXplorer/snp_annotation.txt`
- `data/external/gProfiler_intersections.csv`
- `data/interim/clusters/revigo_lin_distance.txt`


**Command (hierarchical clustering + dynamic cut):**
```bash
Rscript scripts/40_pathways/ClusterCreation.R \
  --snpx data/external/RESULTS_SnpXplorer/snp_annotation.txt \
  --gprof data/external/gProfiler_intersections.csv \
  --revigo data/interim/clusters/revigo_lin_distance.txt \
  --p 0.01 \
  --k 4 \
  --out data/interim/clusters/snp_weighted_mapping_noAPOE.csv
```
**Outputs**
-`data/interim/clusters/snp_weighted_mapping_noAPOE.csv` — per-SNP weights across the 4 clusters 

**Notes**
--p 0.01 enforces the g:Profiler significance threshold (FDR < 0.01).
--k 4 fixes the final number of clusters to 4 after dynamic search (2–15) as in Methods.


### 6.3 Wiring SNP–Cluster Weights into PRS
These weights feed into the score-building step to produce weighted effect sizes per SNP.

**Inputs**
- `data/interim/scores/global.score` — base 3-column score file (SNP, A1, BETA) from §5
- `data/interim/clusters/snp_weighted_mapping_noAPOE.csv` — SNP→cluster weights (this section)

**Command (annotate scores with cluster weights):**
```bash
Rscript scripts/30_scores/build_plink_scores.R \
  --score data/interim/scores/global.score \
  --weights data/interim/clusters/snp_weighted_mapping_noAPOE.csv \
  --out-prefix data/interim/scores/
```

**Outputs**
- `data/interim/scores/global_weighted.score` — global PRS weighted across clusters
- `data/interim/scores/weighted_cluster1.score`
- `data/interim/scores/weighted_cluster2.score`
- `data/interim/scores/weighted_cluster3.score`
- `data/interim/scores/weighted_cluster4.score`

These weighted score files are the direct inputs for PLINK2 --score in §5.1.

---

## 7. Pathway-PRS Calculation with PLINK
_See narrative in [Methods §7](METHODS.md#7-pathway-prs-calculation-with-plink)._

This stage turns weighted SNP–cluster mappings into per-individual PRS values.

### 7.1 Weighted Score File Creation
Combine the three-column score file (`global.score`) with SNP–cluster weights:

```bash
Rscript scripts/30_scores/build_plink_scores.R \
  --score   data/interim/scores/global.score \
  --weights data/interim/clusters/snp_weighted_mapping_noAPOE.csv \
  --out-prefix data/interim/scores/
  ```

**Outputs:**
- `data/interim/scores/global_weighted.score`
- `data/interim/scores/weighted_cluster1.score`
- `data/interim/scores/weighted_cluster2.score`
- `data/interim/scores/weighted_cluster3.score`
- `data/interim/scores/weighted_cluster4.score`


***FinalEffectWeighted*** is computed as:
- **FinalEffectWeighted** = FinalEffect * WeightedFactor

### 7.2 PLINK2 Scoring
Compute PRS for global and cluster-specific scores.

1. **Global PRS:**
```bash
plink2 \
  --pfile data/raw/COHORT/merged_cohort.hg19 \
  --score data/interim/scores/global_weighted.score 1 2 3 \
  cols=+scoresums,+scoreavgs list-variants header \
  --out data/processed/plink_scores/global_noAPOE
```

2. **Cluster-specific PRS** (repeat for each cluster):
```bash
plink2 \
  --pfile data/raw/COHORT/merged_cohort.hg19 \
  --score data/interim/scores/weighted_cluster{1..4}.score 1 2 3 \
  cols=+scoresums,+scoreavgs list-variants header \
  --out data/processed/plink_scores/cluster{1..4}_noAPOE
```

**Outputs:**
- `data/processed/plink_scores/global_noAPOE.sscore`
- `data/processed/plink_scores/cluster1_noAPOE.sscore … cluster4_noAPOE.sscore`

### 7.3 QC and Visualisation
Process `.sscore` outputs in R for QC and integration:

```bash
Rscript scripts/60_qc_and_viz/ADNI_PLINK_input_output.R \
  --sscore data/processed/plink_scores/global_noAPOE.sscore \
  --outdir data/processed/figures/
  ```

- Generates:
1. Scaled histograms per cluster (e.g. hist_scaled_cluster1.png)
2. Boxplots of SCORE1_AVG distributions
3. Combined long-format PRS dataset for downstream analysis