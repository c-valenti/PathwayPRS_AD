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
