# Parameters for Pathway-Specific Polygenic Risk Score (Pathway-PRS) Analysis

See the overall narrative in [docs/METHODS.md](METHODS.md#3-clumping-of-independent-genome-wide-significant-variants).

## 3. PLINK Clumping

See narrative in [Methods §3](METHODS.md#3-clumping-of-independent-genome-wide-significant-variants)._

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

Filter dbSNP to genome-wide significant variants:

```bash
bash scripts/20_match/extract_snps_from_vcf.sh \
    common_all_20180423.vcf.gz \
    meta-GWAS-Wightman-Bellenguez-raw.tbl \
    extracted_snps.vcf
```

## 5. Construction of PLINK Score Files

See narrative in [Methods §5](METHODS.md#5-construction-of-plink-score-files).

Generate the three-column PLINK score file used for PRS calculation:

```bash
Rscript scripts/30_scores/get_score.R \
    --matched data/interim/matched/matched.tsv \
    --sumstats data/raw/meta-GWAS-Wightman-Bellenguez-raw.tbl \
    --out data/interim/scores/global.score \
    --exclude-apoe true
