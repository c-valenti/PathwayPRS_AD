# Data directory

This repository ships with **example placeholder files** under `data/derived/` to illustrate required formats for running the pipeline.  
All rows are synthetic or reduced to a single example SNP and **do not contain real data**.

| Path                                                  | Purpose                       | Format                                        |
|-------------------------------------------------------|-------------------------------|-----------------------------------------------|
| data/derived/clump/clumped.tsv                        | Example PLINK clump output    | Columns: SNP, P |
| data/derived/matched/matched.tsv                      | Example dbSNP→cohort mapping  | RSID, ADNI_VARIANT_ID, REF, ALT, CHR, POS |
| data/derived/scores/global.score                      | Example PLINK score template  | ADNI_ID, A1, FINAL_EFFECT |
| data/derived/clusters/snp_weighted_mapping_noAPOE.csv | Example cluster weights       | SNP, cluster, weighted_factor |
