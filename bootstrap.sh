#!/usr/bin/env bash
set -euo pipefail

# --- folders ---
mkdir -p .github/workflows \
         config envs scripts/{00_meta,10_clump,20_match,30_scores,40_pathways,50_plink_runs,60_qc_and_viz} \
         workflow data/{raw,external,interim,processed} data/processed/{figures,logs} \
         docs/diagrams

# --- gitignore ---
cat > .gitignore <<'EOF'
.DS_Store
.Rproj.user/
.Rhistory
.RData
.Ruserdata
.vscode/
.idea/
__pycache__/
.env
.env.*
venv/
.venv/
renv/
renv/library/
.Renviron
data/raw/
data/external/*
data/interim/
data/processed/
**/*.log
**/*.sscore
**/*.pgen
**/*.pvar
**/*.psam
**/*.bed
**/*.bim
**/*.fam
EOF

# --- SECURITY.md ---
cat > SECURITY.md <<'EOF'
Do not commit raw or identifiable genetic data. Keep all files in data/* out of version control. 
Scripts operate on local files you place in data/raw and data/external.
EOF

# --- README.md ---
cat > README.md <<'EOF'
# Pathway-specific PRS (APOE-excluded scoring; union meta-GWAS)

Reproducible pipeline to:
- meta-analyze (union) Wightman + Bellenguez sumstats,
- clump to independent SNPs,
- match rsID ↔ chr:pos:ref:alt,
- SNP→gene (SNPXplorer) & GO enrichment (g:Profiler),
- REVIGO similarity + dynamic clustering (k=4),
- PLINK2 scoring (global + clusters), APOE excluded from scoring.

Start: edit config/config.yaml and drop required inputs into data/raw and data/external.  
Run: `make run` (uses Snakemake).
EOF

# --- LICENSE (MIT placeholder) ---
cat > LICENSE <<'EOF'
MIT License
Copyright (c)
Permission is hereby granted, free of charge, to any person obtaining a copy...
EOF

# --- CITATION.cff ---
cat > CITATION.cff <<'EOF'
cff-version: 1.2.0
title: "Pathway-specific PRS pipeline"
message: "If you use this repository, please cite it as below."
authors:
  - family-names: Canelas
    given-names: Carolina Isabel Valentim
date-released: 2025-09-26
license: MIT
EOF

# --- config files ---
cat > config/config.yaml <<'EOF'
paths:
  meta_sumstats_union: data/raw/meta-GWAS-Wightman-Bellenguez-raw.tbl
  meta_sumstats_overlap: data/raw/meta-GWightman-Bellenguez.tbl
  thousand_genomes_eur_plink_prefix: data/raw/1KG_EUR
  dbsnp_b151_vcf_gz: data/raw/common_all_20180423.vcf.gz
  cohort_pfile_prefix: data/raw/COHORT/merged_cohort.hg19
  snpxplorer_annotations: data/external/RESULTS_SnpXplorer/snp_annotation.txt
  gprofiler_intersections: data/external/gProfiler_intersections.csv
  go_terms_bp_all: data/external/go_terms_BP_all.txt
out:
  clump_dir: data/interim/clump/
  matched_dir: data/interim/matched/
  scores_dir: data/interim/scores/
  clusters_dir: data/interim/clusters/
  plink_out_dir: data/processed/plink_scores/
  figures_dir: data/processed/figures/
  logs_dir: data/processed/logs/
params:
  clump_p1: 5e-8
  clump_r2: 0.1
  clump_kb: 250
  maf_min: 0.01
  imputation_r2_min: 0.5
  allele_freq_diff_max: 0.2
  palindromic_maf_max: 0.4
  go_p_adj_threshold_strict: 0.01
  dynamic_tree_min_clusters: 2
  dynamic_tree_max_clusters: 15
  n_clusters_final: 4
tools:
  plink2: plink2
  metal: METAL
  python: python3
  rscript: Rscript
flags:
  exclude_apoe_in_prs: true
  include_apoe_in_mapping: true
EOF

cat > config/env.example <<'EOF'
# Example env; copy to .env and adjust if you really need env vars.
PLINK2=plink2
EOF

# --- envs/conda-env.yaml ---
cat > envs/conda-env.yaml <<'EOF'
name: pathway-prs
channels: [conda-forge, bioconda, defaults]
dependencies:
  - python=3.11
  - snakemake-minimal>=8.0
  - pandas
  - numpy
  - r-base>=4.3
  - r-tidyverse
  - r-data.table
  - r-janitor
  - r-tidytext
  - r-wordcloud2
  - r-dendextend
  - r-dynamicTreeCut
  - r-ggplot2
  - r-cowplot
  - r-glue
  - r-optparse
  - plink2
  - unzip
EOF

# --- workflow/Snakefile ---
cat > workflow/Snakefile <<'EOF'
import yaml
cfg = yaml.safe_load(open("config/config.yaml"))
P, O, T, PM, F = cfg["paths"], cfg["out"], cfg["tools"], cfg["params"], cfg["flags"]

rule all:
    input:
        O["plink_out_dir"] + "global_noAPOE.sscore",
        O["figures_dir"] + "prs_histograms.png"

rule clump:
    input: P["meta_sumstats_union"]
    output: O["clump_dir"] + "clumped.tsv"
    shell:
        """
        {T[plink2]} --bfile {P[thousand_genomes_eur_plink_prefix]} \
          --clump {input} \
          --clump-p1 {PM[clump_p1]} --clump-r2 {PM[clump_r2]} --clump-kb {PM[clump_kb]} \
          --clump-id-field MarkerName --clump-p-field P-value \
          --out {O[clump_dir]}clumped
        mv {O[clump_dir]}clumped.clumped {output}
        """

rule match_to_cohort:
    input:
        clumped=O["clump_dir"] + "clumped.tsv",
        dbsnp=P["dbsnp_b151_vcf_gz"],
        pvar=P["cohort_pfile_prefix"] + ".pvar"
    output: O["matched_dir"] + "matched.tsv"
    shell:
        """
        bash scripts/20_match/extract_snps_from_vcf.sh {input.dbsnp} {input.clumped} {O[matched_dir]}extracted.vcf
        bash scripts/20_match/match_dbsnp_with_cohort.sh {O[matched_dir]}extracted.vcf {input.pvar} {output}
        """

rule build_scores:
    input:
        matched=O["matched_dir"] + "matched.tsv",
        sumstats=P["meta_sumstats_union"]
    output: O["scores_dir"] + "global.score"
    shell:
        """
        {T[rscript]} scripts/30_scores/get_score.R \
          --matched {input.matched} --sumstats {input.sumstats} --out {output} \
          --exclude-apoe {F[exclude_apoe_in_prs]}
        """

rule clusters:
    input:
        snpx=P["snpxplorer_annotations"],
        gprof=P["gprofiler_intersections"]
    output: O["clusters_dir"] + "snp_weighted_mapping_noAPOE_4clust.csv"
    shell:
        """
        {T[python]} scripts/40_pathways/Alternative_REVIGO.py \
          {P[go_terms_bp_all]} {O[clusters_dir]}revigo_lin_distance.txt
        {T[rscript]} scripts/40_pathways/ClusterCreation.R \
          --snpx {input.snpx} --gprof {input.gprof} \
          --revigo {O[clusters_dir]}revigo_lin_distance.txt \
          --p {PM[go_p_adj_threshold_strict]} --k {PM[n_clusters_final]} \
          --out {output}
        """

rule annotate_scores_with_clusters:
    input:
        score=O["scores_dir"] + "global.score",
        weights=O["clusters_dir"] + "snp_weighted_mapping_noAPOE_4clust.csv"
    output: O["scores_dir"] + "weighted_4clust.score"
    shell:
        """
        {T[rscript]} scripts/30_scores/build_plink_scores.R \
          --score {input.score} --weights {input.weights} --out {output}
        """

rule plink_score_global:
    input:
        score=O["scores_dir"] + "weighted_4clust.score",
        pfile_prefix=P["cohort_pfile_prefix"]
    output: O["plink_out_dir"] + "global_noAPOE.sscore"
    shell:
        """
        {T[plink2]} --pfile {input.pfile_prefix} \
          --score {input.score} 1 2 3 cols=+scoresums,+scoreavgs list-variants header \
          --out {O[plink_out_dir]}global_noAPOE
        """

rule qc_and_viz:
    input: O["plink_out_dir"] + "global_noAPOE.sscore"
    output: O["figures_dir"] + "prs_histograms.png"
    shell:
        """
        {T[rscript]} scripts/60_qc_and_viz/ADNI_PLINK_input_output.R \
          --sscore {input} --outdir {O[figures_dir]}
        """
EOF

# --- Makefile ---
cat > workflow/Makefile <<'EOF'
.PHONY: env run clean
env:
	conda env create -f envs/conda-env.yaml || conda env update -f envs/conda-env.yaml
run:
	snakemake -j 8 --use-conda --conda-frontend conda
clean:
	rm -rf data/interim/* data/processed/*
EOF

# --- CI ---
mkdir -p .github/workflows
cat > .github/workflows/ci.yml <<'EOF'
name: CI
on: [push, pull_request]
jobs:
  dryrun:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4
      - uses: conda-incubator/setup-miniconda@v3
        with:
          activate-environment: pathway-prs
          environment-file: envs/conda-env.yaml
          auto-activate-base: false
      - name: Validate Snakefile syntax
        shell: bash -l {0}
        run: snakemake -nq || true
EOF

# --- placeholder scripts so Snakefile doesn't fail immediately ---
for f in scripts/20_match/extract_snps_from_vcf.sh \
         scripts/20_match/match_dbsnp_with_cohort.sh \
         scripts/30_scores/get_score.R \
         scripts/30_scores/build_plink_scores.R \
         scripts/40_pathways/Alternative_REVIGO.py \
         scripts/40_pathways/ClusterCreation.R \
         scripts/60_qc_and_viz/ADNI_PLINK_input_output.R \
         scripts/50_plink_runs/run_plink2_scores.sh \
         scripts/10_clump/run_plink2_clump.sh \
         scripts/00_meta/run_metal_template.txt; do
  if [[ $f == *.sh ]]; then
    cat > "$f" <<'E'
#!/usr/bin/env bash
echo "Placeholder: implement me."
E
    chmod +x "$f"
  elif [[ $f == *.R ]]; then
    cat > "$f" <<'E'
#!/usr/bin/env Rscript
cat("Placeholder: implement me.\n")
E
    chmod +x "$f"
  elif [[ $f == *.py ]]; then
    cat > "$f" <<'E'
#!/usr/bin/env python3
print("Placeholder: implement me.")
E
    chmod +x "$f"
  else
    echo "Template for METAL union meta-analysis (fill with your commands)" > "$f"
  fi
done

# --- docs ---
cat > docs/METHODS.md <<'EOF'
Methods summary lives here. No absolute paths; all parameters sourced from config/config.yaml.
EOF

cat > docs/PIPELINE.md <<'EOF'
Run:
  conda env create -f envs/conda-env.yaml
  conda activate pathway-prs
  make -C workflow run
EOF

cat > docs/PARAMETERS.md <<'EOF'
Key thresholds (clump r2, p<5e-8, GO p-adj<=0.01, dynamicTreeCut 2-15; final k=4). Rationale here.
EOF

cat > docs/diagrams/pipeline.mmd <<'EOF'
flowchart TD
  A[Union meta-GWAS] --> B[PLINK2 clump]
  B --> C[dbSNP match to cohort]
  C --> D[Score construction]
  D --> E[SNP→gene & g:Profiler]
  E --> F[REVIGO Lin similarity]
  F --> G[Dynamic clustering (k=4)]
  G --> H[Weighted PLINK scores]
  H --> I[QC & viz]
EOF

echo "✅ Scaffolding complete."
