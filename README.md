### Pathway-specific PRS (APOE-excluded scoring; union meta-GWAS)

Reproducible pipeline to:
- meta-analyze (union) Wightman + Bellenguez sumstats,
- clump to independent SNPs,
- match rsID ↔ chr:pos:ref:alt,
- SNP→gene (SNPXplorer) & GO enrichment (g:Profiler),
- REVIGO similarity + dynamic clustering (k=4),
- PLINK2 scoring (global + clusters), APOE excluded from scoring.

Start: edit config/config.yaml and drop required inputs into data/raw and data/external.  
Run: `make run` (uses Snakemake).

# Quick links
- [Methods (narrative)](docs/METHODS.md)
- [Parameters & commands](docs/PARAMETERS.md)
- [Data formats (schemas only)](docs/DATA_FORMATS.md)


