.PHONY: env run clean
env:
	conda env create -f envs/conda-env.yaml || conda env update -f envs/conda-env.yaml
run:
	snakemake -j 8 --use-conda --conda-frontend conda
clean:
	rm -rf data/interim/* data/processed/*
