.PHONY: all clean 00_scraping_illumina 01_scraping_brainarray_genes 02_selecting_gene_ids

all: out/00_scraping_illumina.dummy out/01_scraping_brainarray_genes.dummy out/02_selecting_gene_ids.dummy
00_scraping_illumina: out/00_scraping_illumina.dummy
01_scraping_brainarray_genes: out/01_scraping_brainarray_genes.dummy
02_selecting_gene_ids: out/02_selecting_gene_ids.dummy

clean:
	rm -rf out


#
# Use dummy files as the output targets to simplify the rules, but still make
# sure that steps are re-run if they fail. If I made one of the TSVs the output
# of each rule, conceivably the R script could fail after writing that file, but
# the pipeline would be complete from the point of view of make so the script
# would not get re-run next time.
#

out/00_scraping_illumina.dummy: 00_scraping_illumina/*
	mkdir -p out/00_scraping_illumina
	cd 00_scraping_illumina && docker build -t illumina-refinery/00_scraping_illumina .
	docker run -t \
		--mount type=bind,source="$(PWD)/out",target=/out \
		illumina-refinery/00_scraping_illumina

	touch out/00_scraping_illumina.dummy


out/01_scraping_brainarray_genes.dummy: 01_scraping_brainarray_genes/*
	mkdir -p out/01_scraping_brainarray_genes
	cd 01_scraping_brainarray_genes && docker build -t illumina-refinery/01_scraping_brainarray_genes .
	docker run -t \
		--mount type=bind,source="$(PWD)/out",target=/out \
		illumina-refinery/01_scraping_brainarray_genes

	touch out/01_scraping_brainarray_genes.dummy


out/02_selecting_gene_ids.dummy: out/00_scraping_illumina.dummy out/01_scraping_brainarray_genes 02_selecting_gene_ids/*
	mkdir -p out/02_selecting_gene_ids
	cd 02_selecting_gene_ids && docker build -t illumina-refinery/02_selecting_gene_ids .
	docker run -t \
		--mount type=bind,source="$(PWD)/out",target=/out \
		illumina-refinery/02_selecting_gene_ids

	touch out/02_selecting_gene_ids.dummy
