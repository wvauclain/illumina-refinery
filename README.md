# Illumina Refinery

A workflow for selecting which gene ID to map a given probe ID to in the case
where a given probe maps to multiple genes.

## Methods

To determine which gene to map a given probe to, we select the genes in order of priority as follows:

- Pick the gene with the most appearances in Brainarray packages for Affymetrix platforms on the same species as the input Illumina platform
- If none of the associated gene IDs appear in any Brainarray platform, an NA is emitted
- If two or more of the associated gene IDs appear an equal number of times in
  Brainarray platforms, the gene ID with the lower number is selected to break the tie.

The workflow is split up into the following steps:

### 00 Scraping Illumina

For each Illumina platform we are interested in, we do the following:

- Load relevant Bioconductor package
- Load the probe-to-gene mappings into a data frame
- Filter the data frame so that the remaining columns are all probe IDs which appear more than once
- Output the data frame as a TSV file for each platform

## 01 Scraping Brainarray

This section of the workflow counts the number of times each gene ID appears among all the brainarray packages for a given species, and outputs this count as a TSV file for the three species with Illumina platforms we are interested in.

### 02 Selecting Gene IDs

This step reads in the outputs from the previous steps and applies the prioritization from above to determine how to map each probe ID. The output is given as a TSV file for each platform so that we can use the results on refine.bio.
