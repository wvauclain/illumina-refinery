"""
For each Illumina platform that we support:
  - Read the dataframe of duplicate probes and the different gene IDs
  - Read the brainarray occurrence counts for the current species
  - Group the duplicate probes dataframe by probe ID
  - For each group, pick the gene ID with the highest occurrence count. If there
    are ties, then we want to pick the gene ID that is on the main Ensembl
    assembly. If there are multiple genes on the main assembly, pick the one with
    the lower gene ID to break ties.
"""

import pandas as pd
import requests

ILLUMINA_PLATFORMS = [
    ("illuminaHumanv1", "Homo sapiens"),
    ("illuminaHumanv2", "Homo sapiens"),
    ("illuminaHumanv3", "Homo sapiens"),
    ("illuminaHumanv4", "Homo sapiens"),
    ("illuminaMousev1", "Mus musculus"),
    ("illuminaMousev1p1", "Mus musculus"),
    ("illuminaMousev2", "Mus musculus"),
    ("illuminaRatv1", "Rattus norvegicus"),
]


def pick_gene_id(species: str):
    brainarray_occurrence_counts = pd.read_csv(
        f"/out/01_scraping_brainarray_genes/{species.replace(' ', '_')}.tsv",
        sep="\t",
        index_col="mapped_gene_ids",
    )

    def brainarray_occurrences(gene: str) -> int:
        try:
            return brainarray_occurrence_counts.loc[gene]["n"]
        except KeyError:
            return 0

    def aggregator(args):
        occurrences = [(gene, brainarray_occurrences(gene)) for gene in args.astype(str)]

        nonzero_occurrences = [o for o in occurrences if o[1] > 0]

        if len(nonzero_occurrences) == 0:
            return None
        elif len(nonzero_occurrences) == 1:
            return nonzero_occurrences[0][0]

        # # Now, we know that we have multiple nonzero occurrences, so we need to
        # # check for ties. In this case, I will define a tie as two gene IDs
        # # having at least half the # of occurrences as the max, because there
        # # are some occurrences like [('ENSG00000143248', 29), ('ENSG00000232995', 3)]
        # # which definitely are not ties, but others like
        # # [('ENSG00000230876', 11), ('ENSG00000236854', 7)] which plausibly are
        max_num_occurrences = max(nonzero_occurrences, key=lambda x: x[1])[1]
        # close_occurrences = [o for o in occurrences if 2 * o[1] > max_num_occurrences]

        # if len(close_occurrences) == 1:
        #     # No ties
        #     return close_occurrences[0]

        max_occurrences = [o for o in occurrences if o[1] == max_num_occurrences]
        if len(max_occurrences) == 1:
            return max_occurrences[0][0]

        # Now we have a tie. The first thing to check is that the gene IDs are
        # valid, because there is at least one case of a tie ( ENSG00000157654,
        # ENSG00000243444) where one of the gene IDs is now retired.
        def is_valid(gene):
            response = requests.get(
                f"https://rest.ensembl.org/lookup/id/{gene}?expand=1;content-type=application/json"
            )

            if set(response.json().keys()) == {"error"}:
                print(f"  > Filtered tied gene {gene} because it has an invalid/outdated ID.")
                return False

            return True

        max_occurrences = [o for o in max_occurrences if is_valid(o[0])]
        if len(max_occurrences) == 1:
            return max_occurrences[0][0]

        # Now if we have a tie, we pick the one with the lower gene ID. When I
        # manually inspected all of the ties, none of them were allele pairs and
        # both appeared on the main assembly, so this is an arbitrary but
        # consistent way to pick one gene.
        lowest_id = min(max_occurrences, key=lambda x: x[0])[0]
        print(
            f"  > We have a tie: {', '.join(x[0] for x in max_occurrences)}.\n"
            f"    Picking the one with lowest ID: {lowest_id}."
        )
        return lowest_id

    return aggregator


for platform, species in ILLUMINA_PLATFORMS:
    print(f"Processing platform {platform}...")
    duplicate_probe_ids = pd.read_csv(f"/out/00_scraping_illumina/{platform}.tsv", sep="\t")

    (
        duplicate_probe_ids.groupby("probe_id", as_index=False)
        .aggregate({"ensembl_id": pick_gene_id(species)})
        .to_csv(f"/out/02_selecting_gene_ids/{platform}.tsv", index=False, sep="\t")
    )
