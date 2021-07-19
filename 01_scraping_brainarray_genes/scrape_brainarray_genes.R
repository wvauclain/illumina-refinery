suppressPackageStartupMessages(library(xml2))
suppressPackageStartupMessages(library(lazyeval))
suppressPackageStartupMessages(library(dplyr))

ensg_url <- "http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/22.0.0/ensg.asp"
html_content <- read_html(ensg_url)

data_table <- xml_find_all(html_content, ".//table")[[2]]
data_rows <- xml_children(data_table)[3:xml_length(data_table)]

species <- c()
chips <- c()
packages <- c()

# Get the chip and species out of every row
invisible(lapply(data_rows, function(row) {
    # Column #2: "Species"
    species_name <- xml_text(xml_child(row, 2))

    # Column #3: "Chip"
    chip_name <- xml_text(xml_child(row, 3))
    # Convert chip_name to lower case, then normalize it by removing
    # "-", "_" and " " characters.
    chip_name <- gsub("[-_ ]", "", tolower(chip_name))

    # Column #16: "R Source Packages (C and P)".
    # We need the URL of "P", which is the second URL.
    source_pkgs <- xml_child(row, 16)
    probe_pkg_url <- xml_attr(xml_child(source_pkgs, 2), "href")
    pkg_filename <- tail(strsplit(probe_pkg_url, "/")[[1]], n=1)
    pkg_name <- strsplit(pkg_filename, '_')[[1]][1]

    if (nzchar(species_name) && nzchar(chip_name) && nzchar(pkg_name)) {
        species <<- c(species, species_name)
        chips <<- c(chips, chip_name)
        packages <<- c(packages, pkg_name)
    }
}))

species_and_chips <- data.frame(species, chip=chips, package=packages)

for (species in c("Homo_sapiens", "Mus_musculus", "Rattus_norvegicus")) {
    write(paste0("Scraping common gene ids for ", species, "..."), file=stdout())

    mapped_gene_ids = c()

    invisible(apply(species_and_chips[species_and_chips$species == species,], 1, function(row) {
        species <- row[1]
        chip <- row[2]
        package <- row[3]

        write(paste0("  > Processing chip ", chip, "..."), file=stdout())

        suppressPackageStartupMessages(library(package, character.only=TRUE))

        brainarray_df <- as.data.frame(get(package), stringsAsFactors = FALSE)
        brainarray_df <- brainarray_df %>%
            dplyr::mutate(Probe.Set.Name = sub("_at", "", Probe.Set.Name)) %>%
            dplyr::filter(startsWith(Probe.Set.Name, "ENSG"))

        mapped_gene_ids <<- c(mapped_gene_ids, unique(brainarray_df$Probe.Set.Name))
    }))

    data.frame(mapped_gene_ids) %>%
        dplyr::count(mapped_gene_ids) %>%
        write.table(paste0("/out/01_scraping_brainarray_genes/", species, ".tsv"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}
