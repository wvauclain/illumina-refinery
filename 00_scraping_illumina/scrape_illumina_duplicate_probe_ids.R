suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(lazyeval))
suppressPackageStartupMessages(library(dplyr))


platforms <- c(
  'illuminaHumanv1',
  'illuminaHumanv2',
  'illuminaHumanv3',
  'illuminaHumanv4',
  'illuminaMousev1',
  'illuminaMousev1p1',
  'illuminaMousev2',
  'illuminaRatv1'
)

for (platform in platforms) {
  write(paste0("Scraping probe ids for ", platform, "..."), file=stdout())
  suppressPackageStartupMessages(library(paste0(platform, ".db"), character.only=TRUE))

  # Extract platform-specific data
  # See here for more discussion:
  # https://github.com/AlexsLemonade/refinebio/pull/212#discussion_r184864928
  probeQualityRef <- lazy_eval(paste0(platform, ".db::", platform, "PROBEQUALITY"))
  probeQuality <- as.data.frame(probeQualityRef[mappedkeys(probeQualityRef)])

  probeGeneRef <-  lazy_eval(paste0(platform, ".db::", platform, "ENSEMBL"))
  probeGene <- as.data.frame(probeGeneRef[mappedkeys(probeGeneRef)])

  # Let's limit our focus to good and perfect probes, since that's all we need for refine.bio
  goodProbeIndices <- which(grepl("Good", probeQuality$ProbeQuality))
  perfectProbeIndices <- which(grepl("Perfect", probeQuality$ProbeQuality))
  probesToKeep <- probeQuality[c(goodProbeIndices, perfectProbeIndices),1]

  isMappedToMultipleGenes <- sapply(probesToKeep, function(probe) {
    return(sum(probeGene$probe_id == probe) > 1)
  })

  probesToKeep <- probesToKeep[isMappedToMultipleGenes]

  # Map probes to genes
  probeGene <- probeGene[probeGene$probe_id %in% probesToKeep,]

  # Save to output file
  write.table(probeGene, paste0("/out/00_scraping_illumina/", platform, ".tsv"), row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
}


