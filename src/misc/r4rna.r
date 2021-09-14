###############################################################################
# R4RNA - Sample code
###############################################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require(R4RNA)) {
  BiocManager::install("R4RNA") # , version = "3.8")
  library(R4RNA)
}

if (!require(Biostrings)) {
  BiocManager::install("Biostrings") # , version = "3.8")
  library(Biostrings)
}

known_file <-
  system.file("extdata", "vienna.txt", package = "R4RNA")
known <- readVienna(known_file)
known <- expandHelix(known)


message("Multiple sequence alignment of interest")
fasta_file <- system.file("extdata", "fasta.txt", package = "R4RNA")
fasta <- as.character(readBStringSet(fasta_file))
message("Plot covariance in alignment")
plotCovariance(fasta, known, cex = 0.5)
