# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("R4RNA", version = "3.8")
# 
# library(R4RNA)
# library(Biostrings)
# 
# orig <- laslo[laslo$Serie=="Original", ]
# 
# for(i in 1:nrow(orig)){
#   line <- paste(">",i, laslo[i,]$Gen)
#   write(line,file="fasta.txt",append=TRUE)
#   write(laslo[i,]$StemLoopSequence ,file="fasta.txt",append=TRUE)
# }
# 
# known_file <- system.file("extdata", "vienna.txt", package = "R4RNA")
# known <- readVienna(known_file)
# known <- expandHelix(known)
# 
# 
# message("Multiple sequence alignment of interest")
# fasta_file <- "fasta.txt" #system.file("extdata", "fasta.txt", package = "R4RNA")
# fasta <- as.character(readBStringSet(fasta_file))
# message("Plot covariance in alignment")
# plotCovariance(fasta, known, cex = 0.5)

###########################################
# glm

laslo_model <- glm(Serie  ~ RNAFoldMFE + LoopPattern + N.1 + N2 
                      + Bulges + InternalLoops + StemLength, 
                      data = laslo, family = "binomial")
