

###############################################################################
# NCBI API - Descargar secuencias de GenBank
###############################################################################


###########################################################
# 1. Libraries required
###########################################################

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!require(annotate)) {
  BiocManager::install("annotate")
  library(annotate)
}

if (!require(dplyr)) {
  install.packages('dplyr')
  library(dplyr)
}

if (!require(tidyverse)) {
  install.packages('tidyverse')
  library(tidyverse)
}

if (!require(stringr)) {
  install.packages('stringr')
  library(stringr)
}

###########################################################
# 2. Functions and datatypes
###########################################################

# https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly

refSeq <- c(
  "AC_",
  "NC_",
  "NG_",
  "NT_",
  "NW_",
  "NZ_b",
  "NM_",
  "NR_",
  "XM_c",
  "XR_c",
  "AP_",
  "NP_",
  "YP_c",
  "XP_c",
  "WP_"
)

names(refSeq) <- c(
  "Genomic Complete genomic molecule, alternate assembly",
  "Genomic Complete genomic molecule, reference assembly",
  "Genomic Incomplete genomic region",
  "Genomic Contig or scaffold, clone-based or WGSa",
  "Genomic Contig or scaffold, primarily WGSa",
  "Genomic Complete genomes and unfinished WGS data",
  "mRNA Protein-coding transcripts (usually curated)",
  "RNA Non-protein-coding transcripts",
  "mRNA Predicted model protein-coding transcript",
  "RNA Predicted model non-protein-coding transcript",
  "Protein Annotated on AC_ alternate assembly",
  "Protein Associated with an NM_ or NC_ accession",
  "Protein Annotated on genomic molecules without an instantiated transcript record",
  "Protein Predicted model, associated with an XM_ accession",
  "Protein Non-redundant across multiple strains and species"
)


# Get transcripts list in NCBI RefSeq code
###########################################################
showError <- function(err, type = "Error:") {
  print(paste(type, err$message))
  print(paste("Call:", err$call))
}

getNCBIRefSeq <- function(lst, prefix, prmType = "accession") {
  # Clear screen
  cat("\014")
  
  print("Starting process...")
  
  for (i in 1:nrow(lst)) {
    #print("Calling GenBank")
    
    tryCatch(
      xdoc <- genbank(lst$GeneID[i], type = prmType, disp = "data"),
      warning = function(err) {
        showError(err, "Warning")
        print(xdoc)
      },
      error = function(err) {
        showError(err)
        print(xdoc)
      }
    )
    
    
    # Check that id exists and check integrity of xdoc
    if (!is.null(xdoc) && class(xdoc) != "logical") {
      xmlDoc <- xmlToList(xdoc)
      xmlDoc <- unlist(xmlDoc)
      ix <- str_detect(xmlDoc, prefix)
      aux <- xmlDoc[ix]
      aux <- unlist(aux)
      
      
      if (i == 1) {
        res <- aux %>% unique()
      } else {
        aux <- aux %>% unique()
        print(paste("Fila:",i,"Cont:",aux))
        res <- c(res, aux)
      }
    }
    # pause between API requests
    Sys.sleep(8)
  }
  
  print("Finished.")
  res <- res %>% unique()
  return(res)
  
}

# [1] "Fila: 166 Cont: NM_164834"    "Fila: 166 Cont: NM_057539"   
# [3] "Fila: 166 Cont: NM_001103656" "Fila: 166 Cont: NM_001298829"
# [5] "Fila: 166 Cont: NM_164833" 
###########################################################
# 3. Use example
###########################################################

lst <- c("", "", "")
lst <- read.table("./data/fly_repressed.txt", header = T)

# Get transcripts (NM_x code)
res <- getNCBIRefSeq(lst, prefix = refSeq[7])

print(res)
write(res, "./data/fly_repressed.lst")
