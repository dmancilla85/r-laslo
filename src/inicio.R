# Cargar framework
source("./src/configuration.R")

######

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("annotate")

library(annotate)

xdoc <- genbank("FBgn0023537", type="accession",disp="data")
class(xdoc)

lst <- xmlToList(xdoc)
lst <- unlist(lst)

ix <- str_detect(lst, "NM_")
res <- lst[ix]
res %>% unique()

######