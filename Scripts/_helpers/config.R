##--------------------------------------------
## required packages
message("Load packages")
suppressPackageStartupMessages({
    library(markdown)
    library(knitr)
    library(devtools)
    library(yaml)
    library(BBmisc)
    library(GenomicAlignments)
    library(tidyr)
    library(plotly)
    library(DelayedMatrixStats)
    library(FRASER)
})

# load the FraseR package with devtools
# suppressPackageStartupMessages({
#     devtools::load_all("/data/ouga/home/ag_gagneur/schuller/MLL_Thesis_LS/workspace/FraseR")
# })

## helper functions
write_tsv <- function(x, file, row.names = FALSE, ...){
  write.table(x=x, file=file, quote=FALSE, sep='\t', row.names= row.names, ...)
}

extract_params <- function(params) {
    unlist(params)[1]
}

