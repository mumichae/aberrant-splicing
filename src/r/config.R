##--------------------------------------------
## required packages
message("Load needed packages")
suppressPackageStartupMessages({
    library(R.utils)
    library(markdown)
    library(knitr)
    library(devtools)
    library(yaml)
    library(BBmisc)
    library(GenomicAlignments)
    library(tidyr)
    library(plotly)
    library(DelayedMatrixStats)
    #library(FraseR)
})

# load the FraseR package with devtools
suppressPackageStartupMessages({
  devtools::load_all("../FraseR")
})


# source function directory for extra functionality
sourceDirectory("./src/r/functions", pattern=".*\\.R")

## helper functions
write_tsv <- function(x, file, row.names = FALSE, ...){
  write.table(x=x, file=file, quote=FALSE, sep='\t', row.names= row.names, ...)
}

