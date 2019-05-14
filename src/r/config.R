##--------------------------------------------
## required packages
message("Load needed packages")
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
    #library(FraseR)
})

# load the FraseR package with devtools
devtools::load_all("../FraseR")
