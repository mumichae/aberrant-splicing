#'---
#' title: Calculate PSI values
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 10
#'   - threads: 10
#'   - progress: FALSE
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - countsJ:  '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsSS: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'  output:
#'  - psiSS:     '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/psiSite.h5"`'
#'  - dPsiSS:    '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/delta_psiSite.h5"`'
#'  - wBhtml:    '`sm config["htmlOutputPath"] + "/FraseR/{dataset}_psi_value_calculation.html"`'
#'  type: noindex
#'---

if(FALSE){
    snakemake <- readRDS("./tmp/snakemake.RDS")
    source(".wBuild/wBuildParser.R")
    parseWBHeader("./Scripts/FraseR/03_filter_expression_FraseR.R", dataset="example")
}

#+ echo=FALSE
source("./src/r/config.R")

#+ input
dataset    <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir <- snakemake@params$workingDir
bpWorkers   <- min(max(extract_params(bpworkers()), 1),
                   as.integer(extract_params(snakemake@params$workers)))
bpThreads   <- as.integer(extract_params(snakemake@params$threads))
bpProgress  <- as.logical(extract_params(snakemake@params$progress))

#'
#' # Load count data
#+ echo=TRUE
dataset

#+ echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))
register(MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress))

#'
#' Calculating PSI values
#'
fds <- calculatePSIValues(fds)

#'
#' FraseR object after PSI value calculation
#'
fds <- saveFraseRDataSet(fds)
fds

