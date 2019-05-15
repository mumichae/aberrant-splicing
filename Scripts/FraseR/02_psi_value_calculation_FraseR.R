#'---
#' title: Calculate PSI values
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 10
#'   - threads: 10
#'   - progress: FALSE
#'  input:
#'   - countsJ:  '`sm config["PROC_DATA"] + "/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsSS: '`sm config["PROC_DATA"] + "/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'  output:
#'  - psiSS:     '`sm config["PROC_DATA"] + "/datasets/savedObjects/raw-{dataset}/psiSite.h5"`'
#'  - dPsiSS:    '`sm config["PROC_DATA"] + "/datasets/savedObjects/raw-{dataset}/delta_psiSite.h5"`'
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
dataset     <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir  <- dirname(dirname(dirname(snakemake@input$countsJ)))
bpWorkers   <- min(bpworkers(), as.integer(snakemake@params[[1]]$workers))
bpThreads   <- as.integer(snakemake@params[[1]]$threads)
bpProgress  <- as.logical(snakemake@params[[1]]$progress)

#'
#' # Load count data
#+ echo=TRUE
dataset

#+ echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))
parallel(fds) <- MulticoreParam(bpWorkers, bpThreads,
        progressbar=bpProgress)

#'
#' Calculating PSI values
#'
fds <- calculatePSIValues(fds)

#'
#' FraseR object after PSI value calculation
#'
fds <- saveFraseRDataSet(fds)
fds

