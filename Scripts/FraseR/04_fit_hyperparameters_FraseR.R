#'---
#' title: Hyper parameter optimization
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 10
#'   - threads: 10
#'   - progress: FALSE
#'  input:
#'   - inFile: '`sm config["htmlOutputPath"] + "/aberrant_splicing/FraseR/{dataset}_filterExpression.html"`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/aberrant_splicing/FraseR/{dataset}_hyper_parameter_optimization.html"`'
#'  type: noindex
#'---
##
## TODO:
##   Add a link to the fraser object to have a proper chain of events in wbuild
##

if(FALSE){
    snakemake <- readRDS("./tmp/snakemake.RDS")
    source(".wBuild/wBuildParser.R")
    parseWBHeader("./Scripts/FraseR/04_fit_hyperparameters_FraseR.R", dataset="example")
}

#+ echo=FALSE
source("./src/r/config.R")

#+ input
dataset    <- snakemake@wildcards$dataset
workingDir <- file.path(snakemake@config$root, "processed_data", 
    "aberrant_splicing", "datasets")
bpWorkers  <- min(bpworkers(), as.integer(snakemake@params$workers))
bpThreads  <- min(bpworkers(), as.integer(snakemake@params$threads))
bpProgress <- snakemake@params$progress


#'
#' # Load PSI data
#+ echo=TRUE
dataset

#+ echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)
parallel(fds) <- MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress)
dim(fds)

#'
#' # Run hyper parameter optimization
#'
for(type in psiTypes){
    message(date(), ": ", type)
    fds <- optimHyperParams(fds, type=type, BPPARAM=parallel(fds))
    fds <- saveFraseRDataSet(fds)
}

#'
#' FraseR object
#'
fds <- saveFraseRDataSet(fds)
fds

