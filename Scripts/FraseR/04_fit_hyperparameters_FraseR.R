#'---
#' title: Hyper parameter optimization
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 10
#'   - threads: 10
#'   - progress: FALSE
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
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
workingDir <- snakemake@params$workingDir
bpWorkers   <- min(max(extract_params(bpworkers()), 1),
                   as.integer(extract_params(snakemake@params$workers)))
bpThreads   <- as.integer(extract_params(snakemake@params$threads))
bpProgress  <- as.logical(extract_params(snakemake@params$progress))

#'
#' # Load PSI data
#+ echo=TRUE
dataset

#+ echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)
register(MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress))
dim(fds)

#'
#' # Run hyper parameter optimization
#'
correction <- snakemake@config$aberrantSplicing$correction
for(type in psiTypes){
    message(date(), ": ", type)
    fds <- optimHyperParams(fds, type=type, correction=correction)
    fds <- saveFraseRDataSet(fds)
}

#'
#' FraseR object
#'
fds <- saveFraseRDataSet(fds)
fds

