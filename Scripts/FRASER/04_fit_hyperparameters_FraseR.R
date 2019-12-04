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
#'   - filter: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/filter.done" `'
#'  output:
#'   - hyper: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/hyper.done" `'
#'  type: script
#'---

source("./src/r/config.R")

#+ input
dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
bpWorkers   <- min(max(extract_params(bpworkers()), 1),
                   as.integer(extract_params(snakemake@params$workers)))
bpThreads   <- as.integer(extract_params(snakemake@params$threads))
bpProgress  <- as.logical(extract_params(snakemake@params$progress))
register(MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress))

# Load PSI data
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)

# Run hyper parameter optimization
correction <- snakemake@config$aberrantSplicing$correction
for(type in FraseR::psiTypes){
    message(date(), ": ", type)
    fds <- optimHyperParams(fds, type=type, correction=correction)
    fds <- saveFraseRDataSet(fds)
}
fds <- saveFraseRDataSet(fds)

