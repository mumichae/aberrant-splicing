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
#'   - countsJ:  '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsSS: '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'  output:
#'  - psiSS:     '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/psiSite.h5"`'
#'  type: script
#'---

source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir <- snakemake@params$workingDir
bpWorkers   <- min(max(extract_params(bpworkers()), 1),
                   as.integer(extract_params(snakemake@params$workers)))
bpThreads   <- as.integer(extract_params(snakemake@params$threads))
bpProgress  <- as.logical(extract_params(snakemake@params$progress))


fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))
register(MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress))

# Calculating PSI values
fds <- calculatePSIValues(fds)

# FraseR object after PSI value calculation
fds <- saveFraseRDataSet(fds)
