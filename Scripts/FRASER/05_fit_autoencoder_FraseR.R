#'---
#' title: Fitting the autoencoder
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 20
#'   - threads: 20
#'   - progress: FALSE
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - hyper: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/hyper.done" `'
#'  output:
#'   - fdsout: '`sm parser.getProcDataDir() + 
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/predictedMeans_psiSite.h5"`'
#'  type: script
#'---

source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
bpWorkers  <- min(max(extract_params(bpworkers()), 1),
                  as.integer(extract_params(snakemake@params$workers)))
bpThreads  <- as.integer(extract_params(snakemake@params$threads))
bpProgress <- as.logical(extract_params(snakemake@params$progress))
register(MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress))

# Load PSI data
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)

# Fit autoencoder
# run it for every type
correction <- snakemake@config$aberrantSplicing$correction

for(type in psiTypes){

    currentType(fds) <- type
    q <- bestQ(fds, type)
    fds <- fit(fds, q=q, type=type, verbose=TRUE, iterations=15, 
               correction=correction)

    fds <- saveFraseRDataSet(fds)
}

