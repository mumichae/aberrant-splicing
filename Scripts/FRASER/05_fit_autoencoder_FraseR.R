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

source("./src/r/config.R")

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

    # set current type
    currentType(fds) <- type
    curDims <- dim(K(fds, type))
    
    q <- bestQ(fds, type)
    probE <- max(0.001, min(1,30000/curDims[1]))

    # subset fitting
    featureExclusionMask(fds) <- sample(c(TRUE, FALSE), curDims[1],
            replace=TRUE, prob=c(probE, 1-probE))

    # run autoencoder
    fds <- fit(fds, q=q, type=type, verbose=TRUE, iterations=15, 
               correction=correction)

    # save autoencoder fit
    fds <- saveFraseRDataSet(fds)
}

