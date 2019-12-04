#'---
#' title: Calculate P values
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 20
#'   - threads: 20
#'   - progress: FALSE
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - fdsin:  '`sm parser.getProcDataDir() + 
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/" +
#'                  "predictedMeans_psiSite.h5"`'
#'  output:
#'   - fdsout: '`sm parser.getProcDataDir() + 
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/" +
#'                  "pajdBetaBinomial_psiSite.h5"`'
#'  type: script
#'---

source("./src/r/config.R")

dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir
bpWorkers   <- min(max(extract_params(bpworkers()), 1),
                   as.integer(extract_params(snakemake@params$workers)))
bpThreads   <- as.integer(extract_params(snakemake@params$threads))
bpProgress  <- as.logical(extract_params(snakemake@params$progress))
register(MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress))

# Load Zscores data
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)

# Calculate stats
for (type in psiTypes) {
    # Zscores
    fds <- calculateZscore(fds, type=type)
    # Pvalues
    fds <- calculatePvalues(fds, type=type)
    # Adjust Pvalues
    fds <- calculatePadjValues(fds, type=type)
}

fds <- saveFraseRDataSet(fds)

