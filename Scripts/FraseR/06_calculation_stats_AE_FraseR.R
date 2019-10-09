#'---
#' title: Calculate P values
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 20
#'   - threads: 20
#'   - progress: FALSE
#'  input:
#'   - fdsin:  '`sm parser.getProcDataDir()+ "/aberrant_splicing/datasets/savedObjects/{dataset}/predictedMeans_psiSite.h5"`'
#'  output:
#'   - fdsout: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/{dataset}/pajdBinomial_psiSite.h5"`'
#'   - wBhtml: '`sm parser.getProcDataDir() + "/aberrant_splicing/FraseR/{dataset}_stat_calculation.html"`'
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
fdsFile    <- snakemake@input$fdsin
workingDir <- dirname(dirname(dirname(fdsFile)))
bpWorkers  <- min(bpworkers(), as.integer(snakemake@params$workers))
bpThreads  <- min(bpworkers(), as.integer(snakemake@params$threads))
bpProgress <- snakemake@params$progress

#'
#' # Load Zscores data
#+ echo=TRUE
dataset
workingDir

#+ echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)
bpparam <- MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress)
parallel(fds) <- bpparam


#'
#' # Calculate stats
#'
#' ## Zscores
#
for(type in psiTypes){
    fds <- calculateZscore(fds, type=type)
}

#'
#' ## Pvalues
for(type in psiTypes){
    fds <- calculatePvalues(fds, type=type)
}

#'
#' ## Adjust Pvalues
for(type in psiTypes){
    fds <- calculatePadjValues(fds, type=type)
}

#'
#' ## Annotate ranges
fds <- annotateRanges(fds)


#'
#' # Save results
#'
fds <- saveFraseRDataSet(fds)
fds

