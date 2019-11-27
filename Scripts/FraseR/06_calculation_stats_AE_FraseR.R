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
#'   - fdsin:  '`sm parser.getProcDataDir()+ "/aberrant_splicing/datasets/savedObjects/{dataset}/predictedMeans_psiSite.h5"`'
#'  output:
#'   - fdsout: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/{dataset}/pajdBetaBinomial_psiSite.h5"`'
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
workingDir <- snakemake@params$workingDir
bpWorkers   <- min(max(extract_params(bpworkers()), 1),
                   as.integer(extract_params(snakemake@params$workers)))
bpThreads   <- as.integer(extract_params(snakemake@params$threads))
bpProgress  <- as.logical(extract_params(snakemake@params$progress))

#'
#' # Load Zscores data
#+ echo=TRUE
dataset
workingDir

#+ echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)
register(MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress))


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

