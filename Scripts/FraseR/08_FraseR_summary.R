#'---
#' title: "FraseR Summary: `r gsub('_', ' ', snakemake@wildcards$dataset)`"
#' author: mumichae, vyepez, ischeller
#' wb:
#'  input:
#'   - fdsin: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/{dataset}/pajdBetaBinomial_psiSite.h5"`'
#'   - results: '`sm parser.getProcDataDir() + "/aberrant_splicing/results/{dataset}_results.tsv"`'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/aberrant_splicing/FraseR/{dataset}_summary.html"`'
#'  type: noindex
#'---

if(FALSE){
    snakemake <- readRDS("./tmp/snakemake.RDS")
}

#+ input
dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- dirname(dirname(dirname(fdsFile)))

#+ load config and setup, echo=FALSE
source("./src/r/config.R")

#'
#' # Load data
#'
#' Dataset:
#+ echo=TRUE
dataset
workingDir

#+ echo=FALSE
fds_raw <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)
#' Number of samples: `r ncol(fds)`
#' Number of junctions: `r nrow(fds)`


# used for most plots
dataset_title <- paste("Dataset:", snakemake@wildcards$dataset)

#'
#' ## Visualize
#' ### Filter expression
plotFilterExpression(fds_raw, bins=100)

#' ### Hyper parameter optimization
for(type in psiTypes){
    print(plotEncDimSearch(fds, type=type))
}

#' ### Aberrant genes per sample
for(type in psiTypes){
    plotAberrantPerSample(fds, type=type, aggregate=TRUE, main=paste(dataset_title, " (", type, ")"))
}

#' ### Batch Correlation
topN <- 30000
topJ <- 10000
#' #### Samples x samples
for(type in psiTypes){
    before <- plotCountCorHeatmap(
        fds = fds,
        type = type,
        logit = TRUE,
        topN = topN,
        topJ = topJ,
        plotType = "sampleCorrelation",
        normalized = FALSE,
        annotation_col = NA,
        annotation_row = NA,
        sampleCluster = NA,
        plotMeanPsi=FALSE,
        plotCov = FALSE,
        annotation_legend = TRUE
    )
    before
    after <- plotCountCorHeatmap(
        fds = fds,
        type = type,
        logit = TRUE,
        topN = topN,
        topJ = topJ,
        plotType = "sampleCorrelation",
        normalized = TRUE,
        annotation_col = NA,
        annotation_row = NA,
        sampleCluster = NA,
        plotMeanPsi=FALSE,
        plotCov = FALSE,
        annotation_legend = TRUE
    )
    after
}
#' #### Junctions x samples
for(type in psiTypes){
    before <- plotCountCorHeatmap(
        fds = fds,
        type = type,
        logit = TRUE,
        topN = topN,
        topJ = topJ,
        plotType = "junctionSample",
        normalized = FALSE,
        annotation_col = NA,
        annotation_row = NA,
        sampleCluster = NA,
        plotMeanPsi=FALSE,
        plotCov = FALSE,
        annotation_legend = TRUE
    )
    before
    after <- plotCountCorHeatmap(
        fds = fds,
        type = type,
        logit = TRUE,
        topN = topN,
        topJ = topJ,
        plotType = "junctionSample",
        normalized = TRUE,
        annotation_col = NA,
        annotation_row = NA,
        sampleCluster = NA,
        plotMeanPsi=FALSE,
        plotCov = FALSE,
        annotation_legend = TRUE
    )
    after
}


