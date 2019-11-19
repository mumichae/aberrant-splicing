#'---
#' title: Filter and clean dataset
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 1
#'  input:
#'   - psiSS:  '`sm parser.getProcDataDir()+ "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/psiSite.h5"`'
#'   - dPsiSS: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/delta_psiSite.h5"`'
#'  output:
#'   - dPsiSS: '`sm parser.getProcDataDir()+ "/aberrant_splicing/datasets/savedObjects/{dataset}/fds-object.RDS"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/aberrant_splicing/FraseR/{dataset}_filterExpression.html"`'
#'  type: noindex
#'---

if(FALSE){
    snakemake <- readRDS("./tmp/snakemake.RDS")
    source(".wBuild/wBuildParser.R")
    parseWBHeader("./Scripts/FraseR/03_filter_expression_FraseR.R", dataset="Kremer")
}

#+ echo=FALSE
source("./src/r/config.R")
opts_chunk$set(fig.width=12, fig.height=8)

#+ input
dataset    <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir <- dirname(dirname(dirname(snakemake@output$countsJ)))
bpWorkers   <- min(max(extract_params(bpworkers()), 1),
                   as.integer(extract_params(snakemake@params$workers)))

#'
#' # Load count data
#+ echo=TRUE
dataset

#+ echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))
dim(fds)


#'
#' Filter FraseR object based on standard values
#'
fds <- filterExpression(fds, filter=FALSE)
table(mcols(fds, type="j")[,"passed"])

#+ filter bad junctions, echo=FALSE
devNull <- saveFraseRDataSet(fds)
name(fds) <- dataset
fds <- saveFraseRDataSet(fds[mcols(fds, type="j")[,"passed"]])

#'
#' Correlation of counts after filtering
#'
plots <- lapply(psiTypes, plotCountCorHeatmap, fds=fds, logit=TRUE, topN=100000)

