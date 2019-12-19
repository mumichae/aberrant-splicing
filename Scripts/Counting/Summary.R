#'---
#' title: "Count Summary: `r gsub('_', ' ', snakemake@wildcards$dataset)`"
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 1
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - filter: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/filter.done" `'
#'  output:
#'   - wBhtml: '`sm config["htmlOutputPath"] + 
#'                  "/aberrant_splicing/FraseR/{dataset}_countSummary.html"`'
#'  type: noindex
#'---

#+ echo=FALSE
source("Scripts/_helpers/config.R")
opts_chunk$set(fig.width=12, fig.height=8)


#+ input
dataset    <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir <- snakemake@params$workingDir
bpWorkers   <- min(max(extract_params(bpworkers()), 1),
                   as.integer(extract_params(snakemake@params$workers)))


fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))

#' Number of samples: `r nrow(colData(fds))`
#' Number of introns (psi5): `r nrow(rowRanges(fds, type = "psi5"))`
#' Number of introns (psi3): `r nrow(rowRanges(fds, type = "psi3"))`
#' Number of splice sites (psiSite): `r nrow(rowRanges(fds, type = "psiSite"))`

#' Introns that passed filter
table(mcols(fds, type="j")[,"passed"])

#' ## Expression filtering
plotFilterExpression(fds) + theme_cowplot(font_size = 16)

#' ## Correlation between samples
plots <- lapply(FraseR::psiTypes, plotCountCorHeatmap, fds=fds, logit=TRUE, 
                topN=100000)
