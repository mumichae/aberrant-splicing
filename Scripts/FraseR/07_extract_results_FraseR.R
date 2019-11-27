#'---
#' title: Results of FraseR analysis
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - fdsin: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/{dataset}/pajdBetaBinomial_psiSite.h5"`'
#'  output:
#'   - resultTable: '`sm parser.getProcDataDir() + "/aberrant_splicing/results/{dataset}_results.tsv"`'
#'   - wBhtml: '`sm parser.getProcDataDir() + "/aberrant_splicing/FraseR/{dataset}_results.html"`'
#'  type: noindex
#'---

#+ input
dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir

#+ load config and setup, echo=FALSE
source("./src/r/config.R")

#
# TODO: how much do we do here for a standard analysis?
#       And how much do we leave it up to the user?
#
opts_chunk$set(fig.width=12, fig.height=8)


#'
#' # Load data
#'
#' Dataset:
#+ echo=TRUE
dataset
workingDir

#+ echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)
register(MulticoreParam(3, 3))

#'
#' ## Extract results
#'
config_params <- snakemake@config$aberrantSplicing
resgr <- results(fds,
                 padjCutoff=config_params$padjCutoff, 
                 zScoreCutoff=config_params$zScoreCutoff,
                 deltaPsiCutoff=config_params$deltaPsiCutoff)
res   <- as.data.table(resgr)
saveFraseRDataSet(fds)

#'
#' * Add features
#'     * number of samples per gene and variant
res[padjust<=0.1, numSamplesPerGene:=length(unique(sampleID)), by=hgncSymbol]
res[padjust<=0.1, numEventsPerGene:=.N, by="hgncSymbol,sampleID"]
res[padjust<=0.1, numSamplesPerJunc:=length(unique(sampleID)), by="seqnames,start,end"]

#'
#'     * add colData at the end
res <- merge(res, as.data.table(colData(fds)), by="sampleID")

#'
#' # Results
#'
write_tsv(res, file=snakemake@output$resultTable)
file <- gsub(".html$", ".tsv", snakemake@output$wBhtml)
write_tsv(res, file=file)

#'
#' The result table can also be downloaded with the link below.
#'
#+ echo=FALSE, results='asis'
cat(paste0("<a href='./", basename(file), "'>Download result table</a>"))

# round numbers
res[,padjust:=signif(padjust, 3)]
res[,deltaPsi:=signif(deltaPsi, 2)]
res[,zscore:=signif(zScore, 2)]
res[,psiValue:=signif(psiValue, 2)]

# set correct order
colOrders <- unique(c("sampleID", "genecards", "padjust", "deltaPsi", "type",
                      "numSamplesPerGene", "numEventsPerGene", "numSamplesPerJunc",
                      "isMitoVIP", "omim", "PMIM", "PINH", "locus", "hgnc", colnames(res)))
colOrders <- colOrders[colOrders %in% colnames(res)]
setcolorder(res, colOrders)

#'
#' * Result table
DT::datatable(res, options=list(scrollX=TRUE), escape=FALSE)

#'
#' * Sample table
DT::datatable(as.data.table(colData(fds)), options=list(scrollX=TRUE))

#'
#' * Sample correlation
plots <- lapply(psiTypes, plotCountCorHeatmap, fds=fds, logit=TRUE, topN=100000, norm=FALSE)
plots <- lapply(psiTypes, plotCountCorHeatmap, fds=fds, logit=TRUE, topN=100000, norm=TRUE)

