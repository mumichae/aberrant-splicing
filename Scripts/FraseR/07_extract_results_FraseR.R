#'---
#' title: Results of FraseR analysis
#' author: Christian Mertes
#' wb:
#'  input:
#'   - fdsin: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/{dataset}/pajdBinomial_psiSite.h5"`'
#'  output:
#'   - resultTable: '`sm parser.getProcDataDir() + "/aberrant_splicing/results/{dataset}_results.tsv"`'
#'   - wBhtml: '`sm parser.getProcDataDir() + "/aberrant_splicing/FraseR/{dataset}_results.html"`'
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

# 
# TODO: how much do we do here for a standard analysis?
#       And how much do we leave it up to the user?
# 
# R.utils::sourceDirectory("../gagneurlab_shared/r/disease")
# R.utils::sourceDirectory("../gagneurlab_shared/r/go_enrichment")
# MGSA_GO_FULL <- load_mgsaset_for_organism('human')
# mimTable <- getFullOmimTable()
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
bpparam <- MulticoreParam(3, 3)
parallel(fds) <- bpparam

#'
#' ## Extract results
#' You can adjust the cutoffs to your needs
#' the current defaults are:
#' padj   <= 0.1
#' dpsi   <= abs(0.1)
#' N      >= 10
#' zScore >= 0
#' 
resgr <- results(fds, zscoreCut=0)
res   <- as.data.table(resgr)
saveFraseRDataSet(fds)

#'
#' * Add features
#'     * number of samples per gene and variant
res[p.adj<=0.1, numSamplesPerGene:=length(unique(sampleID)), by=hgnc_symbol]
res[p.adj<=0.1, numEventsPerGene:=.N, by="hgnc_symbol,sampleID"]
res[p.adj<=0.1, numSamplesPerJunc:=length(unique(sampleID)), by="seqnames,start,end"]

#'
#'     * MitoVIP genes
# vip_genes <- get_vip_info_table()[,.(hgnc_symbol=gene,
#         isMitoVIP=ifelse(causal,"causal", "vip"))]
# res <- merge(res, vip_genes, all.x=TRUE, "hgnc_symbol")

#'
#'     * OMIM phenotypes
#'
# res <- merge(res, mimTable, all.x=TRUE, by.x="hgnc_symbol", by.y="SYMBOL")

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

# get links
# res[,genecards:=get_html_link(hgnc_symbol, website="genecards", TRUE)]
res[,hgnc:=get_html_link(hgnc_symbol, website="hgnc", TRUE)]
# res[,omim:=get_html_link(GMIM, website="omim", TRUE)]
res[,entrez:=get_html_link(hgnc_symbol, website="entrez", TRUE)]
res[,locus:=get_html_link(paste0(seqnames, ":", start, "-", end), website="locus", TRUE)]

# round numbers
res[,p.adj:=signif(p.adj, 3)]
res[,deltaPsi:=signif(deltaPsi, 2)]
res[,zscore:=signif(zscore, 2)]
res[,psiValue:=signif(psiValue, 2)]

# set correct order
colOrders <- unique(c("sampleID", "genecards", "p.adj", "deltaPsi", "type",
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

