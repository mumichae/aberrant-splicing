#'---
#' title: Create Muscle dataset combined with GTEx
#' author: Christian Mertes
#' wb:
#'  input:
#'    - html_out: '`sm config["htmlOutputPath"] + "/FraseR/MUSCLE_GTEx_results.html"`'
#'    - fdsobj:   '`sm config["PROC_DATA"] + "/datasets/savedObjects/MUSCLE_GTEx/fds-object.RDS"`'
#'    - resFile:  '`sm config["PROC_DATA"] + "/processedData/results/MUSCLE_GTEx_results.tsv"`'
#'    - anno:     "/s/project/crg_seq_data/raw_data/URDCAT_sample_annotation.csv"
#'    - vars:     "/s/project/crg_seq_data/processed_results/process_vcf/variant_dt.Rds"
#'  threads: 10
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    parseWBHeader2("./Scripts/Analysis/Muscle_with_GTEx.R")
}

#+ load main config, echo=FALSE
source("./src/r/config.R", echo=FALSE)
library(gridExtra)

#+ input
fdsFile  <- snakemake@input$fdsobj
annoFile <- snakemake@input$anno
name     <- basename(dirname(fdsFile))
wdir     <- dirname(dirname(dirname(fdsFile)))
fds      <- loadFraseRDataSet(wdir, name)
resFile  <- snakemake@input$resFile

#'
#' Load annotation
#'
anno <- fread(annoFile)[RNA_fi %in% samples(fds)]
DT::datatable(anno)


#'
#' # Clustering with GTEx
#'
colData(fds)$Patient <- ifelse(grepl("GTEX-", samples(fds)), "GTEx", "Patient")
for(type in psiTypes){
    p1 <- plotCountCorHeatmap(fds, type, logit=TRUE, normalized=FALSE, topN=30000, annotation_col="Patient")
    p2 <- plotCountCorHeatmap(fds, type, logit=TRUE, normalized=TRUE,  topN=30000, annotation_col="Patient")
    grid.arrange(arrangeGrob(grobs=list(p1[[4]], p2[[4]]),ncol=2))
}


#'
#' # Load results
#'
res <- fread(resFile)
res <- res[signif==TRUE & abs(deltaPsi) > 0.3]
plotNBySample(res[type != "psiSite"], samples(fds),
        main="Aberrant Events by sample\nOnly alternative splicing")
plotNBySample(res, samples(fds),
        main="Aberrant Events by sample\nIncluding intron retention")

#'
#' ## Load variants
#'
vars <- readRDS(snakemake@input$vars)[!is.na(hgncid)]


#' ## Cleanup results
res[,bamFile:=NULL]
res[,condition:=NULL]
res <- res[!is.na(hgnc_symbol) & !grepl("GTEX-", sampleID)]
res <- res[order(-abs(deltaPsi))]

#'
#' # Results
#'
DT::datatable(res)

#'
#' ## Cherry Picking
#'
cherries <- data.table(do.call(rbind, list(
    c(hgncid=NA, sampleID=NA, resIdx=NA),
    c("ABCB7",   "EPR283656.allchr.genotyped", 1),
    c("BLCAP",   "EPR117314.allchr.genotyped", 1),
    c("ARFGAP2", "EPR117314.allchr.genotyped", 1),
    c("FEZ2",    "EPR283656.allchr.genotyped", 1),
    c("ELMOD3",  "EPR117314.allchr.genotyped", 1)
        )))[!is.na(hgncid)]

DT::datatable(cherries)

for(i in nrow(cherries)){
    gene <- cherries[i, hgncid]
    sa   <- cherries[i, sampleID]
    idx  <- as.integer(cherries[i, resIdx])

    res[hgnc_symbol == gene]
    vars[hgncid == gene & sample == sa]
    plotJunctionDistribution(fds, gr=res[hgnc_symbol == gene][idx], plotLegend=FALSE)
}

