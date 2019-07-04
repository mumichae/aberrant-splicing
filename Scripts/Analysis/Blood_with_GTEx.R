#'---
#' title: Analysis of the Blood/GTEx samples
#' author: Christian Mertes
#' wb:
#'  input:
#'    - html_out: '`sm config["htmlOutputPath"] + "/FraseR/BLOOD_GTEx_results.html"`'
#'    - fdsobj:   '`sm config["PROC_DATA"] + "/datasets/savedObjects/BLOOD_GTEx/fds-object.RDS"`'
#'    - resFile:  '`sm config["PROC_DATA"] + "/processedData/results/BLOOD_GTEx_results.tsv"`'
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
  parseWBHeader2("./Scripts/Analysis/Blood_with_GTEx.R")
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
res <- res[numSamplesPerGene < 4]
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
  c(hgncid=NA, sampleID=NA, resIdx=NA, comment=NA),
  c("ALDH3A2", "EPR513463.allchr.genotyped", 1, "Align error"),
  c("SELPLG",  "EPR697596.allchr.genotyped", 1, ""),
  c("HCFC1R1", "EPR133443.allchr.genotyped", 1, ""),
  c("LILRB2",  "EPR133443.allchr.genotyped", 1, ""),
  c("COMMD1",  "EPR133443.allchr.genotyped", 1, ""),
  c("ACO1",    "RNA_AF6412.1868AB", 1, ""),
  c("SHKBP1",  "RNA_AF6412.1868AB", 1, ""),
  c("CISD1",   "EPR697596.allchr.genotyped", 1, "")
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

