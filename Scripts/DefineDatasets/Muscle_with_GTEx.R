#'---
#' title: Create Muscle dataset combined with GTEx
#' author: Christian Mertes
#' wb:
#'  input:
#'    - sampleAnno: '`sm config["SAMPLE_ANNOTATION"]`'
#'    - gtexAnno:   '/s/project/crg_seq_data/resource/sample_anno_gtex.RDS'
#'  output:
#'   - colData: '`sm config["PROC_DATA"] + "/annotations/MUSCLE_GTEx.tsv"`'
#'   - wBhtml:  '`sm config["htmlOutputPath"] + "/annotations/MUSCLE_GTEx.html"`'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

if(FALSE){
  snakemake <- readRDS("tmp/snakemake.RDS")
  source(".wBuild/wBuildParser.R")
  parseWBHeader("./Scripts/DefineDatasets/gtex_tissue.R", tissue="Lung")
  annoFile <- "/s/project/crg_seq_data/raw_data/URDCAT_sample_annotation.csv"
  outFile  <-  "Data/annotations/BLOOD_GTEx.tsv"
  gtexFile <- "/s/project/crg_seq_data/resource/sample_anno_gtex.RDS"
}

saveRDS(snakemake, "tmp/snakemake.RDS")

#+ load main config, echo=FALSE
source("./src/r/config.R", echo=FALSE)

#+ input
outFile       <- snakemake@output$colData
annoFile      <- snakemake@input$sampleAnno
gtexFile      <- snakemake@input$gtexAnno

#+ dataset name
name <- gsub(".tsv$", "", basename(outFile))

#'
#' # Load and merge Annotations
#'
name
anno   <- fread(annoFile)
gtexDT <- readRDS(gtexFile)

#' Clean input variables

#' GTEx
gtexDT[,body_site:=gsub("_$", "", gsub("[\\s-()]+", "_", body_site, perl=TRUE))]
gtexDT <- gtexDT[analyte_type == "RNA:Total RNA"]
gtexDT2merge <- gtexDT[body_site=="Muscle_Skeletal", .(SAMPID=submitted_sample_id, run)]

gtexDT2merge <- gtexDT2merge[!grepl('_rep[0-9]+$', SAMPID)] # remove replicates
gtexDT2merge[,sampleID:=gsub("^([^-]+-[^-]+)-.*", "\\1", SAMPID, perl=TRUE)]
gtexFinal <- gtexDT2merge[,.(sampleID=sampleID, bamFile=paste0("/s/project/sra-download/bamfiles/", run, ".bam"))]
gtexFinal <- gtexFinal[!duplicated(sampleID)][file.exists(bamFile)]

#' CRG dataset
annoFinal <- anno[TISSUE == "Muscle - Skeletal", .(sampleID=RNA_fi, bamFile=paste0("/s/project/crg_seq_data/raw_data/RNA_seq_bams/", RNA_file))]

#' Merge final dataset
#' 
#' TODO: remove subsetting later
#' 
set.seed(42)
colData <- rbind(annoFinal, gtexFinal[sample(c(T, F), .N, replace=TRUE, prob=c(0.3, 0.7))])

#'
#' ## Dataset: `r name`
#'
name

#+ setting-gtex-ids, echo=FALSE
colData[,condition:=sampleID]

#+ echo=FALSE
finalTable <- colData


#'
#' ## Final sample table `r name`
#'

#+ savetable
setcolorder(finalTable, unique(c(
  "sampleID", "condition", "bamFile", colnames(finalTable))))

DT::datatable(finalTable, options=list(scrollX=TRUE))

dim(finalTable)
write_tsv(finalTable, file=outFile)
