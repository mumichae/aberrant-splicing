#'---
#' title: Create datasets from annotation file
#' author: Christian Mertes
#' wb:
#'  input:
#'    - sampleAnnoFile: '`sm config["SAMPLE_ANNOTATION"]`'
#'    - filemappingFile: '`sm config["SAMPLE_FILE_MAPPTING"]`'
#'  output:
#'   - colData: '`sm parser.getProcDataDir() + "/annotations/{dataset}.tsv"`'
#'   - wBhtml:  '`sm parser.getProcDataDir() + "/annotations/{dataset}.html"`'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

if(FALSE){
  snakemake <- readRDS("tmp/snakemake.RDS")
  source("FraseR-analysis/.wBuild/wBuildParser.R")
  parseWBHeader("./Scripts/DefineDatasets/define_datasets_from_anno.R", dataset="small")
  annoFile <- "/s/project/crg_seq_data/raw_data/URDCAT_sample_annotation.csv"
  outFile  <-  "Data/annotations/BLOOD_GTEx.tsv"
  gtexFile <- "/s/project/crg_seq_data/resource/sample_anno_gtex.RDS"
}

#+ load main config, echo=FALSE
source("./src/r/config.R", echo=FALSE)

#+ input
outFile       <- snakemake@output$colData
annoFile      <- snakemake@input$sampleAnnoFile
mappingFile   <- snakemake@input$filemappingFile

#+ dataset name
name <- gsub(".tsv$", "", basename(outFile))

#'
#' # Load and merge Annotations
#'
name
anno    <- fread(annoFile)
mapping <- fread(mappingFile)

#' 
#' Prepare input data
#' 
annoSub <- anno[grepl(paste0("^(.*,)?", name, "(,.*)?$"), ANALYSIS_GROUP)]

#' 
#' Create FraseR annotation for given dataset
#' 
colData <- merge(
    annoSub[,.(sampleID=RNA_ID)], 
    mapping[TYPE == "RNA", .(sampleID=ID, bamFile=FILE)])

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
