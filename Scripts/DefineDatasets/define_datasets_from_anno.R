#'---
#' title: Create datasets from annotation file
#' author: Christian Mertes
#' wb:
#'  input:
#'    - sampleAnnoFile: '`sm config["SAMPLE_ANNOTATION"]`'
#'    - filemappingFile: '`sm config["SAMPLE_FILE_MAPPING"]`'
#'  output:
#'    - colData: '`sm parser.getProcDataDir() + "/aberrant_splicing/annotations/{dataset}.tsv"`'
#'    - wBhtml:  '`sm parser.getProcDataDir() + "/aberrant_splicing/annotations/{dataset}.html"`'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

saveRDS(snakemake, paste0(snakemake@config$tmpdir, "/AberrantSplicing/FraseR_00.snakemake") )
# snakemake <- readRDS(paste0(snakemake@config$tmpdir, "/AberrantSplicing/FraseR_00.snakemake"))

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
name <- snakemake@wildcards$dataset

#'
#' # Load and merge Annotations
#'
name
anno    <- fread(annoFile)
mapping <- fread(mappingFile)

#' 
#' Prepare input data
#' 
subset_ids <- snakemake@config$fraser_filtered[[name]]
annoSub <- anno[RNA_ASSAY %in% subset_ids]
#annoSub <- anno[anno[, name %in% unlist(strsplit(OUTRIDER_GROUP, split = ',')), by = 1:nrow(anno)]$V1,]
#annoSub <- anno[grepl(paste0("^(.*,)?", name, "?(,.*)$"), snakemake@config$outrider_group)]

#' 
#' Create FraseR annotation for given dataset
#' 
colData <- merge(
    annoSub[,.(sampleID = RNA_ASSAY)],
    mapping[ASSAY == "RNA_ASSAY", .(sampleID=ID, bamFile=FILE)])

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
