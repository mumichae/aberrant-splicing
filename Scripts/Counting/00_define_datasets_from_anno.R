#'---
#' title: Create datasets from annotation file
#' author: Christian Mertes
#' wb:
#'  params:
#'    - ids: '`sm parser.fraser_ids`'
#'    - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'    - fileMappingFile: '`sm parser.getProcDataDir() + "/file_mapping.csv"`'
#'  input:
#'    - sampleAnnoFile: '`sm config["sampleAnnotation"]`'
#'  output:
#'    - colData: '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/annotations/{dataset}.tsv"`'
#'    - wBhtml:  '`sm config["htmlOutputPath"] + 
#'                    "/aberrant_splicing/annotations/{dataset}.html"`'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: hide
#'   code_download: TRUE
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FraseR_00_0.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FraseR_00.snakemake")

#+ load main config, echo=FALSE
source("Scripts/_helpers/config.R", echo=FALSE)

#+ input
outFile       <- snakemake@output$colData
annoFile      <- snakemake@input$sampleAnnoFile
fileMapFile   <- snakemake@params$fileMapping

#+ dataset name
name <- snakemake@wildcards$dataset
anno    <- fread(annoFile)
mapping <- fread(fileMapFile)

subset_ids <- snakemake@params$ids[[name]]
annoSub <- anno[RNA_ID %in% subset_ids]
colData <- merge(
    annoSub[,.(sampleID = RNA_ID)],
    mapping[FILE_TYPE == "RNA_BAM_FILE", .(sampleID=ID, bamFile=FILE_PATH)])

#'
#' ## Dataset: `r name`
#'
#+ echo=FALSE
finalTable <- colData

#'
#' ## Final sample table `r name`
#'
#+ savetable
DT::datatable(finalTable, options=list(scrollX=TRUE))

dim(finalTable)
write_tsv(finalTable, file=outFile)
