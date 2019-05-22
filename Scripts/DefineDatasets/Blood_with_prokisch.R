#'---
#' title: Create Blood dataset combined with Prokisch
#' author: Christian Mertes
#' wb:
#'  input:
#'    - sampleAnno: '`sm config["SAMPLE_ANNOTATION"]`'
#'  output:
#'   - colData: '`sm config["PROC_DATA"] + "/annotations/BLOOD_Prokisch.tsv"`'
#'   - wBhtml:  '`sm config["htmlOutputPath"] + "/annotations/BLOOD_Prokisch.html"`'
#'  type: noindex
#' output:
#'  html_document:
#'   code_folding: show
#'   code_download: TRUE
#'---

#+ load main config, echo=FALSE
source("./src/r/config.R", echo=FALSE)

#+ input
outFile       <- snakemake@output$colData
annoFile      <- snakemake@input$sampleAnno

#+ dataset name
name <- gsub(".tsv$", "", basename(outFile))

#'
#' # Load and merge Annotations
#'
name
anno   <- fread(annoFile)

#' ## CRG dataset
annoFinal <- anno[TISSUE == "Whole Blood", .(sampleID=RNA_fi, condition=RNA_fi, bamFile=paste0("/s/project/crg_seq_data/raw_data/RNA_seq_bams/", RNA_file))]

#' ## Define the prokisch blood samples
prokisch_blood_samples <- c(
	"110460R", "110313R",   "70476", "112230R", "105376R", "105377R",
	"105378R", "110311R", "111573R", "112227R", "112233R",
	paste0("MUC", c(c(1449:1544)[-c(81,83)],  c(1650:1743)[-c(28, 64, 72, 81)])))

#' get colData for Prokisch dataset
prokisch_dt <- data.table(
	sampleID=prokisch_blood_samples,
	condition=prokisch_blood_samples,
	bamFile=paste0("/s/project/mitoMultiOmics/raw_data/helmholtz/", prokisch_blood_samples, "/RNAout/paired-endout/stdFilenames/", prokisch_blood_samples, ".bam"))

#' # Merge final dataset
#'
colData <- rbind(annoFinal, prokisch_dt)

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

