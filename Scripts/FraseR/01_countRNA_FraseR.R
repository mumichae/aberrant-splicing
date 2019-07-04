#'---
#' title: Count RNA data with FraseR
#' author: Christian Mertes
#' wb:
#'  params:
#'   - internalThreads: 3
#'   - progress: FALSE
#'  input:
#'   - colData: '`sm config["PROC_DATA"] + "/annotations/{dataset}.tsv"`'
#'  output:
#'   - fdsobj:  '`sm config["PROC_DATA"] + "/datasets/savedObjects/raw-{dataset}/fds-object.RDS"`'
#'   - countsJ: '`sm config["PROC_DATA"] + "/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsS: '`sm config["PROC_DATA"] + "/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'   - wBhtml:  '`sm config["htmlOutputPath"] + "/FraseR/{dataset}_counting.html"`'
#'  threads: 20
#'  type: noindex
#'---

if(FALSE){
    source(".wBuild/wBuildParser2.R")
    wildcards <- list(dataset="BLOOD_Prokisch")
    parseWBHeader2("./Scripts/FraseR/01_countRNA_FraseR.R", wildcards=wildcards)
    slot(snakemake, "wildcards", check=FALSE) <- wildcards
}

#+ echo=FALSE
source("./src/r/config.R")

#+ input,
dataset     <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir  <- dirname(dirname(dirname(snakemake@output$countsJ)))
bpWorkers   <- min(bpworkers(), as.integer(snakemake@threads))
bpProgress  <- as.logical(snakemake@params[[1]]$progress)
iThreads    <- min(as.integer(bpworkers() / 5), snakemake@params[[1]]$internalThreads)

#'
#' # Dataset
#+ echo=TRUE
dataset

#+ echo=FALSE
colData <- fread(colDataFile)
DT::datatable(colData, options=list(scrollX=TRUE))

#'
#' Counting the dataset
#'
fds <- FraseRDataSet(colData,
        workingDir = workingDir,
        name       = paste0("raw-", dataset),
        parallel   = MulticoreParam(bpWorkers, bpWorkers*3,
                progressbar=bpProgress))
fds <- countRNAData(fds, NcpuPerSample=iThreads, recount=TRUE, minAnchor=5)
fds <- saveFraseRDataSet(fds)

#'
#' FraseR object after counting
#'
fds

