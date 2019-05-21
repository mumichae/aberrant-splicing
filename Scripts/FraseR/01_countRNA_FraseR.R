#'---
#' title: Count RNA data with FraseR
#' author: Christian Mertes
#' wb:
#'  py:
#'   - |
#'     def get_input_bam_files(wildcards):
#'       list(pd.read_cvs(config["PROC_DATA"] + "/annotations/" + wildcards.datasets + ".tsv"
#'  params:
#'   - workers: 20
#'   - threads: 60
#'   - internalThreads: 3
#'   - progress: FALSE
#'  input:
#'   - colData: '`sm config["PROC_DATA"] + "/annotations/{dataset}.tsv"`'
#'   - bamFiles: '`sm get_input_bam_files`'
#'  output:
#'   - fdsobj:  '`sm config["PROC_DATA"] + "/datasets/savedObjects/raw-{dataset}/fds-object.RDS"`'
#'   - countsJ: '`sm config["PROC_DATA"] + "/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsS: '`sm config["PROC_DATA"] + "/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'   - wBhtml:  '`sm config["htmlOutputPath"] + "/FraseR/{dataset}_counting.html"`'
#'  type: noindex
#'---

if(FALSE){
    snakemake <- readRDS("./tmp/snakemake.RDS")
    source(".wBuild/wBuildParser.R")
    parseWBHeader("./Scripts/FraseR/01_countRNA_FraseR.R", dataset="Lung")
    bpWorkers <- min(bpworkers(), 30)
    bpThreads <- 60
    bpProgress <- TRUE
    iThreads  <- min(ceiling(bpworkers())/5, 3)
}

#+ echo=FALSE
source("./src/r/config.R")

#+ input,
dataset     <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir  <- dirname(dirname(dirname(snakemake@output$countsJ)))
bpWorkers   <- min(bpworkers(), as.integer(snakemake@params[[1]]$workers))
bpThreads   <- as.integer(snakemake@params[[1]]$threads)
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
        parallel   = MulticoreParam(bpWorkers, bpThreads,
                progressbar=bpProgress))
fds <- countRNAData(fds, NcpuPerSample=iThreads, recount=TRUE, minAnchor=5)
fds <- saveFraseRDataSet(fds)

#'
#' FraseR object after counting
#'
fds

