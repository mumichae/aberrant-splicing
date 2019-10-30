#'---
#' title: Count RNA data with FraseR
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 20
#'   - threads: 60
#'   - internalThreads: 3
#'   - progress: FALSE
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - colData: '`sm parser.getProcDataDir() + "/aberrant_splicing/annotations/{dataset}.tsv"`'
#'  output:
#'   - fdsobj:  '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/fds-object.RDS"`'
#'   - countsJ: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsS: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'   - wBhtml:  '`sm config["htmlOutputPath"] + "/aberrant_splicing/FraseR/{dataset}_counting.html"`'
#'  type: noindex
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FraseR_01.snakemake") )
# snakemake <- readRDS(".drop/tmp/AE/FraseR_01.snakemake")

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
bpWorkers   <- min(bpworkers(), as.integer(snakemake@params$workers))
bpThreads   <- as.integer(snakemake@params$threads)
bpProgress  <- as.logical(snakemake@params$progress)
iThreads    <- min(as.integer(bpworkers() / 5),
                   as.integer(snakemake@params$internalThreads))

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
fds <- countRNAData(fds, NcpuPerSample=iThreads, recount=FALSE, minAnchor=5)
fds <- saveFraseRDataSet(fds)

#'
#' FraseR object after counting
#'
fds

