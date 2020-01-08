#'---
#' title: Count RNA data with FraseR (Part 1)
#' author: Luise Schuller
#' wb:
#'  params:
#'   - workers: 20
#'   - threads: 60
#'   - internalThreads: 3
#'   - progress: FALSE
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'   - spliceSites: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/spliceSites_splitCounts.rds"`'
#'  output:
#'   - nonSplicedCount_sample : '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/nonSplicedCounts/raw-{dataset}/nonSplicedCounts-{sample_id}.RDS"`' 
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FraseR_01_1.snakemake") )
# snakemake <- readRDS(".drop/tmp/AS/FraseR_01_1.snakemake")

source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir <- snakemake@params$workingDir
bpWorkers   <- min(max(extract_params(bpworkers()), 1),
                   as.integer(extract_params(snakemake@params$workers)))
bpThreads   <- as.integer(extract_params(snakemake@params$threads))
bpProgress  <- as.logical(extract_params(snakemake@params$progress))
iThreads    <- min(max(as.integer(bpWorkers / 5), 1),
                   as.integer(extract_params(snakemake@params$internalThreads)))
params <- snakemake@config$aberrantSplicing
datasetname <- snakemake@wildcards$dataset
ids <- snakemake@config$fraser_ids[[datasetname]]


# Load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

# Create FraseR dataset
register(MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress))
fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))


sample_id <- snakemake@wildcards[["sample_id"]]


## RDS erzeugen
spliceSiteCoords <- readRDS(snakemake@input$spliceSites)

sample_result <- countNonSplicedReads(sample_id, splitCounts = NULL, fds = fds,
                                 NcpuPerSample=1, minAnchor=5, recount=FALSE,
                                 spliceSiteCoords=spliceSiteCoords,
                                 longRead=FALSE)
