#'---
#' title: Count RNA data with FraseR (Part 3)
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
#'   - gRanges_splitCount: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_splitCounts.rds"`'
#'   - nonSplitCounts_tsv: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/nonSplitCounts.tsv.gz"`'
#'  output:
#'   - countsJ: '`sm parser.getProcDataDir() +
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsS: '`sm parser.getProcDataDir() +
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FraseR_01_3.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FraseR_01_3.snakemake")


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

# Load libraries
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))

splitCounts <- readRDS(snakemake@input$gRanges_splitCount)
nonSplitCounts <- fread(snakemake@input$nonSplitCounts_tsv)
nonSplitCounts <- makeGRangesFromDataFrame(nonSplitCounts, keep.extra.columns = TRUE)


fds <- addCountsToFraseRDataSet(fds=fds, splitCounts=splitCounts, nonSplitCounts=nonSplitCounts)

fds <- saveFraseRDataSet(fds)