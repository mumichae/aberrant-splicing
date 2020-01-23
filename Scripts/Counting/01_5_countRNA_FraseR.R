#'---
#' title: Count RNA data with FRASER (Part 5)
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
#'   - splitCounts_tsv: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/splitCounts.tsv.gz"`'
#'   - nonSplitCounts_tsv: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/nonSplitCounts.tsv.gz"`'
#'  output:
#'   - countsJ: '`sm parser.getProcDataDir() +
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsS: '`sm parser.getProcDataDir() +
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_01_5.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_01_5.snakemake")


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

# Read FRASER object
fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))

# Get splitReads and nonSplitRead counts in order to store them in FRASER object
splitCounts <- fread(snakemake@input$splitCounts_tsv)
splitCounts <- makeGRangesFromDataFrame(splitCounts, keep.extra.columns = TRUE)
nonSplitCounts <- fread(snakemake@input$nonSplitCounts_tsv)
nonSplitCounts <- makeGRangesFromDataFrame(nonSplitCounts, keep.extra.columns = TRUE)

# Add Counts to FRASER dataset
fds <- addCountsToFraseRDataSet(fds=fds, splitCounts=splitCounts, nonSplitCounts=nonSplitCounts)

# Save final FRASER object 
fds <- saveFraseRDataSet(fds)
