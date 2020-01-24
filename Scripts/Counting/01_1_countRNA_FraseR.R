#'---
#' title: Count RNA data with FRASER (Part 1)
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
#'   - done_fds: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/cache/raw-{dataset}/fds.done" `'
#'  output:
#'   - done_sample: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/cache/raw-{dataset}/sample_tmp/sample_{sample_id}.done"`'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_01_1.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_01_1.snakemake")

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

# Setting path for saving the counts produced from FRASER
countDir <- file.path(workingDir(fds), "savedObjects", 
                      paste0("raw-", dataset))

# Get sample id from wildcard
sample_id <- snakemake@wildcards[["sample_id"]]

# Count splitReads for given sample id
sample_result <- countSplitReads(sampleID=sample_id, 
                                 fds=fds,
                                 NcpuPerSample=iThreads,
                                 genome=NULL,
                                 recount=params$recount)

message(date(), sample_id, ": length = ", length(sample_result))

file.create(snakemake@output$done_sample)