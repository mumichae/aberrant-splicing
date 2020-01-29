#'---
#' title: Count RNA data with FRASER (Part 4)
#' author: Luise Schuller
#' wb:
#'  py:
#'  - |
#'   def getNonSplitCountFiles(dataset):
#'       ids = parser.fraser_ids[dataset]
#'       file_stump = parser.getProcDataDir() + f"/aberrant_splicing/datasets/cache/raw-{dataset}/sample_tmp/nonSplitCounts/"
#'       return expand(file_stump + "sample_{sample_id}.done", sample_id=ids) 
#'  params:
#'   - workers: 20
#'   - threads: 60
#'   - internalThreads: 3
#'   - progress: FALSE
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'   - sample_counts:  '`sm lambda wildcards: getNonSplitCountFiles(wildcards.dataset)`'
#'   - gRanges_only: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_splitCounts_only.rds"`'
#'  output:
#'   - countsSS: '`sm parser.getProcDataDir() +
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_01_4.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_01_4.snakemake")


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

# Read FRASER object
fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))

# Read splice site coordinates from RDS
splitCounts_gRanges <- readRDS(snakemake@input$gRanges_only)

# Directory where splitCounts.tsv.gz will be saved 
countDir <- file.path(workingDir(fds), "savedObjects", 
                      paste0("raw-", dataset))

# Get and merge nonSplitReads for all sample ids
nonSplitCounts <- getNonSplitReadCountsForAllSamples(fds=fds, 
                                                     splitCountRanges=splitCounts_gRanges, 
                                                     NcpuPerSample=iThreads, 
                                                     minAnchor=5, 
                                                     recount=FALSE, 
                                                     longRead=params$longRead,
                                                     outFile=file.path(countDir, 
                                                                       "nonSplitCounts.tsv.gz"))

message(date(), ":", dataset, " nonSplit counts done")
