#'---
#' title: Count RNA data with FRASER (Part 2)
#' author: Luise Schuller
#' wb:
#'  py:
#'  - |
#'   def getSplitCountFiles(dataset):
#'       ids = parser.fraser_ids[dataset]
#'       file_stump = parser.getProcDataDir() + f"/aberrant_splicing/datasets/cache/raw-{dataset}/sample_tmp/splitCounts/"
#'       return expand(file_stump + "sample_{sample_id}.done", sample_id=ids) 
#'  params:
#'   - workers: 20
#'   - threads: 60
#'   - internalThreads: 3
#'   - progress: FALSE
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'   - sample_counts:  '`sm lambda wildcards: getSplitCountFiles(wildcards.dataset)`'
#'  output:
 #'   - countsJ: '`sm parser.getProcDataDir() +
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - assayJ: '`sm parser.getProcDataDir() +
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/splitCounts/assays.h5"`'
#'   - seJ: '`sm parser.getProcDataDir() +
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/splitCounts/se.rds"`'
#'   - gRanges_only: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_splitCounts_only.rds"`'
#'   - spliceSites: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/spliceSites_splitCounts.rds"`'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_01_2.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_01_2.snakemake")

source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir <- snakemake@params$workingDir
recount <- snakemake@params$recount
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

# Directory where splitCounts.tsv.gz will be saved 
countDir <- file.path(workingDir(fds), "savedObjects", 
                      paste0("raw-", dataset))


# Get and merge splitReads for all sample ids
splitCounts <- getSplitReadCountsForAllSamples(fds=fds,
                                               NcpuPerSample=iThreads,
                                               junctionMap=NULL,
                                               recount=params$recount,
                                               BPPARAM=bpparam(),
                                               genome=NULL,
                                               outFile=file.path(countDir,
                                                                 "splitCounts.tsv.gz"))

message(date(), ": Split counts: length = ", length(splitCounts))

# Annoate of granges from the split counts
splitCounts_gRanges <- FRASER:::annotateSpliceSite(rowRanges(splitCounts))
saveRDS(splitCounts_gRanges, snakemake@output$gRanges_only)

message(date(), ": splitCounts_gRanges: length = ", length(splitCounts_gRanges))


# Extracte splitSiteCoodinates: Extract donor and acceptor sites
spliceSiteCoords <- FRASER:::extractSpliceSiteCoordinates(splitCounts_gRanges, fds)
saveRDS(spliceSiteCoords, snakemake@output$spliceSites)

message(date(), ": spliceSiteCoords: length = ", length(spliceSiteCoords))

message(date(), ": In total ", length(spliceSiteCoords),
        " splice sites (acceptor/donor) will be counted ...")