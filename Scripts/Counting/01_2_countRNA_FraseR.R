#'---
#' title: Count RNA data with FRASER (Part 2)
#' author: Luise Schuller
#' wb:
#'  py:
#'  - |
#'   def getSplitCountFiles(dataset):
#'       ids = parser.fraser_ids[dataset]
#'       file_stump = parser.getProcDataDir() + f"/aberrant_splicing/datasets/cache/raw-{dataset}/sample_tmp/"
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
#'   - splitCounts_tsv: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/splitCounts.tsv.gz"`'
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
datasetname <- snakemake@wildcards$dataset


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
                                               recount=recount,
                                               BPPARAM=bpparam(),
                                               genome=NULL,
                                               outFile=file.path(countDir,
                                                                 "splitCounts.tsv.gz"))

# Annoate of granges from the split counts
splitCounts_gRanges <- FRASER:::annotateSpliceSite(granges(splitCounts))

saveRDS(splitCounts_gRanges, snakemake@output$gRanges_only)


# Extracte splitSiteCoodinates: Extract donor and acceptor sites
spliceSiteCoords <- FRASER:::extractSpliceSiteCoordinates(splitCounts_gRanges, fds)
saveRDS(spliceSiteCoords, snakemake@output$spliceSites)

message(date(), ": In total ", length(spliceSiteCoords),
        " splice sites (acceptor/donor) will be counted ...")

#  py:
#  - |
#   def getSplitCountFiles(dataset):
#       ids = parser.fraser_ids[dataset]
#       file_stump = parser.getProcDataDir() + f"/aberrant_splicing/datasets/cache/splicedCounts/"
#       return expand(file_stump + "splicedCounts-{sample_id}.RDS", sample_id=ids) 