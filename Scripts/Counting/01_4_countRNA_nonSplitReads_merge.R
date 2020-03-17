#'---
#' title: Merge Nonsplit Counts
#' author: Luise Schuller
#' wb:
#'  py:
#'  - |
#'   def getNonSplitCountFiles(dataset):
#'       ids = parser.fraser_ids[dataset]
#'       file_stump = parser.getProcDataDir() + f"/aberrant_splicing/datasets/cache/raw-{dataset}/sample_tmp/nonSplitCounts/"
#'       return expand(file_stump + "sample_{sample_id}.done", sample_id=ids) 
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  threads: 20
#'  input:
#'   - sample_counts:  '`sm lambda wildcards: getNonSplitCountFiles(wildcards.dataset)`'
#'   - gRangesNonSplitCounts: '`sm parser.getProcDataDir() + 
#'                          "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_NonSplitCounts.rds"`'
#'  output:
#'   - countsSS: '`sm parser.getProcDataDir() +
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_01_4.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_01_4.snakemake")


source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Read FRASER object
fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))

# Read splice site coordinates from RDS
splitCounts_gRanges <- readRDS(snakemake@input$gRangesNonSplitCounts)

# Directory where splitCounts.tsv.gz will be saved 
countDir <- file.path(workingDir(fds), "savedObjects", 
                      paste0("raw-", dataset))

# Get and merge nonSplitReads for all sample ids
nonSplitCounts <- getNonSplitReadCountsForAllSamples(fds=fds, 
                                                     splitCountRanges=splitCounts_gRanges, 
                                                     minAnchor=5, 
                                                     recount=FALSE, 
                                                     longRead=params$longRead,
                                                     outFile=file.path(countDir, 
                                                                       "nonSplitCounts.tsv.gz"))

message(date(), ":", dataset, " nonSplit counts done")
