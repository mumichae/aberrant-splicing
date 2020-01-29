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
#'   - countsJ:  '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#'   - countsSS: '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#'   - gRanges_only: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_splitCounts_only.rds"`'
#'   - spliceSites: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/spliceSites_splitCounts.rds"`'
#'  output:
#'   - counting_done: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/counting.done" `'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_01_5.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_01_5.snakemake")


source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
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
splitCounts_gRanges <- readRDS(snakemake@input$gRanges_only)
spliceSiteCoords <- readRDS(snakemake@input$spliceSites)

# Get splitReads and nonSplitRead counts in order to store them in FRASER object
splitCounts_h5 <- HDF5Array::HDF5Array(snakemake@input$countsJ, "rawCountsJ")
splitCounts_se <- SummarizedExperiment(
  colData = colData(fds),
  rowRanges = splitCounts_gRanges,
  assays = list(rawCountsJ=splitCounts_h5)
)


nonSplitCounts_h5 <- HDF5Array::HDF5Array(snakemake@input$countsSS, "rawCountsSS")
print(dim(nonSplitCounts_h5))
print(length(spliceSiteCoords))
nonSplitCounts_se <- SummarizedExperiment(
  colData = colData(fds),
  rowRanges = spliceSiteCoords,
  assays = list(rawCountsSS=nonSplitCounts_h5)
)

# Add Counts to FRASER dataset
fds <- addCountsToFraseRDataSet(fds=fds, splitCounts=splitCounts_se,
                                nonSplitCounts=nonSplitCounts_se)

# Save final FRASER object 
fds <- saveFraseRDataSet(fds)

file.create(snakemake@output$counting_done)
