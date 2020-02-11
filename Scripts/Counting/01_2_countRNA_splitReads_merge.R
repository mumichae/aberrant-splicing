#'---
#' title: Merge Split Counts
#' author: Luise Schuller
#' wb:
#'  py:
#'  - |
#'   def getSplitCountFiles(dataset):
#'       ids = parser.fraser_ids[dataset]
#'       file_stump = parser.getProcDataDir() + f"/aberrant_splicing/datasets/cache/raw-{dataset}/sample_tmp/splitCounts/"
#'       return expand(file_stump + "sample_{sample_id}.done", sample_id=ids) 
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  threads: 20
#'  input:
#'   - sample_counts: '`sm lambda wildcards: getSplitCountFiles(wildcards.dataset)`'
#'  output:
 #'   - countsJ: '`sm parser.getProcDataDir() +
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
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
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))

# Read FRASER object
fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))

# Directory where splitCounts.tsv.gz will be saved 
countDir <- file.path(workingDir, "savedObjects", paste0("raw-", dataset))


# Get and merge splitReads for all sample ids
splitCounts <- getSplitReadCountsForAllSamples(fds=fds,
                                               recount=FALSE,
                                               outFile=file.path(countDir,
                                                                 "splitCounts.tsv.gz"))

# Annotate of granges from the split counts
splitCounts_gRanges <- FRASER:::annotateSpliceSite(rowRanges(splitCounts))
saveRDS(splitCounts_gRanges, snakemake@output$gRanges_only)


# Extract splitSiteCoodinates: Extract donor and acceptor sites
spliceSiteCoords <- FRASER:::extractSpliceSiteCoordinates(splitCounts_gRanges, fds)
saveRDS(spliceSiteCoords, snakemake@output$spliceSites)


message(date(), ": ", dataset, " total no. splice junctions = ", 
        length(splitCounts))
