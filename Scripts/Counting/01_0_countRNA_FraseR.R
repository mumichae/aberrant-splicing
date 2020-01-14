#'---
#' title: Count RNA data with FraseR (Part 0)
#' author: Luise Schuller
#' wb:
#'  params:
#'   - workers: 5
#'   - threads: 5
#'   - internalThreads: 3
#'   - progress: FALSE
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'   - colData: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/annotations/{dataset}.tsv"`'
#'  output:
#'   - splitCounts_tsv: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/splitCounts.tsv.gz"`'
#'   - gRanges_only: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_splitCounts_only.rds"`'
#'   - spliceSites: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/cache/raw-{dataset}/spliceSites_splitCounts.rds"`'
#'   - fdsobj:  '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/fds-object.RDS"`'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FraseR_01_0.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FraseR_01_0.snakemake")

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

# Create FraseR dataset
register(MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress))
colData <- fread(colDataFile)
fds <- FraseRDataSet(colData,
                     workingDir = workingDir,
                     name       = paste0("raw-", dataset))

countDir <- file.path(workingDir(fds), "savedObjects", 
                      nameNoSpace(name(fds)))

# Count split reads
splitCounts <- getSplitReadCountsForAllSamples(fds=fds,
                                               NcpuPerSample=iThreads,
                                               junctionMap=NULL,
                                               recount=params$recount,
                                               BPPARAM=bpparam(),
                                               genome=NULL,
                                               outFile=file.path(countDir,
                                                                 "splitCounts.tsv.gz"))





splitCounts_gRanges <- granges(splitCounts) %>% annotateSpliceSite

saveRDS(splitCounts_gRanges, snakemake@output$gRanges_only)


### Extracting splitSiteCoodinates

# extract donor and acceptor sites
spliceSiteCoords <- extractSpliceSiteCoordinates(splitCounts_gRanges, fds)

saveRDS(spliceSiteCoords, snakemake@output$spliceSites)


fds <- saveFraseRDataSet(fds)

message(date(), ": In total ", length(spliceSiteCoords),
        " splice sites (acceptor/donor) will be counted ...")

#   - gRanges: '`sm parser.getProcDataDir() + 
#                   "/aberrant_splicing/datasets/cache/raw-{dataset}/gRanges_splitCounts.rds"`'
#saveRDS(splitCounts, snakemake@output$gRanges)
