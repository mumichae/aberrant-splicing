#'---
#' title: Count RNA data with FraseR (Part 2)
#' author: Luise Schuller
#' wb:
#'  py:
#'  - |
#'   def getNonSplitCountFiles(dataset):
#'       ids = parser.fraser_ids[dataset]
#'       file_stump = parser.getProcDataDir() + f"/aberrant_splicing/datasets/cache/nonSplicedCounts/raw-{dataset}/"
#'       return expand(file_stump + "nonSplicedCounts-{sample_id}.RDS", sample_id=ids) 
#'  params:
#'   - workers: 20
#'   - threads: 60
#'   - internalThreads: 3
#'   - progress: FALSE
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'   - sample_counts:  '`sm lambda wildcards: getNonSplitCountFiles(wildcards.dataset)`'
#'  output:
#'   - nonSplitCounts_tsv: '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/nonSplitCounts.tsv.gz"`'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FraseR_01_2.snakemake") )
# snakemake <- readRDS(".drop/tmp/AS/FraseR_01_2.snakemake")


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

#count non spliced reads for every samples
countList <- lapply(file.path(snakemake@input$sample_counts),
                    FUN=readRDS)
names(countList) <- samples(fds)
 countList
siteCounts <- mergeCounts(countList, assumeEqual=TRUE)
 siteCounts
mcols(siteCounts)$type <- factor(countList[[1]]$type,
                                 levels = c("Acceptor", "Donor"))

# write tsv
writeCountsToTsv(siteCounts, file=snakemake@output$nonSplitCounts_tsv)
