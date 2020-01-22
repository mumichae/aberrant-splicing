#'---
#' title: Filter and clean dataset
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 1
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - psiSS:  '`sm parser.getProcDataDir()+ 
#'                  "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/psiSite.h5"`'
#'  output:
#'   - fds: '`sm parser.getProcDataDir() +
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/fds-object.RDS"`'
#'   - done: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/filter.done" `'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_03.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_03.snakemake")

source("Scripts/_helpers/config.R")
opts_chunk$set(fig.width=12, fig.height=8)

#+ input
dataset    <- snakemake@wildcards$dataset
colDataFile <- snakemake@input$colData
workingDir <- snakemake@params$workingDir
bpWorkers   <- min(max(extract_params(bpworkers()), 1),
                   as.integer(extract_params(snakemake@params$workers)))
params <- snakemake@config$aberrantSplicing


fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))

# Apply filter
minExpressionInOneSample <- params$minExpressionInOneSample
fds <- filterExpression(fds, 
                        minExpressionInOneSample = minExpressionInOneSample,
                        filter=FALSE)
devNull <- saveFraseRDataSet(fds)

# Keep junctions that pass filter
name(fds) <- dataset
if (params$filter == TRUE) {
    filtered <- mcols(fds, type="j")[,"passed"]
    fds <- fds[filtered,]
    message(paste("filtered to", nrow(fds), "junctions"))
}
fds <- saveFraseRDataSet(fds)
file.create(snakemake@output$done)
