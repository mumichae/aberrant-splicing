#'---
#' title: Calculate PSI values
#' author: Christian Mertes
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  threads: 20
#'  input:
#'   - counting_done: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/counting.done" `'
#'  output:
#'  - psiSS:     '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/psiSite.h5"`'
#'  type: script
#'--- 

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_02.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_02.snakemake")
source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))

fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))

# Calculating PSI values
fds <- calculatePSIValues(fds)

# FraseR object after PSI value calculation
fds <- saveFraseRDataSet(fds)
