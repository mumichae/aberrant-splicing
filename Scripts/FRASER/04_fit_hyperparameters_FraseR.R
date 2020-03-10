#'---
#' title: Hyper parameter optimization
#' author: Christian Mertes
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  threads: 20
#'  input:
#'   - filter: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/filter.done" `'
#'  output:
#'   - hyper: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/hyper.done" `'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_04.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_04.snakemake")

source("Scripts/_helpers/config.R")

#+ input
dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

# Load PSI data
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)

# Run hyper parameter optimization
correction <- snakemake@config$aberrantSplicing$correction
for(type in psiTypes){
    message(date(), ": ", type)
    fds <- optimHyperParams(fds, type=type, 
                            correction=correction, 
                            plot = FALSE)
    fds <- saveFraseRDataSet(fds)
}
fds <- saveFraseRDataSet(fds)
file.create(snakemake@output$hyper)

