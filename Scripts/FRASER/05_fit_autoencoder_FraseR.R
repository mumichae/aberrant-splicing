#'---
#' title: Fitting the autoencoder
#' author: Christian Mertes
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'   
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  threads: 20
#'  input:
#'   - hyper: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/savedObjects/{dataset}/hyper.done" `'
#'  output:
#'   - fdsout: '`sm parser.getProcDataDir() + 
#'                  "/aberrant_splicing/datasets/savedObjects/{dataset}/predictedMeans_psiSite.h5"`'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_05.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_05.snakemake")

source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

fds <- loadFraseRDataSet(dir=workingDir, name=dataset)

# Fit autoencoder
# run it for every type
correction <- snakemake@config$aberrantSplicing$correction

for(type in psiTypes){
    currentType(fds) <- type
    q <- bestQ(fds, type)
    fds <- fit(fds, q=q, type=type, verbose=TRUE, iterations=15, 
               correction=correction)
    fds <- saveFraseRDataSet(fds)
}

