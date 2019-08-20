#'---
#' title: Fitting the autoencoder
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 20
#'   - threads: 20
#'   - progress: FALSE
#'  input:
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/aberrant_splicing/FraseR/{dataset}_hyper_parameter_optimization.html"`'
#'  output:
#'   - fdsout: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/savedObjects/{dataset}/predictedMeans_psiSite.h5"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/aberrant_splicing/FraseR/{dataset}_autoencoder_fit.html"`'
#'  type: noindex
#'---
##
## TODO:
##   Add a link to the fraser object to have a proper chain of events in wbuild
##

if(FALSE){
    snakemake <- readRDS("./tmp/snakemake.RDS")
    source(".wBuild/wBuildParser.R")
    parseWBHeader("./Scripts/FraseR/05_fit_autoencoder_FraseR.R", dataset="SimulationDM")
    parseWBHeader("./Scripts/FraseR/05_fit_autoencoder_FraseR.R", dataset="example")
}

#+ echo=FALSE
source("./src/r/config.R")

#+ input
dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@output$fdsout
workingDir <- dirname(dirname(dirname(fdsFile)))
bpWorkers  <- min(bpworkers(), as.integer(snakemake@params[[1]]$workers))
bpThreads  <- min(bpworkers(), as.integer(snakemake@params[[1]]$threads))
bpProgress <- snakemake@params[[1]]$progress


#'
#' # Load PSI data
#+ echo=TRUE
dataset

#+ echo=FALSE
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)
bpparam <- MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress)
parallel(fds) <- bpparam
dim(fds)

#'
#' # Fit autoencoder
#'

#'
#' run it for every type
#'
for(type in psiTypes){

    # set current type
    currentType(fds) <- type
    curDims <- dim(K(fds, type))
    q <- bestQ(fds, type)
    probE <- max(0.001, min(1,30000/curDims[1]))

    # subset fitting
    featureExclusionMask(fds) <- sample(c(TRUE, FALSE), curDims[1],
            replace=TRUE, prob=c(probE, 1-probE))
    print(table(featureExclusionMask(fds)))

    # run autoencoder
    fds <- fitAutoencoder(fds, q=q, type=type, verbose=TRUE, BPPARAM=bpparam, iterations=15)

    # save autoencoder fit
    fds <- saveFraseRDataSet(fds)
}


#'
#' # Save results
#'
fds <- saveFraseRDataSet(fds)

