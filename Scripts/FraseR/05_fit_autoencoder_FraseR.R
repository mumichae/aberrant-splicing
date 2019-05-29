#'---
#' title: Fitting the autoencoder
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 20
#'   - threads: 20
#'   - progress: FALSE
#'  input:
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/FraseR/{dataset}_hyper_parameter_optimization.html"`'
#'  output:
#'   - fdsout: '`sm config["PROC_DATA"] + "/datasets/savedObjects/{dataset}/predictedMeans_psiSite.h5"`'
#'   - wBhtml: '`sm config["htmlOutputPath"] + "/FraseR/{dataset}_autoencoder_fit.html"`'
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
BPPARAM <- MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress)
parallel(fds) <- BPPARAM
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

    # subset fitting
    curX <- x(fds, type=type, all=TRUE, center=FALSE, noiseAlpha=NULL)
    xsd <- colSds(as.matrix(curX))
    nMostVarJuncs <- which(xsd > sort(xsd, TRUE)[min(length(xsd), 30000)])
    exMask <- logical(length(xsd))
    exMask[sample(nMostVarJuncs, min(length(xsd), 15000))] <- TRUE

    featureExclusionMask(fds) <- exMask
    print(table(featureExclusionMask(fds)))

    # run autoencoder
    fds <- fitAutoencoder(fds, q=q, type=type, verbose=TRUE, BPPARAM=BPPARAM, iterations=15, nrDecoderBatches=5)

    # save autoencoder fit
    fds <- saveFraseRDataSet(fds)
}


#'
#' # Save results
#'
fds <- saveFraseRDataSet(fds)

