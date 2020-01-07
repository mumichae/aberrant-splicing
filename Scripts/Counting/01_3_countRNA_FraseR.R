#' #'---
#' #' title: Count RNA data with FraseR
#' #' author: Luise Schuller
#' #' wb:
#' #'  params:
#' #'   - workers: 20
#' #'   - threads: 60
#' #'   - internalThreads: 3
#' #'   - progress: FALSE
#' #'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#' #'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#' #'  input:
#' #'   - colData: '`sm parser.getProcDataDir() + 
#' #'                   "/aberrant_splicing/annotations/{dataset}.tsv"`'
#' #'  output:
#' #'   - fdsobj:  '`sm parser.getProcDataDir() + 
#' #'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/fds-object.RDS"`'
#' #'   - countsJ: '`sm parser.getProcDataDir() + 
#' #'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsJ.h5"`'
#' #'   - countsS: '`sm parser.getProcDataDir() + 
#' #'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/rawCountsSS.h5"`'
#' #'  type: script
#' #'---
#' saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FraseR_01.snakemake") )
#' # snakemake <- readRDS(".drop/tmp/AS/FraseR_01.snakemake")
#' 
#' source("./src/r/config.R")
#' 
#' dataset    <- snakemake@wildcards$dataset
#' colDataFile <- snakemake@input$colData
#' workingDir <- snakemake@params$workingDir
#' bpWorkers   <- min(max(extract_params(bpworkers()), 1),
#'                    as.integer(extract_params(snakemake@params$workers)))
#' bpThreads   <- as.integer(extract_params(snakemake@params$threads))
#' bpProgress  <- as.logical(extract_params(snakemake@params$progress))
#' iThreads    <- min(max(as.integer(bpWorkers / 5), 1),
#'                    as.integer(extract_params(snakemake@params$internalThreads)))
#' params <- snakemake@config$aberrantSplicing
#' 


#fds <- saveFraseRDataSet(fds)