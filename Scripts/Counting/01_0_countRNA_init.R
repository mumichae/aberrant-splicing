#'---
#' title: Count RNA data with FRASER (Part 0)
#' author: Luise Schuller
#' wb:
#'  params:
#'   - workers: 20
#'   - threads: 60
#'   - internalThreads: 3
#'   - progress: FALSE
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'    - colData: '`sm parser.getProcDataDir() + 
#'                    "/aberrant_splicing/annotations/{dataset}.tsv"`'
#'  output:
#'   - fdsobj:  '`sm parser.getProcDataDir() + 
#'                   "/aberrant_splicing/datasets/savedObjects/raw-{dataset}/fds-object.RDS"`'
#'   - done_fds: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/cache/raw-{dataset}/fds.done" `'
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_01_0.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_01_0.snakemake")

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


# Create initial FRASER object
register(MulticoreParam(bpWorkers, bpThreads, progressbar=bpProgress))
colData <- fread(colDataFile)
fds <- FraseRDataSet(colData,
                     workingDir = workingDir,
                     name       = paste0("raw-", dataset))

# Save initial FRASER dataset
fds <- saveFraseRDataSet(fds)

message(date(), ": FRASER object initialized for ", dataset)

file.create(snakemake@output$done_fds)
