#'---
#' title: Count Split Reads
#' author: Luise Schuller
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets"`'
#'  input:
#'   - done_fds: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/cache/raw-{dataset}/fds.done" `'
#'  output:
#'   - done_sample_splitCounts: '`sm parser.getProcDataDir() + 
#'                "/aberrant_splicing/datasets/cache/raw-{dataset}"
#'                +"/sample_tmp/splitCounts/sample_{sample_id}.done"`'
#'  threads: 1
#'  type: script
#'---
saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_01_1.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_01_1.snakemake")

source("Scripts/_helpers/config.R")

dataset    <- snakemake@wildcards$dataset
workingDir <- snakemake@params$workingDir
params <- snakemake@config$aberrantSplicing


# Read FRASER object
fds <- loadFraseRDataSet(dir=workingDir, name=paste0("raw-", dataset))

# Get sample id from wildcard
sample_id <- snakemake@wildcards[["sample_id"]]

# Count splitReads for given sample id
sample_result <- countSplitReads(sampleID=sample_id, 
                                 fds=fds,
                                 NcpuPerSample = snakemake@threads,
                                 recount=params$recount)

message(date(), ": ", dataset, ", ", sample_id,
        " no. splice junctions (split counts) = ", length(sample_result))

file.create(snakemake@output$done_sample_splitCounts)
