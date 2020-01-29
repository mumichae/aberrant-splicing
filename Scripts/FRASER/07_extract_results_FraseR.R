#'---
#' title: Results of FraseR analysis
#' author: Christian Mertes
#' wb:
#'  params:
#'   - workers: 3
#'   - threads: 3
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  input:
#'   - fdsin: '`sm parser.getProcDataDir() +
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}/" +
#'                 "pajdBetaBinomial_psiSite.h5"`'
#'  output:
#'   - resultTable: '`sm parser.getProcDataDir() + 
#'                       "/aberrant_splicing/results/{dataset}_results.tsv"`'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_07.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_07.snakemake")

source("Scripts/_helpers/config.R")
opts_chunk$set(fig.width=12, fig.height=8)

dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir
bpWorkers   <- min(max(extract_params(bpworkers()), 1),
                   as.integer(extract_params(snakemake@params$workers)))
bpThreads   <- as.integer(extract_params(snakemake@params$threads))
register(MulticoreParam(bpWorkers, bpThreads))

params <- snakemake@config$aberrantSplicing


# Load data
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)

# Annotate ranges
fds <- annotateRanges(fds)

# Extract results
resgr <- results(fds,
                 padjCutoff=params$padjCutoff, 
                 zScoreCutoff=params$zScoreCutoff,
                 deltaPsiCutoff=params$deltaPsiCutoff)
res   <- as.data.table(resgr)
saveFraseRDataSet(fds)

# Add features
if(nrow(res) > 0){
  # number of samples per gene and variant  
  res[padjust <= params$padjCutoff,
    numSamplesPerGene := length(unique(sampleID)), by=hgncSymbol]
  res[padjust <= params$padjCutoff, 
    numEventsPerGene :=.N, by="hgncSymbol,sampleID"]
  res[padjust <= params$padjCutoff, 
    numSamplesPerJunc:=length(unique(sampleID)), by="seqnames,start,end"]
  
  # add colData at the end
  res <- merge(res, as.data.table(colData(fds)), by="sampleID")
} else{
  warning(paste0("The aberrant splicing pipeline gave 0 results for the ", dataset, " dataset."))
}

# Results
write_tsv(res, file=snakemake@output$resultTable)

