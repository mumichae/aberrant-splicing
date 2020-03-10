#'---
#' title: Results of FraseR analysis
#' author: Christian Mertes
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'   - workingDir: '`sm parser.getProcDataDir() + "/aberrant_splicing/datasets/"`'
#'  threads: 10
#'  input:
#'   - fdsin: '`sm parser.getProcDataDir() +
#'                 "/aberrant_splicing/datasets/savedObjects/{dataset}/" +
#'                 "pajdBetaBinomial_psiSite.h5"`'
#'  output:
#'   - resultTableJunc: '`sm parser.getProcDataDir() + 
#'                          "/aberrant_splicing/results/{dataset}_results_per_junction.tsv"`'
#'   - resultTableGene: '`sm parser.getProcDataDir() + 
#'                          "/aberrant_splicing/results/{dataset}_results.tsv"`'
#'  type: script
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FRASER_07.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FRASER_07.snakemake")

source("Scripts/_helpers/config.R")
opts_chunk$set(fig.width=12, fig.height=8)

dataset    <- snakemake@wildcards$dataset
fdsFile    <- snakemake@input$fdsin
workingDir <- snakemake@params$workingDir

register(MulticoreParam(snakemake@threads))
# Limit number of threads for DelayedArray operations
setAutoBPPARAM(MulticoreParam(snakemake@threads))

params <- snakemake@config$aberrantSplicing

# Load data and annotate ranges with gene names
fds <- loadFraseRDataSet(dir=workingDir, name=dataset)
fds <- annotateRanges(fds)

# Extract results per junction
res_junc <- results(fds,
                 padjCutoff=params$padjCutoff, 
                 zScoreCutoff=params$zScoreCutoff,
                 deltaPsiCutoff=params$deltaPsiCutoff)
res_junc_dt   <- as.data.table(res_junc)
print('Results per junction extracted')
saveFraseRDataSet(fds)

# Add features 
if(nrow(res_junc_dt) > 0){
  
  # number of samples per gene and variant  
  res_junc_dt[, numSamplesPerGene := uniqueN(sampleID), by = hgncSymbol]
  res_junc_dt[, numEventsPerGene := .N, by = "hgncSymbol,sampleID"]
  res_junc_dt[, numSamplesPerJunc := uniqueN(sampleID), by = "seqnames,start,end"]
  
  # add colData to the results
  res_junc_dt <- merge(res_junc_dt, as.data.table(colData(fds)), by = "sampleID")
  res_junc_dt[, bamFile := NULL]
} else{
  warning("The aberrant splicing pipeline gave 0 results for the ", dataset, " dataset.")
}

# Aggregate results by gene
if(length(res_junc) > 0){
  res_genes_dt <- resultsByGenes(res_junc) %>% as.data.table
  res_genes_dt <- merge(res_genes_dt, as.data.table(colData(fds)), by = "sampleID")
  res_genes_dt[, bamFile := NULL]
  
  # add HPO overlap information
  sa <- fread(snakemake@config$sampleAnnotation)
  if(!is.null(sa$HPO_TERMS)){
    if(!all(is.na(sa$HPO_TERMS))){
      hpo_dt <- fread('https://i12g-gagneurweb.in.tum.de/public/paper/drop_analysis/resource/hpo_genes.tsv.gz')
      f2 <- merge(res_genes_dt[, .(sampleID, hgncSymbol)], 
                hpo_dt[,.(hgncSymbol, HPO_id, HPO_label)], by = 'hgncSymbol')
      if(nrow(f2) > 0){
        f3 <- merge(f2, sa[,.(RNA_ID, HPO_TERMS)], by.x = 'sampleID', by.y = 'RNA_ID')
        f3[, HPO_match := HPO_id %in% unlist(strsplit(HPO_TERMS, split = ',')), by = 1:nrow(f3)]
        f3 <- f3[HPO_match == TRUE]
        if(nrow(f3) > 0){
          f4 <- f3[, .(HPO_id_overlap = paste(HPO_id, collapse = ', '), 
                       HPO_label_overlap = paste(HPO_label, collapse = ', ')), 
                   by = .(sampleID, hgncSymbol)]
          res_genes_dt <- merge(res_genes_dt, f4, by = c('sampleID', 'hgncSymbol'), all.x = TRUE)
        }
      }
    }
  }
} else res_genes_dt <- data.table()

# Results
write_tsv(res_junc_dt, file=snakemake@output$resultTableJunc)
write_tsv(res_genes_dt, file=snakemake@output$resultTableGene)

