#'---
#' title: Full FraseR analysis over all datasets
#' wb:
#'  params:
#'   - tmpdir: '`sm drop.getMethodPath(METHOD, "tmp_dir")`'
#'  input:
#'   - summaryPlots: '`sm expand(config["htmlOutputPath"] + "/aberrant_splicing/FraseR/{dataset}_summary.html", dataset=config["aberrantSplicing"]["groups"])`'
#'   - resultTable: '`sm expand(parser.getProcDataDir()+ "/aberrant_splicing/results/{dataset}_results.tsv", dataset=config["aberrantSplicing"]["groups"])`'
#' output:
#'  html_document
#'---

saveRDS(snakemake, file.path(snakemake@params$tmpdir, "FraseR_99.snakemake"))
# snakemake <- readRDS(".drop/tmp/AS/FraseR_99.snakemake")

#+ input
allResults <- snakemake@input$resultTable
datasets <- gsub("_results.tsv$", "", basename(allResults))

#+ echo=FALSE, results="asis"
devNull <- sapply(datasets, function(name){
    cat(paste0(
        "<h1>Dataset: ", name, "</h1>",
        "<p>",
        "</br>", "<a href='FraseR/", name, "_countSummary.html'   >Count Summary</a>",
        "</br>", "<a href='FraseR/", name, "_summary.html'        >FRASER Summary</a>",
        "</br>", "</p>"
    ))
})
