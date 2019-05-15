#'---
#'
#' title: Full FraseR analysis over all datasets
#' author: Christian Mertes
#' wb:
#'  input:
#'   - resultTable: '`sm expand(config["PROC_DATA"] + "/processedData/results/{dataset}_results.tsv", dataset=config["DATASETS"])`'
#' output:
#'  html_document
#'---

#+ input
allResults <- snakemake@input$resultTable
datasets <- gsub("_results.tsv$", "", basename(allResults))

#+ echo=FALSE, results="asis"
devNull <- sapply(datasets, function(name){
    cat(paste0(
        "<h1>Dataset: ", name, "</h1>",
        "<p>",
        "</br>", "<a href='FraseR/", name, "_counting.html'                     >01. Counting</a>",
        "</br>", "<a href='FraseR/", name, "_psi_value_calculation.html'        >02. Calculate PSI values</a>",
        "</br>", "<a href='FraseR/", name, "_filterExpression.html'             >03. Filter dataset by expression</a>",
        "</br>", "<a href='FraseR/", name, "_hyper_parameter_optimization.html' >04. Hyper parameter optimization</a>",
        "</br>", "<a href='FraseR/", name, "_autoencoder_fit.html'              >05. Autoencoder fitting</a>",
        "</br>", "<a href='FraseR/", name, "_stat_calculation.html'             >06. Calculate P-values and other stats</a>",
        "</br>", "<a href='FraseR/", name, "_results.html'                      >07. Final results</a>",
        "</br>", "</p>"
    ))
})
