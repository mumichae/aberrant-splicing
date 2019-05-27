# R function
#
# author daniel, mertes
#------------------------

check_for_na_urls <- function(url, input){
	url[is.na(input) | input == ''] <- ""
	return(url)
}


get_omim_url <- function(omimID){
	check_for_na_urls(input = omimID,
			url = paste0("http://omim.org/entry/", omimID)
	)
}


get_entrez_url <- function(entrezID){
	check_for_na_urls(input = entrezID,
			url = paste0("http://www.ncbi.nlm.nih.gov/gene/", entrezID)
	)
}


get_locus_url <- function(locus){
	check_for_na_urls(input = locus,
			url = paste0("http://genome-euro.ucsc.edu/cgi-bin/hgTracks?position=", locus)
	)
}


get_hgnc_search_url= function(gene){
    check_for_na_urls(input = gene,
        url = paste0("http://www.genenames.org/cgi-bin/gene_search?search=",gene)
    )
}


get_jensen_disease_url= function(gene){
    check_for_na_urls(input = gene,
        url = paste0("http://diseases.jensenlab.org/Search?query=",gene)
    )
}


get_genecards_search_url= function(gene){
	check_for_na_urls(input = gene,
			url = paste0("http://www.genecards.org/cgi-bin/carddisp.pl?gene=",gene)
	)
}


print_url_as_html_link= function(url, linkname){
	check_for_na_urls(input = linkname,
			url = paste0('<a target="_blank" href=\"',url, '\">',linkname,'</a>')
	)
}


get_html_links_for_gene <- function(...){
	get_html_link(...)
}


get_html_link <- function(input, website='genecards', asLink=TRUE){
	url <- switch(website,
			genecards = get_genecards_search_url(input),
			hgnc      = get_hgnc_search_url(input),
			omim      = get_omim_url(input),
			entrez    = get_entrez_url(input),
			locus     = get_locus_url(input)
	)

	url <- check_for_na_urls(input=input, url=url)

	if(asLink){
	    url <- print_url_as_html_link(url, input)
	}

	return(url)
}
