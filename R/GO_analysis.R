#' Get gene symbols from ensemble gene names
#' @param geneNames a character vector
#' @param dataset default human: hsapiens_gene_ensembl,
#' for zebrafish drerio_gene_ensembl, for yeast scerevisiae_gene_ensembl
#' @param biomart default "ensembl"
#' @importFrom biomaRt getBM
#' @importFrom biomaRt useEnsembl
#' @return a data.table of geneNames and symbols (2 columns)
getAllORFGeneSymbols <- function(geneNames, dataset, biomart = "ensembl"){
  ensembl <- useEnsembl(biomart = "ensembl", dataset = dataset)
  uniqueGenes <- unique(geneNames)
  geneHits <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                    filters = 'ensembl_gene_id', values = uniqueGenes, mart = ensembl)
  group2 <- data.table::chmatch(geneNames, geneHits$ensembl_gene_id)
  return(data.table(geneNames = geneNames, symbol = geneHits$hgnc_symbol[group2]))
}

#' Get Go terms
#' @param geneNames ensembl gene names
#' @param organism scientific name
#' @importFrom biomartr getGO
#' @return a data.table of geneNames and go terms(2 columns)
getORFsGoTerms <- function(geneNames, organism){
  old <- geneNames
  geneNames <- unique(geneNames)
  Go <- biomartr::getGO(organism = organism,
                        genes    = geneNames,
                        filters  = "ensembl_gene_id")
  desc <- Go$goslim_goa_description
  return(desc[data.table::chmatch(as.character(old), as.character(geneNames))])
}

#' Guess biomart from organism name
#' @param organism scientific name
#' @param biomart character, default "ensembl"
#' @importFrom biomaRt getBM
#' @importFrom biomaRt useEnsembl
#' @return a character with dataset used
getBiomartFromOrganism <- function(organism, biomart="ensembl") {
  ensembl = useEnsembl(biomart = biomart)
  a <- listDatasets(ensembl)
  guess <- a[grep(pattern = p(unlist(strsplit(organism, " ")), collapse = "|"),
       x = a$dataset, value = FALSE, ignore.case = TRUE),][, 1:2]
  if (nrow(guess) == 0) {
    message(p("Did not find biomart candidate for organism", organism))
    message("Set the 'dataset' argument in orfikDirs(dataset = ) from this list:")
    print(a)
    stop()
  }
  if (nrow(guess) > 1) {
    message(p("Found multiple biomart candidates for organism", organism))
    message("Possibilities were")
    print(guess)
    message("Set the 'dataset' argument in orfikDirs(dataset = ) to correct choice this list:")
    print(a)
    stop()
  }

  message(p("Using biomart dataset: ", guess$dataset))
  return(guess$dataset)
}
