#' Get indices of uORFs from genes
#' @param hgncSymbol a character vector of gene symbols
#' @return a integer vector of indices
#' @importFrom biomaRt getBM
#' @importFrom biomaRt useEnsembl
getORFsGeneSymbols <- function(hgncSymbol = "ATF4", refTable = refTable){
  ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
  geneHits <- getBM(attributes = c('ensembl_gene_id', 'ensembl_transcript_id','hgnc_symbol'),
                    filters = 'hgnc_symbol', values = hgncSymbol, mart = ensembl)
  if(nrow(geneHits) == 0) stop(p("could not find any genes with the name: ", hgncSymbol))

  return(which(refTable$geneNames == geneHits[1, 1]))
}

#' Get gene symbols from ensemble gene names
#' @param geneNames a character vector
#' @param dataset default human: hsapiens_gene_ensembl,
#' for zebrafish drerio_gene_ensembl, for yeast scerevisiae_gene_ensembl
#' @importFrom biomaRt getBM
#' @importFrom biomaRt useEnsembl
#' @return a data.table of geneNames and symbols (2 columns)
getAllORFGeneSymbols <- function(geneNames, dataset){

  ensembl <- useEnsembl(biomart = "ensembl", dataset = dataset)
  uniqueGenes <- unique(geneNames)
  geneHits <- getBM(attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                    filters = 'ensembl_gene_id', values = uniqueGenes, mart = ensembl)
  group2 <- data.table::chmatch(geneNames, geneHits$ensembl_gene_id)
  return(data.table(geneNames = geneNames, symbol = geneHits$hgnc_symbol[group2]))
}

#' Get Go terms
#' @param geneNames ensembl gene names
#' @param organism scientifi name
#' @importFrom biomartr getGO
getORFsGoTerms <- function(geneNames, organism = "Homo sapiens"){
  old <- geneNames
  geneNames <- unique(geneNames)
  Go <- biomartr::getGO(organism = organism,
                        genes    = geneNames,
                        filters  = "ensembl_gene_id")
  desc <- Go$goslim_goa_description
  return(desc[data.table::chmatch(as.character(old), as.character(geneNames))])
}
