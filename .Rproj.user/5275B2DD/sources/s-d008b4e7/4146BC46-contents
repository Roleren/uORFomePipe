#' all features of Riboseq and RNAseq
allFeaturesAtlas <- function(){
  #validate ribo and rna
  validateRiboRNAPairs()

  # first sequence features
  getSequenceFeatures()

  # then computeFeatures
  getAllFeaturesFromUorfs()
  # then rfpTables
  riboAtlasFPKMTissue() # fix this is wrong

  # then rna tables
  getRNAFpkms()
  rnaAtlasFPKMTissue()
  # then do te tables
  getTeFeatures()
  teAtlasTissue()

  #cds and 3'
  getCDSFeatures()
  getFeaturesThreeUTRs()
}

#' Main function to fill uORF database
createCatalogueDB <- function() {
  dataBaseFolder <- p(mainFolder,"/dataBase")
  setwd(dataBaseFolder)
  uorfDB <- createDataBase("uorfCatalogue.sqlite")
  createUniqueIDs() # IDs for uORFs
  createGRObjects() # GRanges objects for all uORFs
  createUORFAtlas() # Per CAGE reassigned tx annotation, does uORF exist ?
  getTissueTable()  # Per CAGE tissue, does uORF exist ?
  allFeaturesAtlas() # Ribo-seq and RNA-seq features
}
