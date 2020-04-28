


#' Ribo-seq table grouped by tissue
#' 1st table is filtered on fpkm > 1 per tissue
#' 2nd table is row-mean fpkm per tissue
#' @param riboDbName "Ribofpkm"
#' @param dbOutputNames the 2 output names c("RiboByTissueTF", "RiboByTissueMean")
riboAtlasFPKMTissue <- function(riboDbName = "Ribofpkm",
                                dbOutputNames = c("RiboByTissueTF", "RiboByTissueMean"),
                                onlyMatching = FALSE){
  if(length(dbOutputNames) != 2) stop("dbOutputNames must have 2 character elements")

  # now do per tissue true/false
  rpfSamples <- getRiboRNAInfoTable()
  #1. we have tissues in link

  uniqueTissues <- as.character(unique(rpfSamples$tissue))
  #2. we have all the tables in riboTables
  #3. So for each tissue, find the group, then for each group ->
  riboTable <- readTable(riboDbName)
  idColumns <- getIDColumns(riboTable)
  riboTable <- removeIDColumns(riboTable)
  # if not all
  if (onlyMatching) riboTable <- riboTable[,getRiboMatchedToAll(), with = F]
  # number of id columns used
  riboByTissue <- as.data.table(matrix(nrow = nrow(riboTable),
                                       ncol = length(uniqueTissues)))
  colnames(riboByTissue) <- uniqueTissues
  riboByTissueTemp <- riboByTissue
  #4. rowSum(riboColumns > 1) > 1
  #5. So if at least 2 samples in that tissue have
  # .. fpkm of > 1 on that uorf, it will be true
  for(i in uniqueTissues){
    print(i)
    riboColumns <- riboTable[,rpfSamples$tissue == i, with = F]
    riboByTissue[,i] <- rowSums(riboColumns > 1) > 1
  }
  riboByTissue <- data.table(idColumns, riboByTissue)
  insertTable(Matrix = riboByTissue, tableName = dbOutputNames[1],  rmOld = T)

  #now get mean value instead of true/false
  riboByTissueMean <- riboByTissueTemp
  rm(riboByTissue); rm(riboByTissueTemp)
  # rowMeans(riboColumns) mean per orf in tissue
  for(i in uniqueTissues){
    print(i)
    riboColumns <- riboTable[,rpfSamples$tissue == i, with = F]
    riboByTissueMean[,i] <- rowMeans(riboColumns)
  }
  riboByTissueMean <- data.table(idColumns, riboByTissueMean)
  insertTable(Matrix = riboByTissueMean, tableName = dbOutputNames[2], rmOld = T)
}

#' RNA-seq table grouped by tissue
#' 1st table is filtered on fpkm > 1 per tissue
#' 2nd table is row-mean fpkm per tissue
#' @param nIDColumns 1L transcript names
rnaAtlasFPKMTissue <- function(){
  if(!tableNotExists("RNAByTissueMean")){
    print("rnaAtlasFPKMTissue table exists, remove it if you want to rerun with new values")
    return(NULL)
  }
  # now do per tissue true/false
  rnaSamples <- getRiboRNAInfoTable()
  #1. we have tissues in link

  uniqueTissues <- as.character(unique(rnaSamples$tissue))
  #2. we have all the tables in riboTables
  #3. So for each tissue, find the group, then for each group ->
  rnaTable <- readTable("RNAfpkm")
  idColumns <- getIDColumns(rnaTable)
  rnaTable <- removeIDColumns(rnaTable)

  rnaByTissue <- as.data.table(matrix(nrow = nrow(rnaTable),
                                      ncol = length(uniqueTissues)))
  colnames(rnaByTissue) <- uniqueTissues
  rnaByTissueTemp <- rnaByTissue

  #4. rowSum(riboColumns > 1) > 1
  #5. So if at least 2 samples in that tissue have
  # .. fpkm of > 1 on that uorf, it will be true
  for(i in uniqueTissues){
    print(i)
    rnaColumns <- rnaTable[,rnaSamples$tissue == i, with = F]
    rnaByTissue[,i] <- rowSums(rnaColumns > 1) > 1
  }
  rnaByTissue <- data.table(idColumns, rnaByTissue)
  insertTable(Matrix = rnaByTissue, tableName = "RNAByTissueTF", rmOld = T)

  #now get mean value instead of true/false
  rnaByTissueMean <- rnaByTissueTemp
  rm(rnaByTissue)
  rm(rnaByTissueTemp)
  # rowMeans(riboColumns) mean per orf in tissue
  for(i in uniqueTissues){
    print(i)
    rnaColumns <- rnaTable[,rnaSamples$tissue == i, with = F]
    rnaByTissueMean[,i] <- rowMeans(rnaColumns)
  }
  rnaByTissueMean <- data.table(idColumns, rnaByTissueMean)
  insertTable(Matrix = rnaByTissueMean, tableName = "RNAByTissueMean", rmOld = T)
  return(NULL)
}

teAtlasTissue <- function(inputDT, dbOutputNames = c("TEByTissueMean")) {
  info <- getRiboRNAInfoTable()
  inputDTNon <- removeIDColumns(inputDT)

  tissues <- info$tissue
  uniques <- unique(tissues)

  if(!is.null(inputDT$uorfIDs)){
    dt <- data.table(uorfIDs = inputDT$uorfIDs,
                     txNames = inputDT$txNames)
  } else {
    dt <- data.table(txNames = inputDT$txNames)
  }

  for(tissue in uniques) {
    which <- tissues == tissue
    dt <- data.table(dt, rowMeans(inputDTNon[,which, with = F]))
  }

  if(!is.null(inputDT$uorfIDs)){
    colnames(dt) <- c("uorfIDs", "txNames", as.character(uniques))
  } else {
    colnames(dt) <- c("txNames", as.character(uniques))
  }

  insertTable(dt, dbOutputNames, rmOld = T)
  return(NULL)
}
