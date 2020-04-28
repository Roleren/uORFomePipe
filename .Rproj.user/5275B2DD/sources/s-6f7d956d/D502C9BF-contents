#' Group feature table by tissue
#' 1st table is filtered on fpkm > 1 per tissue
#' 2nd table is row-mean fpkm per tissue
#' @param df.rfp an ORFik experiment
#' @param table "Ribofpkm"
#' @param dbOutputNames the 2 output names c("RiboByTissueTF", "RiboByTissueMean")
atlasTissue <- function(df.rfp,
                        table = readTable("fpkmRFP"),
                        dbOutputNames = c("RiboByTissueTF", "RiboByTissueMean")) {
  if(length(dbOutputNames) != 2) stop("dbOutputNames must have 2 character elements")
  if (nrow(df.rfp) != ncol(table)) stop("Not matching rows and columns in tissue creation!")
  # now do per tissue true/false
  samples <- as.data.table(df.rfp)
  samples$tissue <- samples$stage
  #1. we have tissues in link
  samples[, variants := .N, by = tissue]

  uniqueTissues <- as.character(unique(samples$tissue))
  #2. we have all the tables in riboTables
  #3. So for each tissue, find the group, then for each group ->

  # number of id columns used
  riboByTissue <- as.data.table(matrix(nrow = nrow(table),
                                       ncol = length(uniqueTissues)))
  colnames(riboByTissue) <- uniqueTissues
  riboByTissueTemp <- riboByTissue

  j = 1
  for(i in uniqueTissues){
    print(i)
    riboColumns <- table[,samples$tissue == i, with = F]
    if (samples$variants[j] > 1) {
      riboByTissue[,i] <- rowSums(riboColumns > 1) > 1
    } else {
      riboByTissue[,i] <- rowSums(riboColumns > 1) > 0
    }
    j = j + 1
  }
  riboByTissue <- data.table(riboByTissue)
  insertTable(Matrix = riboByTissue, tableName = dbOutputNames[1],  rmOld = T)
  #now get mean value instead of true/false
  riboByTissueMean <- riboByTissueTemp
  # rowMeans(riboColumns) mean per orf in tissue
  for(i in uniqueTissues){
    print(i)
    riboColumns <- table[,samples$tissue == i, with = F]
    riboByTissueMean[,i] <- rowMeans(riboColumns)
  }
  riboByTissueMean <- data.table(riboByTissueMean)
  insertTable(Matrix = riboByTissueMean, tableName = dbOutputNames[2], rmOld = T)
}
