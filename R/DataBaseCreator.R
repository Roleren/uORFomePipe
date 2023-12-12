
#' Create 1 column of all unique ids from uorfID folder
createUniqueIDs <- function(idFolder) {
  idFiles <- list.files(idFolder, full.names = TRUE)
  if (length(idFiles) == 0) {
    stop("idFiles can not have 0 files!")
  }
  if (tableNotExists("uniqueIDs")) {
    j = 1
    for(i in idFiles){
      uorfID <- unique(readRDS(i))
      if (j == 1) {
        allUniqueIDs <- uorfID
      }else{
        matching <- data.table::`%chin%`(uorfID, allUniqueIDs)
        toAdd <- uorfID[which(matching == FALSE)]
        allUniqueIDs <- c(allUniqueIDs, toAdd)
      }
      j <- j+1
    }
    allUniqueIDs <- allUniqueIDs[data.table::chorder(allUniqueIDs)]

    #save(allUniqueIDs,file = "allUniqueIDs.rdata")
    insertTable(Matrix = allUniqueIDs, tableName = "uniqueIDs")
  } else {
    message("uniqueIDs already exist, skipping remake of them")
  }
}

#' convert to gr from string and filter NB!!! put this  in pipeline!!
createGRObjects <- function(dataFolder, leadersFolder){
  if (!file.exists(p(dataFolder,"/uniqueUorfsAsGRWithTx.rdata"))) {
    message("Creating GRanges objects")
    grl <- toGRFromUniqueID(readTable("uniqueIDs"))
    uniqueIDs <- ORFik:::orfID(grl)
    save(grl, file = p(dataFolder, "/uniqueUorfsAsGR.rdata"))
    insertTable(Matrix = uniqueIDs, tableName =  "uniqueIDs", rmOld = TRUE)
    insertTable(Matrix = grl,tableName = "SplittedByExonsuniqueUORFs", rmOld = TRUE)

    # make all spanning cage leader from cage
    allLeadersSpanningLeader(leadersFolder, dataFolder)
    # find tx matching
    linkORFsToTx(dataFolder)
  } else {
    message("GRObjects already exist, skipping remake of them")
  }
  return(invisible(NULL))
}

#' Create cage ATLAS of uORFs, per CAGE sample, does uORF exist in leaders defined
#' by that CAGE
#'
#' Will create a matrix (1 column per CAGE sample) and save to database
createUORFAtlas <- function(idFolder = idFolder,
                            dataFolder = dataFolder) {
  if (!file.exists(p(dataFolder,"/UORFAtlas.rdata"))) {
    idFiles <- list.files(idFolder, full.names = TRUE)
    uorfIDsAllUnique <- readTable("uniqueIDs")
    colnames(uorfIDsAllUnique) = "uorfID"
    uorfAtlas <- as.data.table(matrix(F, nrow = nrow(uorfIDsAllUnique), ncol = length(idFiles)+1))
    uorfAtlas[,"V1" := uorfIDsAllUnique$uorfID]
    colnames(uorfAtlas) = c("uorfID", as.character(1:(length(idFiles))))
    j = 1
    print(paste("creating uorfAtlas at index of max:", length(idFiles)))
    for(i in idFiles) {
      uorfAtlas[, as.character(j) :=  data.table::`%chin%`(uorfIDsAllUnique$uorfID, readRDS(i))]
      j = j+1
      print(j)
    }

    save(uorfAtlas,file = p(dataFolder,"/UORFAtlas.rdata"))
  } else {
    message("UORFAtlas already exist, skipping remake of them")
  }
  return(invisible(NULL))
}

#' Create tissueTable for cage, 1 row per unique uorf
#' Tissue must have at least 2 CAGE libraries supporting the uORF
#' to declare a hit in that tissue. If only one sample, that sample
#' must include uORF to be a hit.
#' @param cageTable an ORFik experiment of cage, or NULL
#' @param dataFolder a path to folder with data
getTissueTable <- function(cageTable, dataFolder){
  if (tableNotExists("tissueAtlasByCage")) {
    if (is.null(cageTable)) { # No cage, make all TRUE
      message("Running pipeline without CAGE data, set to NULL")
      groups <- readTable("experiment_groups")[[1]]
      uorfIDs <- readTable("uniqueIDs"); colnames(uorfIDs) = "uorfID"
      tissueAtlas <- as.data.table(matrix(TRUE, nrow = nrow(uorfIDs), ncol = length(groups)+1))
      tissueAtlas[, V1 := uorfIDs$uorfID]
      colnames(tissueAtlas) <- c("uorfID", groups)

    } else {
      groups <- bamVarName(cageTable, skip.replicate = TRUE, skip.libtype = TRUE)
      cageTable <- as.data.table(cageTable)
      cageTable$tissue <- groups; cageTable$stage <- NULL
      cageTable[, index := .N, by = tissue]
      uniqueTissues <- as.character(unique(cageTable$tissue))
      # load needed tables, and make tissue atlas of cage
      load(p(dataFolder,"/UORFAtlas.rdata"))
      uorfIDs <- readTable("uniqueIDs"); colnames(uorfIDs) = "uorfID"

      if(any(rowSums(uorfAtlas[,2:ncol(uorfAtlas)]) == 0))
        stop("uorfAtlas and unique uorf IDs does not match!")

      finalMatrix <- as.data.table(matrix(nrow =
        nrow(uorfIDs), ncol = length(uniqueTissues)+1))

      finalMatrix[, V1 := uorfIDs$uorfID]
      colnames(finalMatrix) <- c("uorfID", uniqueTissues)
      onlyOneCAGE <- c()
      for(i in 2:(length(uniqueTissues)+1)){
        cageFilestoCheck <- cageTable[tissue == uniqueTissues[i-1]]$index
        if(length(cageFilestoCheck) == 1) onlyOneCAGE <- c(onlyOneCAGE, i)
      }

      for(i in 2:(length(uniqueTissues)+1)){
        cageFilestoCheck <- cageTable[tissue == uniqueTissues[i-1]]$index
        cageFilestoCheck <- cageFilestoCheck + 1
        makeGroupingForColumn <- rowSums(uorfAtlas[,cageFilestoCheck, with = F])
        if (i %in% onlyOneCAGE) {
          finalMatrix[, uniqueTissues[i-1] := makeGroupingForColumn > 0]
        } else {
          finalMatrix[, uniqueTissues[i-1] := makeGroupingForColumn > 0]
        }
      }
      tissueAtlas <- finalMatrix
    }
    saveRDS(tissueAtlas, file = p(dataFolder,"/tissueAtlas.rds"))
    insertTable(Matrix = tissueAtlas,tableName = "tissueAtlasByCage", rmOld = T)
    return("ok tissueAtlassCage")
  }
  return("tissueAtlasCage already exists, stop if you want new")
}
