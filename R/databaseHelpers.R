#' Create the database for uORFome
#'
#' Will assing to .GlobalEnv
#' @param name full path name of database, existing or fresh.
#' @return invisible(NULL)
createDataBase <- function(name){
  uorfDB <- dbConnect(RSQLite::SQLite(), name)
  assign(x = "uorfDB", uorfDB, envir = .GlobalEnv)
  return(invisible(NULL))
}

#' Delete the database for uORFome
#'
#' Will assing to .GlobalEnv
#' @inheritParams createDataBase
#' @param uorfDB the defined uORFome sql database, default:
#' get("uorfDB", envir = .GlobalEnv)
#' @return invisible(NULL)
deleteDataBase <- function(name, uorfDB = get("uorfDB", envir = .GlobalEnv)){
  dbDisconnect(uorfDB)
  unlink(name)
  return(invisible(NULL))
}

#' Delete the database for uORFome
#'
#' Will assing to .GlobalEnv
#' @param Matrix the data.table / matrix to insert
#' @inheritParams readTable
#' @param appends rowbind instead of insert
#' @param rmOld  (TRUE) allow replacing existing table (if TRUE)
#' @return invisible(NULL)
insertTable <- function(Matrix, tableName, appends = FALSE, rmOld = FALSE,
                        uorfDB = get("uorfDB", envir = .GlobalEnv)){
  if (rmOld){
    if(!tableNotExists(tableName))
      deleteTable(tableName)
  }
  dbWriteTable(uorfDB, tableName, as.data.table(Matrix),append = appends)
  return(invisible(NULL))
}

#' Read a table from the database
#' @param tableName name of table in sql database
#' @param asGR convert to GRanges
#' @param with.IDs include ID column (remove to make calculations easier)
#' @param uorfDB the defined uORFome sql database, default:
#' get("uorfDB", envir = .GlobalEnv)
#' @return the table as data.table or GRanges
#' @importFrom data.table as.data.table
#' @export
readTable <- function(tableName, asGR = FALSE, with.IDs = TRUE,
                      uorfDB = get("uorfDB", envir = .GlobalEnv)) {
  if (asGR){
    grl <- as.data.table(dbReadTable(uorfDB, tableName))
    return(makeGRangesListFromDataFrame(grl, split.field = "group",
                                        names.field = "group_name",
                                        keep.extra.columns = TRUE))

  } else{
    if (!with.IDs) {
      dt <- as.data.table(dbReadTable(uorfDB,tableName))

      return(removeIDColumns(dt))
    }
    return(as.data.table(dbReadTable(uorfDB,tableName)))
  }
}

#' List current tables in database
#' @inheritParams deleteDataBase
#' @export
listTables <- function(uorfDB = get("uorfDB", envir = .GlobalEnv)){
  sort(dbListTables(uorfDB))
}

#' Does table exist
#' @inheritParams readTable
#' @param name Table name
#' @param exact (TRUE), if FALSE do "grep" search instead
#' @return logical (TRUE if not exists)
tableNotExists <- function(name, exact = TRUE,
                           uorfDB = get("uorfDB", envir = .GlobalEnv)){
  if (exact) return(sum(name %in% listTables()) == 0)

  return(sum(grep(pattern = name, x = listTables())) == 0)
}

#' Delete single sql table
#'
#' @inheritParams readTable
deleteTable = function(tableName, uorfDB = get("uorfDB", envir = .GlobalEnv)){
  if (!tableNotExists(tableName)) { dbRemoveTable(uorfDB,tableName)
  } else { print(paste(tableName, "is not a table in the uORF database"))
    }
}

#' Delete uorf tables
#'
#' For rerunning
#' It keeps the RNA seq tables, since they usually not change
#' @param onlyRibo (FALSE), only ribo-seq / RNA seq features
#' @param onlySeq (FALSE), only sequence features
#' @export
deleteUorfTables <- function(onlyRibo = FALSE, onlySeq = FALSE) {
  # Ribo seq features
  if(onlyRibo | !onlySeq) {
    deleteTable("RSS")
    deleteTable("RRS")
    deleteTable("Ribofpkm")
    deleteTable("countRFP")
    deleteTable("fpkmRFP")
    deleteTable("startRegionRelative")
    deleteTable("startRegionCoverage")
    deleteTable("ORFScores")
    deleteTable("disengagementScores")
    deleteTable("ioScore")
    deleteTable("RiboByTissueMean")
    deleteTable("RiboByTissueTF")
    deleteTable("TEByTissueMeanWithInf")
    deleteTable("TEByTissueMeanWithoutInf")
    deleteTable("TEByTissueTF")
    deleteTable("teFiltered")
    deleteTable("teUnfiltered")
    deleteTable("tissueAtlasByCage")
    deleteTable("entropyRFP")
    deleteTable("floss")
    deleteTable("startCodonCoverage")
    if(!tableNotExists("biggestTEVariance")) {
      deleteTable("biggestTEVariance")
      deleteTable("smallestTEVariance")
    }
    if(onlyRibo) return(NULL)
  }

  # sequence features
  if(onlySeq | !onlyRibo) {
    deleteTable("numberOfUorfsPerTx")
    deleteTable("rankInTx")
    deleteTable("fractionLengths")
    deleteTable("kozak")
    deleteTable("inFrameCDS")
    deleteTable("isOverlappingCds")
    deleteTable("distORFCDS")
    deleteTable("distORFTSS")
    deleteTable("linkORFsToTx")
    deleteTable("isOverlappingCds")
    deleteTable("StopCodons")
    deleteTable("StartCodons")
    deleteTable("finalCAGEuORFPrediction")
    deleteTable("gcContent")
    deleteTable("goTerms")
    deleteTable("exon-exonJunctionsLeader")
    deleteTable("exon-exonJunctionsuORFs")
    deleteTable("uORFTxToGene")
    if(onlySeq) return(NULL)
  }

  deleteTable("uorfsAsGRWithTx")
  deleteTable("uniqueIDs")
  deleteTable("SplittedByExonsuniqueUORFs")
  deleteTable("toUniqueOrder")

  return(NULL)
}
