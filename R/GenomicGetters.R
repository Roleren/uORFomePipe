# These functions are used the speed up loading of annotation in parallel loops

#' Get the Genomic transcript format, currently using GRch38 data
getGTF <- function(assignIt = TRUE) {
  if(exists("Gtf", mode = "S4") == FALSE) {
    Gtf = loadTxdb(gtfdb)
    if(assignIt)
      assign("Gtf",Gtf,envir = .GlobalEnv)
  } else if(!dbIsValid(Gtf$conn)) {
    Gtf = loadTxdb(Gtf$conn@dbname)
    assign("Gtf", Gtf, envir = .GlobalEnv)
  }
}

#' Get transcripts from gtf
getTx <- function(assignIt = TRUE,
                  dataFolder = get("dataFolder", envir = .GlobalEnv)){
  if (exists("tx", mode = "S4") == FALSE) {
    if (file.exists(p(dataFolder, "/tx.rds"))) {
      tx <- readRDS(p(dataFolder, "/tx.rds"))
    } else {
      getGTF()
      tx <- loadRegion(Gtf)
    }
    if (assignIt) {
      assign("tx",tx, envir = .GlobalEnv)
      return(tx)
    } else {
      return(tx)
    }
  }
}

#Get the coding sequences from the gtf file
getCDS <- function(assignIt = TRUE,
                   dataFolder = get("dataFolder", envir = .GlobalEnv)) {
  if (exists("cds", mode = "S4") == FALSE) {
    if (file.exists(p(dataFolder, "/cds.rds"))) {
      cds <- readRDS(p(dataFolder, "/cds.rds"))
    } else {
      getGTF()
      cds <- loadRegion(Gtf, "cds")
    }
    if (assignIt) {
      assign("cds", cds, envir = .GlobalEnv)
      return(cds)
    } else {
      return(cds)
    }
  }
}

#Get the 3' sequences from the gtf file
getThreeUTRs <- function(dataFolder = get("dataFolder", envir = .GlobalEnv)){
  if (!exists("threeUTRs", mode = "S4")) {
    if (!file.exists(p(dataFolder, "/threeUTRs.rds"))) {
      getGTF()
      threeUTRs = threeUTRsByTranscript(Gtf, use.names = TRUE)
    } else {
      threeUTRs <- readRDS(p(dataFolder, "/threeUTRs.rds"))
    }
    assign("threeUTRs", threeUTRs, envir = .GlobalEnv)
  }
}

#' Get the 5' leaders
getLeaders <- function(assignIt = TRUE,
                       dataFolder = get("dataFolder", envir = .GlobalEnv)) {
  if(exists("fiveUTRs", mode = "S4") == FALSE) {
    if (file.exists(p(dataFolder,"/fiveUTRs.rds"))) {
      fiveUTRs <- readRDS(p(dataFolder,"/fiveUTRs.rds"))
    } else stop("Could not find pre saved leaders")
  }
  if(assignIt)
    assign("fiveUTRs", fiveUTRs, envir = .GlobalEnv)
  return(invisible(NULL))
}

leaderCage <- function(with.cds = TRUE){
  if(with.cds)
    return(readRDS(p(dataFolder,"/CageFiveUTRsWithCDS.rds")))
  return(readRDS(p(dataFolder,"/CageFiveUTRs.rds")))
}

#' Get the fasta indexed file
#'
#' if assignIt is TRUE, the object is not return to local scope
#' Only assigned to globalenvir
getFasta = function(filePath = NULL, assignIt = T){
  if(exists("fa", mode = "S4") == FALSE){ #index files
    if (is.null(filePath)){
      fa = FaFile(faiName)
    } else {
      fa = FaFile(filePath)
    }
    if (assignIt){
      assign("fa",fa,envir = .GlobalEnv)
    } else {
      return(fa)
    }
  }
}

#' Get all annotation parts
#' leader, cds, threeUTRs and tx
#' @param include.cage logical T
#' @param cdsOnFiveEnd logical F
#' @export
getAll <- function(include.cage = T, cdsOnFiveEnd = F){
  getFasta()
  getLeaders(); getCDS(); getThreeUTRs()
  #or with extension
  if (include.cage) {
    cageFiveUTRs <- leaderCage(cdsOnFiveEnd)
    assign("cageFiveUTRs", cageFiveUTRs,  envir = .GlobalEnv)
    getCageTx()
  } else getTx(T)

  return(invisible(NULL))
}

#' New leaders defined as maximum upstream for each original leader
#'
#' With these leader all found uORFs are guaranteed to overlap within
getCageTx <- function() {
  if (file.exists(p(dataFolder, "/cageTx.rdata"))) {
    load(p(dataFolder, "/cageTx.rdata"), envir = .GlobalEnv)
  } else {
    tx <- getTx()
    cageFiveUTRs <- leaderCage()
    tx[names(cageFiveUTRs)] <- ORFik:::extendLeaders(tx, cageFiveUTRs)
    assign("tx", tx,  envir = .GlobalEnv)
    save(tx, file = p(dataFolder, "/cageTx.rdata"))
  }
  return(NULL)
}

#' Used as negative set for training
#' Defined as stop site of 3' UTRs that has width >= 6,
#' and extend in direction of downstream direction the median
#' distance of 3' UTRs.
getSpecialThreeUTRs <- function() {
  getThreeUTRs()
  threeUTRs <- threeUTRs[widthPerGroup(threeUTRs) > 5]
  threeWidth <- median(widthPerGroup(threeUTRs))
  # TODO: think harder on if there is a better region
  downstreamThree <- GRanges(seqnamesPerGroup(threeUTRs, F),
                       IRanges(stopSites(threeUTRs, is.sorted = T) + 1, width = threeWidth),
                       strand = strandPerGroup(threeUTRs, F))
  names(downstreamThree) <- names(threeUTRs)
  return(groupGRangesBy(downstreamThree))
}