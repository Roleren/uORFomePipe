
#' Get orf names from the orf data-base
#'
#' This is the primary key for most tables in the data-base
#' @param with.transcript a logical(F), should the transcript be included, this makes the
#'  list have duplicated orfs
#' @param only.transcripts a logical(F), should only the transcript and not orfId be included
#' @param asCharacter a logical(T), should it return as character or data.table(F)
getORFNamesDB <- function(with.transcript = F, only.transcripts = F, asCharacter = T, uniques = F){
  if (with.transcript) {
    if(uniques) {
      dt <- readTable("linkORFsToTxUnique")
    } else {
      dt <- readTable("linkORFsToTx")
    }
    if(only.transcripts){
      if (asCharacter) {
        return(as.character(unlist(dt[, 2], use.names = F)))
      }
      return(dt[, 2])
    }
    return(dt)
  }

  dt <-readTable("uniqueIDs")
  if (asCharacter) {
    dt <- as.character(unlist(dt, use.names = F))
  }
 return(dt)
}

#' Takes two tables from the database and extracts the rows of toBeMatched
#' that matches the txNames in referenced.
#' Both must have a column called txNames
#' @return the toBeMatched object matched by txNames
matchByTranscript <- function(toBeMatched, referenced){

  Indices <- data.table(txNames = toBeMatched$txNames, ind = 1:length(toBeMatched$txNames))
  merged <- merge(Indices, data.table(txNames = referenced$txNames),
                  by = "txNames", all.y = T, sort = F)
  return(toBeMatched[merged$ind, ])
}

getIDColumns <- function(dt, allowNull = F){
  nIDs <- 0
  if (!is.numeric(dt[1,1][[1]])) {
    nIDs = nIDs + 1
    if (!is.numeric(dt[1,2][[1]])) {
      nIDs = nIDs + 1
    }
  }
  if(!nIDs){
    if (allowNull) {
      return(NULL)
    } else {
      stop("No id columns found for dt")
    }
  }
  return(dt[, nIDs, with = FALSE])
}
#' fix this to work on string tables
removeIDColumns <- function(dt){
  if (!is.numeric(dt[1,1][[1]])) {
    dt <- dt[, -1]
    if (!is.numeric(dt[1,1][[1]])) {
      dt <- dt[, -1]
    }
  }
  return(dt)
}

#' get the uorfs in the database
#' @param withExons logical, default TRUE, should the uorfs be splitted by exons
#' @param withTranscripts logical, default TRUE, should the uorfs have transcript information,
#' warning, this will duplicate some uorfs.
#' @param uniqueORFs logical, default TRUE, only unique
#' @param tissue character, only applicable when mode is "aCDS"
#' @return a GRangesList or data.table, if(F, F)
#' @export
getUorfsInDb <- function(withExons = T, withTranscripts = T, uniqueORFs = T,
                         mode = "uORF") {
  dataFolder <- get("dataFolder", envir = .GlobalEnv)
  if (withExons && withTranscripts) {
    if(uniqueORFs) {
      file <-  if (mode == "CDS") {
        p(dataFolder, "/uniqueUorfsAsGRWithTx_verify", ".rdata")
      } else {
        file <- p(dataFolder, "/uniqueUorfsAsGRWithTx", ".rdata")
      }
      if (file.exists(file)) {
        load(file)
        return(grl)
      } else stop("unique uorfs with tx does not exists")
    }
    if(file.exists(p(dataFolder, "/uoRFsAsGRAllWithTx.rdata"))) {
      load(p(dataFolder, "/uoRFsAsGRAllWithTx.rdata"))
      return(grl)
    } else if(!tableNotExists("uorfsAsGRWithTx")) {
      grl <- readTable("uorfsAsGRWithTx", asGR = T)
      gr <- unlist(grl, use.names = F)
      names(gr) <- gsub("_[0-9]*", "", names(gr))
      return(groupGRangesBy(gr, gr$names))
    }
    stop("uORFs could not be found, check that they exist")

  } else if (!withExons) {
    return(readTable("uniqueIDs"))
  } else if (withExons && !withTranscripts) {
    if(uniqueORFs) {
      if (file.exists(p(dataFolder, "/uniqueUorfsAsGR.rdata"))) {
        load(p(dataFolder, "/uniqueUorfsAsGR.rdata"))
        return(grl)
      }
    }
    return(readTable("SplittedByExonsuniqueUORFs", asGR = T))
  } else {
    stop("not supported way of getting uorfs")
  }
}
