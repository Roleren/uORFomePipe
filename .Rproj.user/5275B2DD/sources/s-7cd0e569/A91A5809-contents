#' Make a leader that spans the max leader per transcript by cage files
allLeadersSpanningLeader <- function(){

  leadersList = list.files(leadersFolder)
  nLeadersList = length(leadersList)

  # Get all width of cage experiment TSS reassignments
  widths <- foreach(i=1:nLeadersList, .combine = 'cbind') %dopar% {
    # TODO: Fix this path
    setwd("/export/valenfs/projects/uORFome/RCode1/")
    source("./HelperVariables.R")
    leadersList = list.files(leadersFolder)
    load(p(leadersFolder,leadersList[i]))

    return(widthPerGroup(fiveUTRs, keep.names = F))
  }

  maxWidths <- rowMaxs(widths)
  getCDS()
  if( exists("fiveUTRs", mode = "S4")) {
    rm(fiveUTRs)
  }
  getLeaders()

  change <- maxWidths - widthPerGroup(fiveUTRs, F)
  newStarts <- rep.int(0L, length(fiveUTRs))
  outsidePos <- strandBool(fiveUTRs) & (change >= 0)
  outsideMin <- !strandBool(fiveUTRs) & (change >= 0)
  either <- outsidePos | outsideMin
  # Outside Leader
  # pos
  newStarts[outsidePos] <- ORFik:::startSites(fiveUTRs)[outsidePos] -
    change[outsidePos]
  # min
  newStarts[outsideMin] <- ORFik:::startSites(fiveUTRs)[outsideMin] +
    change[outsideMin]
  fOut <- fiveUTRs[either]
  fOut <- ORFik:::downstreamFromPerGroup(fOut, newStarts[either])

  # Inside Leader
  inside <- change < 0
  if(any(inside)) {
    fIn <- pmapToTranscripts(fiveUTRs[inside], fiveUTRs[inside])
    start(fIn) <- start(fIn) - change[inside] # -- = +
    fIn <- unlist(fIn, use.names = F)
    fIn <- ranges(fIn)
    names(fIn) <- which(inside)
    fInNew <- ORFik:::pmapFromTranscriptF(fIn, fiveUTRs)
  } else {
    fInNew <- GRangesList()
  }


  fTot <- fiveUTRs
  fTot[either] <- fOut
  fTot[inside] <- fInNew
  if(!all(ORFik:::widthPerGroup(fTot, F) == maxWidths)) {
    stop("Algorithm is wrong for five extension!")
  }


  # all ok, then save
  setwd("/export/valenfs/projects/uORFome/RCode1/")
  source("./DataBaseSetup.R")
  save(widths, file = "./leaderLengths.rdata")
  CageFiveUTRs <- fTot
  save(CageFiveUTRs, file = "CageFiveUTRs.rdata")
  CageFiveWithCDS <- ORFik:::addCdsOnLeaderEnds(fTot, cds)
  save(CageFiveWithCDS, file = "CageFiveUTRsWithCDS.rdata")
  return(NULL)
}

#' Assign transcriptnames to orfs, and find rank for each orf
#'
#' Given orfs and transcripts, find all transcripts the orfs are within
#' and name them by this. Also the second orf in
linkORFsToTx <- function(){
  if (!file.exists(p(dataBaseFolder,"/uniqueUorfsAsGRWithTx.rdata"))) {
    leaders <- leaderCage()
    grl <- getUorfsInDb(T, F, F)
    overlaps <- findOverlaps(grl, leaders, type = "within")
    if( length(unique(from(overlaps))) != length(grl)) {
      stop("leader is not spanning all uORFs, check for uORFs going into cds.")
    }

    sortedIndeces <- order(to(overlaps))
    from <- from(overlaps)[sortedIndeces]
    to <- to(overlaps)[sortedIndeces]
    txNames <- names(leaders)[to]
    uorfIDs <- ORFik:::orfID(grl)[from]
    dt <- data.table(uorfID = uorfIDs, txNames = txNames)
    insertTable(Matrix = dt, tableName = "linkORFsToTx",rmOld = T)

    # now make grl with transcript mapping
    grlb <- grl[from]
    names(grlb@unlistData) <- NULL
    names(grlb) <- txNames
    asGR <- unlist(grlb, use.names = T)
    names(grlb@unlistData) <- names(asGR)

    grl <- ORFik:::makeORFNames(grlb, F)
    save(grl, file = "uoRFsAsGRAllWithTx.rdata")
    # Unique uORFs
    c <- ORFik:::orfID(grl, with.tx = F)
    d <- which(!duplicated(c))
    grl <- grl[d]
    insertTable(Matrix = data.table(d), tableName = "toUniqueOrder", rmOld = T)
    insertTable(Matrix = dt[d,], tableName = "linkORFsToTxUnique", rmOld = T)
    save(grl, file = "uniqueUorfsAsGRWithTx.rdata")

  } else {
    message("linkORFsToTx already exist, skipping remake of them")
  }
  return(NULL)
}
