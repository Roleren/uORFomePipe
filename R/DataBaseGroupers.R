#' Make a leader that spans the max leader per transcript by cage files
#' @param leadersFolder folder with CAGE leaders
#' @param dataFolder dataFolder with helper objects
#' @importFrom BiocParallel bplapply
#' @importFrom GenomicRanges start
#' @return return(invisible(NULL))
allLeadersSpanningLeader <- function(leadersFolder, dataFolder){
  message("Creating merged leader")
  leadersList = list.files(leadersFolder, full.names = TRUE)

  # Get all width of cage experiment TSS reassignments
  widths <- bplapply(leadersList, function(i) {
    return(widthPerGroup(readRDS(i), keep.names = FALSE))
  })
  widths <- matrix(unlist(widths), ncol = length(widths))

  maxWidths <- rowMaxs(widths)

  getLeaders()
  change <- maxWidths - widthPerGroup(fiveUTRs, FALSE)
  newStarts <- rep.int(0L, length(fiveUTRs))
  outsidePos <- strandBool(fiveUTRs) & (change >= 0)
  outsideMin <- !strandBool(fiveUTRs) & (change >= 0)
  either <- outsidePos | outsideMin
  # Outside Leader
  # pos
  newStarts[outsidePos] <- startSites(fiveUTRs)[outsidePos] -
    change[outsidePos]
  # min
  newStarts[outsideMin] <- ORFik:::startSites(fiveUTRs)[outsideMin] +
    change[outsideMin]
  fOut <- fiveUTRs[either]
  fOut <- ORFik:::downstreamFromPerGroup(fOut, newStarts[either])

  # Inside Leader
  inside <- change < 0
  if(any(inside)) {
    fIn <- ORFik:::pmapToTranscriptF(fiveUTRs[inside], fiveUTRs[inside])
    start(fIn) <- start(fIn) - change[inside] # -- = +
    fIn <- unlist(fIn, use.names = FALSE)
    fIn <- ranges(fIn)
    names(fIn) <- which(inside)
    fInNew <- pmapFromTranscriptF(fIn, fiveUTRs, removeEmpty = TRUE)
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
  saveRDS(widths, file = p(dataFolder, "/leaderLengths.rds"))
  CageFiveUTRs <- fTot
  saveRDS(CageFiveUTRs, file = p(dataFolder, "/CageFiveUTRs.rds"))
  getCDS()
  CageFiveWithCDS <- ORFik:::addCdsOnLeaderEnds(fTot, cds)
  saveRDS(CageFiveWithCDS, file = p(dataFolder, "/CageFiveUTRsWithCDS.rds"))
  return(invisible(NULL))
}


#' Assign transcriptnames to orfs, and find rank for each orf
#'
#' Given orfs and transcripts, find all transcripts the orfs are within
#' and name them by this. Also the second orf in
#' @param dataFolder dataFolder with helper objects
#' @return return(invisible(NULL))
linkORFsToTx <- function(dataFolder){
  if (!file.exists(p(dataFolder,"/uniqueUorfsAsGRWithTx.rdata"))) {
    message("Linking ORFs to transcripts")
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
    save(grl, file = p(dataFolder,"/uoRFsAsGRAllWithTx.rdata"))
    # Unique uORFs
    c <- ORFik:::orfID(grl, with.tx = F)
    d <- which(!duplicated(c))
    grl <- grl[d]
    insertTable(Matrix = data.table(d), tableName = "toUniqueOrder", rmOld = T)
    insertTable(Matrix = dt[d,], tableName = "linkORFsToTxUnique", rmOld = T)
    save(grl, file = p(dataFolder,"/uniqueUorfsAsGRWithTx.rdata"))

  } else {
    message("linkORFsToTx already exist, skipping remake of them")
  }
  return(invisible(NULL))
}
