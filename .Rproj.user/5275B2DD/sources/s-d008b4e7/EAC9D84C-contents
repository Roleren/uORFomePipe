
getLeadersFromCage <- function(nCageList){
  export.list <- c("cageFolder","cageFiles", "regionUORFsFolder",
                   "getLeaders", "p", "dataFolder", "getGTF",
                   "getCDS", "leadersFolder")
  foreach(i=1:nCageList, .inorder = F, .export = export.list,
          .packages = c("ORFik")) %dopar% {

    getLeaders(cageName = p(cageFolder, cageFiles[i]), assignLeader = F , exportUorfRegions = T)
    print("ok")
  }
}

getUorfsFromLeadersNew <- function(folder = regionUORFsFolder,
                                   startCodons = "ATG|CTG|TTG|GTG|AAG|AGG|ACG|ATC|ATA|ATT") {
  leadersList = list.files(folder)
  nLeadersList = length(leadersList)
  export.list <- c("uorfFolder","folder","leadersList","dataFolder","getCDS", "faiName","getFasta", "filterORFs", "p")
  foreach(i=1:nLeadersList, .inorder = F, .export = export.list, .packages = "ORFik") %dopar% {
    saveName = p(uorfFolder, gsub(pattern = "regionUORF.rdata", replacement = "uorf.rdata",
                                  x = leadersList[i]))
    if (!file.exists(saveName)) {
      load(p(folder, leadersList[i]))
      getFasta()
      getCDS()
      rangesOfuORFs <- findUORFs(uORFSeachRegion, fa, startCodon = "ATG|CTG|TTG|GTG|AAG|AGG|ACG|ATC|ATA|ATT",
                                 minimumLength = 0, longestORF = FALSE)
      rangesOfuORFs <- ORFik:::filterUORFs(rangesOfuORFs, get(cds, mode = "S4"))
      save(rangesOfuORFs, file = saveName)
      return(i)
    }
    return(0)
  }
}

getIDsFromUorfs <- function(folder = uorfFolder){
  uorfFiles = list.files(folder)
  nuorfsList <- length(uorfFiles)

  foreach(i=1:nuorfsList, .inorder = F, .export = c("idFolder","folder")) %dopar% {
    load(paste0(folder, list.files(folder)[i]))

    uorfID <- unique(ORFik:::orfID(rangesOfuORFs))
    saveName = paste0(idFolder, gsub("uorf.rdata","",list.files(folder)[i]),"uorfID.rdata")
    save(uorfID, file = saveName)
    print("ok")
  }
}

getAllFeaturesFromUorfs <- function() {
  if (tableNotExists("ioScore")) {
    # if(length(RFPPath) != 1) stop(paste("did not find unique RFP file for:", matching_rna_ribo$study[i],matching_rna_ribo$ribo[i]))
    getAll()
    getGeneralRiboFeatures(grl = getUorfsInDb(), cds = cds, threeUTRs = threeUTRs,
                           tx = tx, df.rfp = df.rfp)


  } else {
    print("AllFeaturesFromUorfs exists in DB (ioScore), delete and run again if you want new")
  }
}

#' Make RNA-seq fpkm values for database
#'
#' Since rna-seq fpkms are normalized by transcript, the size of this
#' tabke might be different than the ribo-seq. i.g. there can be several
#' uORFs per transcript.
getRNAFpkms <- function() {
  if (!tableNotExists("RNAfpkm")) {
    print("RNAfpkm table exists, remove it if you want to rerun with new values")
    return(NULL)
  }

  rnaFPKMs <- foreach(RNAPath = df.rna$filepath, .combine = 'cbind', .export = c("getCageTx"),
                      .packages = c("GenomicFeatures", "ORFik")) %dopar% {
    getCageTx()
    RNA <- fimport(RNAPath)
    rnaFPKM <- ORFik:::fpkm(tx, reads = RNA)
  }
  getCageTx()
  rnaFPKMs <- as.data.table(rnaFPKMs)
  txNames <- names(tx)
  rnaFPKMs <- data.table(txNames, rnaFPKMs)

  insertTable(rnaFPKMs, "RNAfpkm", rmOld = T)
  return(NULL)
}

#' This function uses the fact that 1st col of ribo is connected to 1st col of RNA.
getTeFeatures <- function(riboDbName = "Ribofpkm",
                          dbOutputNames = c("teUnfiltered", "teFiltered")){
  if(length(dbOutputNames) != 2) stop("dbOutputNames must have 2 character elements")

  # load linking and ribo / rna
  RFP <- readTable(riboDbName)
  RNA <- readTable("RNAfpkm")

  RNA <- matchByTranscript(RNA, RFP)
  RNA <- removeIDColumns(RNA)
  RFP <- RFP[,getRiboMatchedToAll(), with = F]

  if(nrow(RFP) != nrow(RNA)) stop("riboseq and rnaseq tables have different # of rows")
  if(ncol(RFP) != ncol(RNA)) stop("riboseq and rnaseq tables have different # of cols")
  # find number of linkings we have
  nTE <- ncol(RFP)

  # unfiltered without pseudoCounts
  teTable <- foreach(i = 1:nTE, .combine = 'cbind') %do% {
    print(i)
    return(RFP[,i, with = F] / RNA[,i, with = F])
  }

  teTable <- data.table(txNames, teTable)
  insertTable(teTable, dbOutputNames[1], rmOld = T)

  # filtered with pseudoCounts

  teTable <- foreach(i = 1:nTE, .combine = 'cbind') %do% {
    print(i)
    return((RFP[,i, with = F] + 1) / (RNA[,i, with = F] + 1))
  }
  teTable <- data.table(txNames, teTable)
  insertTable(teTable, dbOutputNames[2], rmOld = T)
}

getCDSFeatures <- function(){
  if(tableNotExists("cdsKozak")) {
    getCDS()
    getFasta()
    getCageTx()
    cdsKozak <- ORFik:::kozakSequenceScore(cds, tx, fa)
    cdsKozak <- data.table(txNames = names(cds), cdsKozak)
    setwd(codeFolder)
    insertTable(cdsKozak, "cdsKozak")
  }

  if(tableNotExists("cdsIos")) {
    rfpList <- grep(pattern = "merged",x = list.files(rfpFolder), value = T)
    nrfpList <- length(rfpList)
    getCDS()
    # gr <- unlist(cds, use.names = T)
    # cds <- groupGRangesBy(gr)
    getThreeUTRs()
    # gr <- unlist(threeUTRs, use.names = T)
    # threeUTRs <- groupGRangesBy(gr)
    rm(tx)
    getCageTx()
    cageFiveUTRs <- leaderCage()

    getGeneralRiboFeatures(grl = cds, cds = cds, threeUTRs = threeUTRs,
                           tx = tx, name = "three")
  }

  if(tableNotExists("cdsFractionLengths")) {
    getCDS()
    tx <- getTx()
    tx <- tx[names(cds)]
    cageFiveUTRs <- leaderCage()
    tx[names(cageFiveUTRs)] <- ORFik:::extendLeaders(tx[names(cageFiveUTRs)], cageFiveUTRs)
    cdsFrac <- fractionLength(cds, widthPerGroup(tx))

    cdsFrac <- data.table(txNames = names(cds), cdsFrac)
    insertTable(cdsFrac, "cdsFractionLengths")
  }

  ## RNA fpkms already made, so go to TE:
  if(tableNotExists("cdsTEFiltered")) {
    getTeFeatures(riboDbName = "cdsRfpFPKMs",
                  dbOutputNames = c("cdsTEUnfiltered", "cdsTEFiltered"))
    cdsTEs <- readTable("cdsTEFiltered", with.IDs = T)
    teAtlasTissue(cdsTEs, "cdsTETissueMean")
  }
}

#' downstream of three utrs is used as new 3'utre
getFeaturesThreeUTRs <- function(){
  if(tableNotExists("threeIos")) {

    getCDS()
    # gr <- unlist(cds, use.names = T)
    # cds <- groupGRangesBy(gr)
    getThreeUTRs()
    # gr <- unlist(threeUTRs, use.names = T)
    # threeUTRs <- groupGRangesBy(gr)
    threeUTRs <- threeUTRs[widthPerGroup(threeUTRs) > 5]
    threeWidth <- median(widthPerGroup(threeUTRs))
    # this is a crap region to use ->
    fakeThree <- GRanges(seqnamesPerGroup(threeUTRs, F),
                         IRanges(stopSites(threeUTRs, is.sorted = T) + 1, width = threeWidth),
                         strand = strandPerGroup(threeUTRs, F))
    names(fakeThree) <- names(threeUTRs)
    fakeThree <- groupGRangesBy(fakeThree)
    rm(tx)
    getCageTx()
    cageFiveUTRs <- leaderCage()


    getGeneralRiboFeatures(grl = threeUTRs, cds = cds, threeUTRs = fakeThree,
                           tx = tx, name = "three")
  }
  # sequence features
  if (tableNotExists("threeKozak")) {
    getThreeUTRs()
    gr <- unlistGrl(threeUTRs)
    threeUTRs <- groupGRangesBy(gr)
    threeUTRs <- threeUTRs[widthPerGroup(threeUTRs) > 5]
    getFasta()
    getCageTx()
    Kozak <- kozakSequenceScore(threeUTRs, tx, fa)
    insertTable(Kozak, "threeKozak")
  }
}


