#' Get sequence features from orfik
#'@param organism scientifi name, example "Homo sapiens" or "Danio rerio"
#'@param biomart name of biomart for organism, example "hsapiens_gene_ensembl"
#'or "drerio_gene_ensembl"
getSequenceFeatures <- function(organism, biomart) {
  if(!tableNotExists("kozak") & !tableNotExists("gcContent")) {
    print("sequence features exist, delete and run again if you want remake!")
    return(NULL)
  }
  message("Creating sequence features from uORFs")
  grl <- uORFomePipe:::getUorfsInDb()
  uORFomePipe:::getAll(); uORFomePipe:::getGTF()
  # gene transcript connections for naming
  dt <- data.table(txNames = txNames(grl), geneNames = ORFik:::txNamesToGeneNames(txNames(grl), Gtf))
  insertTable(dt, "uORFTxToGene")

  # kozak
  kozak <- kozakSequenceScore(grl, tx, fa)
  # lengths
  lengths <- widthPerGroup(grl, F)
  distORFCDS <- distToCds(grl, cageFiveUTRs, cds)
  distORFTSS <- distToTSS(grl, cageFiveUTRs)
  fractionLengths <- fractionLength(grl, ORFik:::widthPerGroup(tx))
  inFrameCDS <- ORFik:::isInFrame(distORFCDS)
  isOverlappingCds <- isOverlapping(distORFCDS)
  rankInTx <- rankOrder(grl)
  starts <- startCodons(grl, is.sorted = T)
  StartCodons <- ORFik:::txSeqsFromFa(starts, fa, TRUE, FALSE)
  stops <- stopCodons(grl, is.sorted = TRUE)
  StopCodons <- ORFik:::txSeqsFromFa(stops, fa, TRUE, FALSE)
  exonExonJunctions <- numExonsPerGroup(grl, FALSE)
  gcContent <- gcContent(grl, fa)
  goTerms <- uORFomePipe:::getORFsGoTerms(dt$geneNames, organism)

  seqData <- data.table(rankInTx, lengths, kozak, fractionLengths,
                         StartCodons = as.factor(StartCodons),
                         StopCodons = as.factor(StopCodons),
                         distORFCDS, distORFTSS, inFrameCDS, isOverlappingCds,
                         exonExonJunctions, gcContent, goTerms = as.factor(goTerms))
  fwrite(seqData, file = "features/uORFSequenceFeatures.csv")
  for(f in colnames(seqData)) { # Create one table per feature in DB
    featu <-  seqData[, which(colnames(seqData) == f), with = FALSE]
    insertTable(featu, f)
  }

  # Gene information
  # Gene to symbol
  insertTable(getAllORFGeneSymbols(dt$geneNames, biomart), "geneSymbols")

  # exon-exon junctions
  eej <- numExonsPerGroup(fiveUTRs, TRUE)
  link <- readTable("linkORFsToTx")
  eej <- as.integer(eej[link$txNames])
  insertTable(data.table(eej = eej), "exon-exonJunctionsLeader")
  # Stop codon grouping
  insertTable(data.table(stopCodonGrouping = uniqueOrder(stops)), "stopCodonGrouping")
  # number of uorfs per tx
  txNames <- txNames(grl)
  numberOfUorfsPerTx <- S4Vectors::Rle(txNames)
  insertTable(data.table(nUorfs = runLength(numberOfUorfsPerTx)), "numberOfUorfsPerTx")

  return(invisible(NULL))
}

#' Get Ribo-seq features
#'
#' Excluding features that uses RNA-seq normalizations if df.rna is NULL
#' @param df.rfp ORFik experiment
#' @param df.rna NULL
#' @param grl the ORFs as GRangesList
#' @param preName name to pre append in database for each feature
#' @param threeUTRsSpecial default NULL, or a GRangesList of 3' UTRs
#' @return NULL (features saved to database)
getGeneralRiboFeatures <- function(df.rfp, df.rna = NULL,
                                   grl, preName = "",
                                   threeUTRsSpecial = NULL) {
  pr <- ifelse(preName == "", "uORFs", preName)
  if (tableNotExists(p(preName, "ioScore"))) {
    message(paste("Finding Ribo-seq features from", pr))
    uORFomePipe:::getAll()
    if (!is.null(threeUTRsSpecial)) threeUTRs <- threeUTRsSpecial
    startRegion <- startRegion(grl, tx, T, -3, 9)

    paths <- filepath(df.rfp, "pshifted")
    paths.rna <- NULL
    if (!is.null(df.rna)) paths.rna <- filepath(df.rna, "bedo")

    libs <-bplapply(seq(length(paths)),
      function(x, grl, fiveUTRs, threeUTRs, cds, startRegion, paths, paths.rna) {
        X <- paths[[x]]
        Y <- paths.rna[[x]]
        style <- seqlevelsStyle(grl)
        RNA <- if (!is.null(Y)) {
          fimport(Y, style)
        } else NULL
        ORFik:::allFeaturesHelper(grl, RFP = fimport(X, style), RNA = RNA, tx, fiveUTRs, cds ,
                                  threeUTRs,
                                  faFile = NULL, riboStart = 26, riboStop = 34,
                                  sequenceFeatures = FALSE, grl.is.sorted = TRUE,
                                  weight.RFP = "score", weight.RNA = 1L,
                                  st = startRegion)
      }, grl = grl, fiveUTRs = fiveUTRs, threeUTRs = threeUTRs,
      cds = cds, startRegion = startRegion, paths = paths, paths.rna = paths.rna)

    allRiboFeatures <- setDT(unlist(libs, recursive = FALSE))
    for(f in unique(colnames(allRiboFeatures))) { # Create one table per feature in DB
      featu <-  allRiboFeatures[, which(colnames(allRiboFeatures) == f), with = FALSE]
      colnames(featu) <- paste0(f, "_", 1:ncol(featu))
      insertTable(data.table(featu), p(preName, f), rmOld = TRUE)
    }
  } else {
    print(paste0("Ribo-seq features exists in DB for", pr, "delete and run again if you want new"))
  }
}