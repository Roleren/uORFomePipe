#' Get sequence features from orfik
#' @inheritParams checkAndInitPipe
#' @param grl the orfs as GRangesList, default uORFomePipe:::getUorfsInDb()
#' @return return(invisible(NULL)
getSequenceFeatures <- function(organism, biomart,
                                grl = uORFomePipe:::getUorfsInDb(),
                                mode = "uORF") {
  addit <- ifelse(mode == "CDS", "verify", "")
  if(!tableNotExists(p(addit, "kozak")) & !tableNotExists(p(addit,"stopCodonGrouping")) &
     file.exists("features/uORFSequenceFeatures.csv")) {
    print("sequence features exist, delete and run again if you want remake!")
    return(NULL)
  }
  if (!(mode %in% c("uORF", "CDS", "aCDS")))
    stop("mode must be uORF or CDS or aCDS (artificial CDS)")
  message("Creating sequence features from ORFs")
  uORFomePipe:::getAll(); uORFomePipe:::getGTF()
  # gene transcript connections for naming

  # kozak
  if (organism %in% c(c("human", "Homo sapiens"),
                      c("mouse", "Mus musculus"),
                      c("zebrafish", "Danio rerio"))) {
    kozak <- kozakSequenceScore(grl, tx, fa, species = organism)
  } else {
    message("Did not find kozak sequence for organism, checking vs human kozak")
    kozak <- kozakSequenceScore(grl, tx, fa)
  }

  # Start and stop
  starts <- startCodons(grl, is.sorted = T)
  StartCodons <- ORFik:::txSeqsFromFa(starts, fa, TRUE, FALSE)
  stops <- stopCodons(grl, is.sorted = TRUE)
  StopCodons <- ORFik:::txSeqsFromFa(stops, fa, TRUE, FALSE)

  # lengths
  lengths <- widthPerGroup(grl, F)
  if (mode == "uORF") {
    dt <- data.table(txNames = txNames(grl), geneNames = ORFik:::txNamesToGeneNames(txNames(grl), Gtf))
    distORFCDS <- distToCds(grl, cageFiveUTRs, cds)
    distORFTSS <- distToTSS(grl, cageFiveUTRs)
    fractionLengths <- fractionLength(grl, ORFik:::widthPerGroup(tx))
    inFrameCDS <- ORFik:::isInFrame(distORFCDS)
    isOverlappingCds <- isOverlapping(distORFCDS)
    rankInTx <- rankOrder(grl)
    exonExonJunctions <- numExonsPerGroup(grl, FALSE)
    gcContent <- gcContent(grl, fa)
    if (!is.null(biomart)) {
      goTerms <- uORFomePipe:::getORFsGoTerms(dt$geneNames, organism)
    } else {
      goTerms <- rep("NA", length(lengths))
    }


    seqData <- data.table(rankInTx, lengths, kozak, fractionLengths,
                          StartCodons = as.factor(StartCodons),
                          StopCodons = as.factor(StopCodons),
                          distORFCDS, distORFTSS, inFrameCDS, isOverlappingCds,
                          exonExonJunctions, gcContent, goTerms = as.factor(goTerms))



    # Gene information
    insertTable(dt, "uORFTxToGene")
    # Gene to symbol
    if (!is.null(biomart))
      insertTable(getAllORFGeneSymbols(dt$geneNames, biomart), "geneSymbols")

    # exon-exon junctions
    eej <- numExonsPerGroup(fiveUTRs, TRUE)
    txNames <- txNames(grl)
    eej <- as.integer(eej[txNames])
    insertTable(data.table(eej = eej), "exon-exonJunctionsLeader")

    # number of uorfs per tx
    numberOfUorfsPerTx <- S4Vectors::Rle(txNames)
    insertTable(data.table(nUorfs = runLength(numberOfUorfsPerTx)), "numberOfUorfsPerTx")
  } else if (mode != "ORF") {
    #distORFCDS <- -widthPerGroup(grl, FALSE)
    seqData <- data.table(lengths, kozak,
                          StartCodons = as.factor(StartCodons),
                          StopCodons = as.factor(StopCodons))
  }
  fwrite(seqData, file = p("features/", addit, "uORFSequenceFeatures.csv"))
  # Create one table per feature in DB

  for(f in colnames(seqData)) {
    featu <-  seqData[, which(colnames(seqData) == f), with = FALSE]
    insertTable(featu, p(addit, f))
  }
  # Stop codon grouping
  insertTable(data.table(stopCodonGrouping = uniqueOrder(stops)), p(addit, "stopCodonGrouping"))
  return(invisible(NULL))
}

#' Get Ribo-seq features
#'
#' Excluding features that uses RNA-seq normalizations if df.rna is NULL
#' @param df.rna NULL
#' @param grl the ORFs as GRangesList
#' @param preName name to pre append in database for each feature
#' @param threeUTRsSpecial default NULL, or a GRangesList of 3' UTRs
#' @inheritParams find_uORFome
#' @return NULL (features saved to database)
getGeneralRiboFeatures <- function(df.rfp, df.rna = NULL,
                                   grl, preName = "",
                                   threeUTRsSpecial = NULL,
                                   BPPARAM = bpparam()) {
  pr <- ifelse(preName == "", "uORFs", preName)
  if (tableNotExists(p(preName, "ioScore"))) {
    message("--------------------------------------")
    message(paste("Finding Ribo-seq features for", pr))
    uORFomePipe:::getAll()
    if (!is.null(threeUTRsSpecial)) threeUTRs <- threeUTRsSpecial
    startRegion <- startRegion(grl, tx, T, -3, 9)

    paths <- filepath(df.rfp, "pshifted")
    paths.rna <- NULL
    if (!is.null(df.rna)) paths.rna <- filepath(df.rna, "ofst")
    if (!is.null(df.rna) & length(grep("\\.ofst$", paths.rna)) == 0) paths.rna <- filepath(df.rna, "bedo")


    libs <-bplapply(seq(length(paths)),
      function(x, grl, fiveUTRs, threeUTRs, cds, startRegion, paths, paths.rna) {
        message(paths[[x]])
        X <- paths[[x]]
        Y <- paths.rna[[x]]
        style <- seqlevelsStyle(grl)
        RNA <- if (!is.null(Y)) {
          fimport(Y, style)
        } else NULL
        weight.RNA <- 1L
        if (!is.null(RNA)) {
          if (!is.null(mcols(RNA)$score)) weight.RNA <- "score"
        }


        ORFik:::allFeaturesHelper(grl, RFP = fimport(X, style), RNA = RNA, tx, fiveUTRs, cds ,
                                  threeUTRs,
                                  faFile = NULL, riboStart = 26, riboStop = 34,
                                  sequenceFeatures = FALSE, grl.is.sorted = TRUE,
                                  weight.RFP = "score", weight.RNA = weight.RNA,
                                  st = startRegion)
      }, grl = grl, fiveUTRs = fiveUTRs, threeUTRs = threeUTRs,
      cds = cds, startRegion = startRegion, paths = paths, paths.rna = paths.rna,
      BPPARAM = BPPARAM)

    allRiboFeatures <- setDT(unlist(libs, recursive = FALSE))
    for(f in unique(colnames(allRiboFeatures))) { # Create one table per feature in DB
      featu <-  allRiboFeatures[, which(colnames(allRiboFeatures) == f), with = FALSE]
      colnames(featu) <- paste0(f, "_", 1:ncol(featu))
      insertTable(data.table(featu), p(preName, f), rmOld = TRUE)
    }
  } else {
    print(paste("Ribo-seq features exists in DB for:", pr, "delete and run again if you want new"))
  }
}
