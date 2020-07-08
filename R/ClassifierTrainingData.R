#' Get CDS and 3'UTR TrainingData of ribo seq features
#'
#' Positive set is cds, negative is downstream region of 3' UTRs
#' @param tissue Tissue to train on, use "combined" if you want all in one,
#' first run of training it is a ORFik experiment.
#' @param features features to train model on, must exist in database, default:
#' c("countRFP", "disengagementScores", "entropyRFP", "floss",
#' "fpkmRFP","ioScore", "ORFScores", "RRS", "RSS", "startCodonCoverage",
#' "startRegionCoverage","startRegionRelative")
#' @inheritParams makeTrainingAndPredictionData
#' @return invisible(NULL), saved to disc
makeTrainingData <- function(tissues = "combined",
                             features = c("countRFP", "disengagementScores", "entropyRFP", "floss",
                                          "fpkmRFP","ioScore", "ORFScores", "RRS", "RSS", "startCodonCoverage",
                                          "startRegionCoverage","startRegionRelative"),
                             mode = "uORF", max.artificial.length,
                             dataFolder = get("dataFolder", envir = .GlobalEnv)) {

  cds_file <- p(dataFolder, "/uniqueUorfsAsGRWithTx_verify",".rdata")
  if (mode == "CDS" & !file.exists(cds_file)) { # If CDS
    cds_total <- getCDSFiltered()
    cds_training <- getCDSTraining(mode = mode)
    grl <- cds_total[!(names(cds_total) %in% names(cds_training))]
    save(grl, file = cds_file)
  }

  # else for uORF and aCDS (artificial CDS)
  if (length(tissues) == 1) {
    if (file.exists(paste0("features/TrainingData_",tissues,".rds"))) {
      return(readRDS(paste0("features/TrainingData_",tissues,".rds")))
    }
  } else if (is(tissues, "experiment")) tissues <- readTable("experiment_groups")[[1]]

  for (tissue in tissues) {
    if (file.exists(paste0("features/TrainingData_",tissue,".rds"))) next
    posFeatureNames <- grep(pattern = p("cds", features, collapse = "|"), x = listTables(), value = TRUE)
    negFeatureNames <- grep(pattern = p("three", features, collapse = "|"), x = listTables(), value = TRUE)
    if (length(posFeatureNames) != length(negFeatureNames)) stop("Not equal length of pos and neg feature names")
    message("--------------------------------------")

    pos <-bplapply(posFeatureNames, function(x, tissue) {
      uORFomePipe:::getTissueFromFeatureTable(x, tissue = tissue)
    }, tissue = tissue); pos <- as.matrix(setDT(unlist(pos, recursive = FALSE)))

    neg <-bplapply(negFeatureNames, function(x) {
      uORFomePipe:::getTissueFromFeatureTable(tableName = x, tissue = tissue)
    }); neg <- as.matrix(setDT(unlist(neg, recursive = FALSE)))

    # Filter out translating 3' UTRs
    weakThree <- uORFomePipe:::weakTrailer(coverage = neg[,1], fpkm = neg[,5],
                                           startRegionRelative = neg[,12],
                                           ORFScore = neg[,7])
    # Filter out non-translating CDS'
    strongCDS <- uORFomePipe:::strongCDS(coverage = pos[,1], fpkm = pos[,5],
                                         startCodonCoverage = pos[,11],
                                         startRegionRelative = pos[,12],
                                         ORFScore = pos[,7])
    if (sum(strongCDS) < 30) stop("Training on less than 30 valid 'active' CDS is not allowed!")
    # If artificial uORFs split set in two:
    indices <- seq.int(length(strongCDS))
    if (mode == "aCDS" & !file.exists(p(dataFolder, "/uniqueUorfsAsGRWithTx",".rdata"))) {
      valid <- which(strongCDS)
      split_valid <- sample(valid, floor(length(valid) / 2))
      not_valid <- which(!strongCDS)
      if (length(not_valid) < 2) stop("No CDS were not valid, not allowed")
      split_not_valid <- sample(not_valid, floor(length(not_valid) / 2))
      indices <- sort(c(split_valid, split_not_valid))
      cds_total <- getCDSFiltered()

      cds_training <- cds_total[indices]
      saveRDS(cds_training, p(dataFolder, "/cds_training.rds"))
      aCDS <- cds_total[-indices]
      # Make artificial
      grl <- ORFik:::artificial.orfs(aCDS, end5 = max.artificial.length,
                                     start3 = -max.artificial.length)

      if (length(cds_total) != (length(cds_training) + length(grl)))
        stop("Wrong making of artificial CDS")
      save(grl, file = p(dataFolder, "/uniqueUorfsAsGRWithTx",".rdata"))
      strongCDS <- strongCDS[indices]
      pos <- pos[indices,]
    }
    message(p("Valid CDS for positive set: ", sum(strongCDS),
              " as ratio of total: ", round(sum(strongCDS) / length(strongCDS), 3)))
    message(p("Valid trailers for negative set: ", sum(weakThree),
              " as ratio of total: ", round(sum(weakThree) / length(weakThree), 3)))

    if (sum(strongCDS) < 100) warning("Training on less than 100 valid 'active' CDS,",
                                      "you might have low quality data!")
    if (sum(weakThree) < 100) warning("Training on less than 100 valid 'non-active' trailers",
                                      "verify your annotation has enough trailers!")

    if (length(strongCDS) != nrow(pos)) stop("Wrong making of strongCDS!")
    negR <- data.table(rbind(pos[!strongCDS, ], neg[weakThree, ])) # add bad cds to neg
    pos <- data.table(rbind(neg[!weakThree, ], pos[strongCDS, ]))
    neg <- negR
    training <- data.table(rbind(pos, neg)); colnames(training) <- features

    training <- uORFomePipe:::fixNAandINFTable(training)
    y <- as.factor(c(rep(1, nrow(pos)), rep(0, nrow(neg))))
    training <- data.table(y, training)

    ############# Make start codon best grouping
    dCDSThree <- uORFomePipe:::getBestIsoformStartCodonCoverage(tissue, cdsAndThree = TRUE,
                                                                mode = mode)
    pos <- dCDSThree[[1]]
    neg <- dCDSThree[[2]]
    if (length(strongCDS) != nrow(pos)) stop("wrong making of strongCDS for best codon!")
    if (length(weakThree) != nrow(neg)) stop("wrong making of weakThree for best codon!")
    negR <- data.table(rbind(pos[!strongCDS, ], neg[weakThree, ]))      # add bad cds to neg
    pos <- data.table(rbind(neg[!weakThree, ], pos[strongCDS, ]))  #!!! UPDATE if new features
    neg <- negR

    dCDSThree <- data.table(rbind(pos, neg))
    dCDSThree <- dCDSThree[readHits >= quantile(dCDSThree$readHits, 0.85),]
    dInts <- dCDSThree[dCDSThree[, .I[readHits == max(readHits)], by=group]$V1]
    ints <- c(length(strongCDS) + which(!weakThree),
              which(strongCDS), which(!strongCDS),
              length(strongCDS) + which(weakThree))
    if (length(ints) != nrow(training)) stop("wrong making of ints!")
    training$startCodonPerGroupBest <- ints %in% dInts$index
    #analysis
    message(p("Correlation between features before prediction for training set tissue: ", tissue))
    print(cor(data.matrix(training), use = "complete.obs"))

    saveRDS(training, file = paste0("features/TrainingData_",tissue,".rds"))
  }
  message("Training data on CDS and 3' UTRs finished")

  return(invisible(NULL))
}

#' Predict table for uORFs
#'
#' Get ribo-seq features for uORFs
#' Will only load if it already exists
#' @inheritParams makeTrainingData
#' @export
makeORFPredictionData <- function(tissues = "combined",
                                  mode = "uORF",
                                  features = c("countRFP", "disengagementScores", "entropyRFP", "floss",
                                               "fpkmRFP","ioScore", "ORFScores", "RRS", "RSS", "startCodonCoverage",
                                               "startRegionCoverage","startRegionRelative")) {
  if (is(tissues, "experiment")) tissues <- readTable("experiment_groups")[[1]]
  preName <- ifelse(mode == "CDS", "verify", "")
  if (length(tissues) == 1) {
    if(file.exists(paste0("features/PredictionData_",tissues, preName,".rds"))) {
      return(readRDS(paste0("features/PredictionData_",tissues, preName,".rds")))
    }
  }

  for (tissue in tissues) {
    if (file.exists(paste0("features/PredictionData_",tissue, preName,".rds"))) next
    pre_features <- features
    if (preName != "") {
      pre_features <- grep(pattern = p(preName, features, collapse = "|"),
                           x = listTables(), value = TRUE)
    }

    pos <-bplapply(pre_features, function(x) {
      uORFomePipe:::getTissueFromFeatureTable(tableName = x, tissue = tissue)
    }); pos <- setDT(unlist(pos, recursive = FALSE))
    colnames(pos) <- features
    predictors <- fixNAandINFTable(pos)
    # filter on isoforms
    d <- uORFomePipe:::getBestIsoformStartCodonCoverage(tissue, mode = mode)
    # combine filter with ribo-seq prediction
    d <- d[readHits > quantile(d$readHits, 0.973),]
    d <- d[d[, .I[readHits == max(readHits)], by=group]$V1]
    predictors$startCodonPerGroupBest <- seq.int(1, nrow(predictors)) %in% d$index
    message(p("Correlation between features before prediction for prediction set for tissue: ",
              tissue))
    print(cor(data.matrix(predictors), use = "complete.obs"))
    saveRDS(predictors, file = paste0("features/PredictionData_",tissue, preName,".rds"))
  }
  return(invisible(NULL))
}

#' Get all sequence features for ORFs
#'
#' Extracted from the data base
#' @return data.table with sequence features
#' @export
getAllSequenceFeaturesTable <- function() {
  return(fread(file = "features/uORFSequenceFeatures.csv", header = TRUE))
}

#' A filter per stop codon group
#'
#' Uses the ribo-seq libraries for ribo-seq validation
#' get start for each in group, count overlaps, return orf with
#' highest per group
#' @param cdsAndThree logical, default: FALSE, that is run of uORFs, if TRUE
#' run for cds and trailers.
#' @return list of length 1 or 2, (2 if cdsAndThree is TRUE),
#' list contains numeric vectors of best start codons by coverage.
getBestIsoformStartCodonCoverage <- function(tissue, cdsAndThree = FALSE,
                                             mode) {
  # reduce isoform overlaps by highest start codon reads per group
  if (cdsAndThree) {
    if(file.exists(paste0("features/bestStartCodonsCDSTHREE_", tissue, ".rdata"))) {
      load(paste0("features/bestStartCodonsCDSTHREE_", tissue, ".rdata"))
      return(dCDSThree)
    }
    # cds
    match <- names(getCDSFiltered()) %in% names(getCDSTraining(mode = mode))
    sg <- stopCodons(getCDSTraining(mode = mode), is.sorted = TRUE)
    dt <- data.table(stopCodonGrouping = uniqueOrder(sg))
    #insertTable(dt, "cdsstopCodonGrouping")
    grps <- dt$stopCodonGrouping
    counts <- rowMeans(readTable("cdsstartCodonCoverage"))[match]
    dCDS <- data.table(readHits = counts, group = grps, index = seq.int(length(grps)))

    # trailers
    sg <- stopCodons(threeUTRs[widthPerGroup(threeUTRs) > 5], is.sorted = TRUE)
    dt <- data.table(stopCodonGrouping = uniqueOrder(sg))
    #insertTable(dt, "threestopCodonGrouping")
    grps <- dt$stopCodonGrouping
    counts <- rowMeans(readTable("threestartCodonCoverage"))
    dThree <- data.table(readHits = counts, group = grps, index = seq.int(length(grps)))
    dThree$group <- dThree$group + max(dCDS$group)
    dThree$index <- dThree$index + max(dCDS$index)
    dCDSThree <- list(dCDS, dThree)
    save(dCDSThree, file = paste0("features/bestStartCodonsCDSTHREE_", tissue, ".rdata"))
    return(dCDSThree)
  }
  addit <- ifelse(mode == "CDS", "verify", "")
  additional <- p("_", addit)
  if(file.exists(p("features/bestStartCodons_", tissue, additional,".rdata"))) {
    load(p("features/bestStartCodons_", tissue, additional,".rdata"))
    return(d)
  }

  grps <- readTable(p(addit, "stopCodonGrouping"))$stopCodonGrouping
  counts <- rowMeans(readTable(p(addit, "startCodonCoverage"), with.IDs = F))
  d <- data.table(readHits = counts, group = grps, index = seq.int(length(grps)))
  save(d, file = p("features/bestStartCodons_", tissue, additional,".rdata"))
  return(d)
}
