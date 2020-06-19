#' Get CDS and 3'UTR TrainingData of ribo seq features
#'
#' Positive set is cds, negative is downstream region of 3' UTRs
#' @param tissue Tissue to train on, use "combined" if you want all in one,
#' first run of training it is a ORFik experiment.
#' @param features features to train model on, must exist in database
#' @return invisible(NULL), saved to disc
makeTrainingData <- function(tissues = "combined",
                             features = c("countRFP", "disengagementScores", "entropyRFP", "startRegionCoverage","startRegionRelative",
                                          "floss", "ioScore", "ORFScores","fpkmRFP", "RRS", "RSS", "startCodonCoverage")) {
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


    pos <-bplapply(posFeatureNames, function(x, tissue) {
      uORFomePipe:::getTissueFromFeatureTable(x, tissue = tissue)
    }, tissue = tissue); pos <- as.matrix(setDT(unlist(pos, recursive = FALSE)))

    neg <-bplapply(negFeatureNames, function(x) {
      uORFomePipe:::getTissueFromFeatureTable(tableName = x, tissue = tissue)
    }); neg <- as.matrix(setDT(unlist(neg, recursive = FALSE)))

    # Filter out translating 3' UTRs
    weakThree <- (((neg[,1] < quantile(neg[,1], 0.95)) |
                       (neg[,12] < quantile(neg[,12], 0.98))) | neg[,5] < 1.1)
    # Filter out non-translating CDS'
    strongCDS <- uORFomePipe:::strongCDS(pos[,1], pos[,9], pos[,12], pos[,5])

    negR <- data.table(rbind(pos[!strongCDS, ], neg[weakThree, ])) # add bad cds to neg
    pos <- data.table(rbind(neg[!weakThree, ], pos[strongCDS, ]))
    neg <- negR

    training <- data.table(rbind(pos, neg)); colnames(training) <- features


    training <- uORFomePipe:::fixNAandINFTable(training)
    y <- as.factor(c(rep(1, nrow(pos)), rep(0, nrow(neg))))
    training <- data.table(y, training)

    dCDSThree <- uORFomePipe:::getBestIsoformStartCodonCoverage(tissue, cdsAndThree = TRUE)
    pos <- dCDSThree[[1]]
    neg <- dCDSThree[[2]]
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
#' @param tissue Tissue to train on, use "combined" if you want all in one
#' @param features features to train model on, must exist in database
#' @export
makeORFPredictionData <- function(tissues = "combined",
                                  features = c("countRFP", "disengagementScores", "entropyRFP", "startRegionCoverage","startRegionRelative",
                                               "floss", "ioScore", "ORFScores","fpkmRFP", "RRS", "RSS", "startCodonCoverage")) {
  if (is(tissues, "experiment")) tissues <- unique(tissues$stage)
  if (length(tissues) == 1) {
    if(file.exists(paste0("features/PredictionData_",tissues,".rds"))) {
      return(readRDS(paste0("features/PredictionData_",tissues,".rds")))
    }
  }

  for (tissue in tissues) {
    if (file.exists(paste0("features/PredictionData_",tissue,".rds"))) next
    pos <-bplapply(features, function(x) {
      uORFomePipe:::getTissueFromFeatureTable(tableName = x, tissue = tissue)
    }); pos <- setDT(unlist(pos, recursive = FALSE))
    colnames(pos) <- features
    predictors <- fixNAandINFTable(pos)
    # filter on isoforms
    d <- getBestIsoformStartCodonCoverage(tissue)
    # combine filter with ribo-seq prediction
    d <- d[readHits > quantile(d$readHits, 0.973),]
    d <- d[d[, .I[readHits == max(readHits)], by=group]$V1]
    predictors$startCodonPerGroupBest <- seq.int(1, nrow(predictors)) %in% d$index
    message(p("Correlation between features before prediction for prediction set for tissue: ",
              tissue))
    print(cor(data.matrix(predictors), use = "complete.obs"))
    saveRDS(predictors, file = paste0("features/PredictionData_",tissue,".rds"))
  }
  return(invisible(NULL))
}

#' Get all sequence features for ORFs
#'
#' Extracted from the data base
#' @export
getAllSequenceFeaturesTable <- function() {
  return(fread(file = "features/uORFSequenceFeatures.csv", header = TRUE))
}

#' A filter per stop codon group
#'
#' Uses the ribo-seq libraries for ribo-seq validation
#' get start for each in group, count overlaps, return orf with
#' highest per group
getBestIsoformStartCodonCoverage <- function(tissue, cdsAndThree = FALSE) {
  # reduce isoform overlaps by highest start codon reads per group
  if (cdsAndThree) {
    if(file.exists(paste0("features/bestStartCodonsCDSTHREE_", tissue, ".rdata"))) {
      load(paste0("features/bestStartCodonsCDSTHREE_", tissue, ".rdata"))
      return(dCDSThree)
    }

    grps <- readTable("cdsstopCodonGrouping")$stopCodonGrouping
    counts <- rowMeans(readTable("cdsstartCodonCoverage"))
    dCDS <- data.table(readHits = counts, group = grps, index = seq.int(length(grps)))

    grps <- readTable("threestopCodonGrouping")$stopCodonGrouping
    counts <- rowMeans(readTable("threestartCodonCoverage"))
    dThree <- data.table(readHits = counts, group = grps, index = seq.int(length(grps)))
    dThree$group <- dThree$group + max(dCDS$group)
    dThree$index <- dThree$index + max(dCDS$index)
    dCDSThree <- list(dCDS, dThree)
    save(dCDSThree, file = paste0("features/bestStartCodonsCDSTHREE_", tissue, ".rdata"))
    return(dCDSThree)
  }
  if(file.exists(paste0("features/bestStartCodons_", tissue, ".rdata"))) {
    load(paste0("features/bestStartCodons_", tissue, ".rdata"))
    return(d)
  }

  grps <- readTable("stopCodonGrouping")$stopCodonGrouping
  counts <- rowMeans(readTable("startCodonCoverage", with.IDs = F))
  d <- data.table(readHits = counts, group = grps, index = seq.int(length(grps)))
  save(d, file = paste0("features/bestStartCodons_", tissue, ".rdata"))
  return(d)
}
