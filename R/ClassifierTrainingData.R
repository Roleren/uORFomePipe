#' Get CDS and 3'UTR TrainingData of ribo seq features
#'
#' Positive set is cds, negative is downstream region of 3' UTRs
#' @param tissue Tissue to train on, use "all" if you want all in one
makeTrainingData <- function(tissue,
                             features = c("countRFP", "disengagementScores", "entropyRFP", "startRegionCoverage","startRegionRelative",
                                          "floss", "ioScore", "ORFScores","fpkmRFP", "RRS", "RSS", "startCodonCoverage")) {

  if (file.exists(paste0("features/TrainingData_",tissue,".rds"))) {
    return(readRDS(paste0("features/TrainingData_",tissue,".rds")))
  }
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
  filterThree <- (((neg[,1] < quantile(neg[,1], 0.95)) |
                     (neg[,12] < quantile(neg[,12], 0.98))) | neg[,5] < 1.1)
  # Filter out non-translating CDS'
  filterCDS <- uORFomePipe:::strongCDS(pos[,1], pos[,9], pos[,12], pos[,5])

  negR <- data.table(rbind(pos[!filterCDS, ], neg[filterThree, ])) # add bad cds to neg
  pos <- data.table(rbind(neg[!filterThree, ], pos[filterCDS, ]))
  neg <- negR

  training <- data.table(rbind(pos, neg))
  colnames(training) <- features

  training <- uORFomePipe:::fixNAandINFTable(training)
  y <- as.factor(c(rep(1, nrow(pos)), rep(0, nrow(neg))))
  training <- data.table(y, training)

  dCDSThree <- uORFomePipe:::getBestIsoformStartCodonCoverage(cdsAndThree = TRUE)
  pos <- dCDSThree[[1]]
  neg <- dCDSThree[[2]]
  negR <- data.table(rbind(pos[!filterCDS, ], neg[filterThree, ]))      # add bad cds to neg
  pos <- data.table(rbind(neg[!filterThree, ], pos[filterCDS, ]))  #!!! UPDATE if new features
  neg <- negR

  dCDSThree <- data.table(rbind(pos, neg))
  dCDSThree <- dCDSThree[readHits >= quantile(dCDSThree$readHits, 0.85),]
  dInts <- dCDSThree[dCDSThree[, .I[readHits == max(readHits)], by=group]$V1]
  ints <- c(length(filterCDS) + which(!filterThree),
            which(filterCDS),
            which(!filterCDS),
            length(filterCDS) + which(filterThree))
  if (length(ints) != nrow(training)) stop("wrong making of ints!")
  training$startCodonPerGroupBest <- ints %in% dInts$index
  #analysis
  message("Correlation between features before prediction for training set")
  print(cor(data.matrix(training), use = "complete.obs"))

  saveRDS(training, file = paste0("features/TrainingData_",tissue,".rds"))
  return(training)
}

#' Predict table for uORFs
#'
#' Get ribo-seq features for uORFs
#' Will only load if it already exists
#' @param tissue default "all"
#' @export
makeORFPredictionData <- function(tissue = "all",
                                  features = c("countRFP", "disengagementScores", "entropyRFP", "startRegionCoverage","startRegionRelative",
                                               "floss", "ioScore", "ORFScores","fpkmRFP", "RRS", "RSS", "startCodonCoverage")) {
  if(file.exists(paste0("features/PredictionData_",tissue,".rds"))) {
    return(readRDS(paste0("features/PredictionData_",tissue,".rds")))
  }

  pos <-bplapply(features, function(x) {
    uORFomePipe:::getTissueFromFeatureTable(tableName = x, tissue = tissue)
  }); pos <- setDT(unlist(pos, recursive = FALSE))
  colnames(pos) <- features
  predictors <- uORFomePipe:::fixNAandINFTable(pos)
  # filter on isoforms
  d <- uORFomePipe:::getBestIsoformStartCodonCoverage()
  # combine filter with ribo-seq prediction
  d <- d[readHits > quantile(d$readHits, 0.973),]
  d <- d[d[, .I[readHits == max(readHits)], by=group]$V1]
  predictors$startCodonPerGroupBest <- seq.int(1, nrow(predictors)) %in% d$index
  message("Correlation between features before prediction for training set")
  print(cor(data.matrix(predictors), use = "complete.obs"))
  saveRDS(predictors, file = paste0("features/PredictionData_",tissue,".rds"))
  return(predictors)
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
getBestIsoformStartCodonCoverage <- function(cdsAndThree = F) {
  # reduce isoform overlaps by highest start codon reads per group
  if (cdsAndThree) {
    if(file.exists(paste0("features/bestStartCodonsCDSTHREE.rdata"))) {
      load(paste0("features/bestStartCodonsCDSTHREE.rdata"))
      return(dCDSThree)
    }

    getCDS()
    cds <- cds[widthPerGroup(cds) > 5]
    sg <- stopCodons(cds, is.sorted = TRUE)
    uo <- uniqueOrder(sg) # <- grouping
    counts <- rowMeans(readTable("cdsstartCodonCoverage"))
    dCDS <- data.table(readHits = counts, group = uo, index = seq.int(length(uo)))

    getThreeUTRs()
    threeUTRs <- threeUTRs[widthPerGroup(threeUTRs) > 5]
    sg <- stopCodons(threeUTRs, is.sorted = T)
    uo <- uniqueOrder(sg) # <- grouping
    counts <- rowMeans(readTable("threestartCodonCoverage"))
    dThree <- data.table(readHits = counts, group = uo, index = seq.int(length(uo)))
    dThree$group <- dThree$group + max(dCDS$group)
    dThree$index <- dThree$index + max(dCDS$index)
    dCDSThree <- list(dCDS, dThree)
    save(dCDSThree, file = paste0("features/bestStartCodonsCDSTHREE.rdata"))
    return(dCDSThree)
  }
  if(file.exists(paste0("features/bestStartCodons.rdata"))) {
    load(paste0("features/bestStartCodons.rdata"))
    return(d)
  }
  uo <- readTable("stopCodonGrouping")$stopCodonGrouping

  counts <- rowMeans(readTable("startCodonCoverage", with.IDs = F))
  d <- data.table(readHits = counts, group = uo, index = seq.int(length(uo)))
  save(d, file = paste0("features/bestStartCodons.rdata"))
  return(d)
}