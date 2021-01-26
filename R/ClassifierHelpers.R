#' Train h2o rf model.
#' negDT if you want own samples for that
#' @param dt data.table of features to train on
#' @param cv number of cross validations, default 10
#' @inheritParams predictUorfs
#' @import h2o
forest <- function(dt, cv = 10, ntrees = 64, ip,
                   port, nthreads_h2o,  max_mem_size){
  h2o.init(ip = ip, port = port, nthreads = nthreads_h2o, max_mem_size = max_mem_size)
  # new h2o model training
  indices <- sample(1:nrow(dt), 0.6*nrow(dt), replace = F)
  validationIndices <- seq.int(nrow(dt))[-indices]
  trainingFrame <-  as.h2o(dt[indices, ])
  validationFrame <- as.h2o(dt[validationIndices,])
  forestH2o <- h2o.randomForest(y = "y", training_frame = trainingFrame,
                                validation_frame = validationFrame,
                                nfolds = cv, ntrees = ntrees)
  print(data.table(feature = forestH2o@model$variable_importances[,1],
                   round(forestH2o@model$variable_importances[,3:4], 2)))
  print(data.frame(name = rownames(forestH2o@model$cross_validation_metrics_summary)[c(1,17,19)],
                   score = round(as.double(forestH2o@model$cross_validation_metrics_summary[c(1,17,19),1]),
                                 2)))
  return(forestH2o)
}

#' Combine classifier and CAGE data, for final prediction table
#' @param tissue tissue to make combined prediction
#' @return data.table with 3 columns. 1. "prediction" (1 or 0), 2. p0 (probability for FALSE),
#' p1 (probability for TRUE)
makeCombinedPrediction <- function(tissues, dataFolder = get("dataFolder", envir = .GlobalEnv),
                                   mode,
                                   grl = getUorfsInDb(mode = mode)) {
  preName <- ifelse(mode == "CDS", "verify_", "")
  if (!tableNotExists(p(preName, "finalPredWithProb"))) {
    return(readTable(p(preName, "finalPredWithProb")))
  }
  # load data
  predicted_ribo <- data.table()
  for (tissue in tissues) {
    prediction <- readRDS(paste0("prediction_model/prediction_", preName, tissue, ".rds"))
    predicted_ribo <- cbind(predicted_ribo, prediction$predict == 1) # set value
  }
  colnames(predicted_ribo) <- tissues
  if (file.exists(paste0(dataFolder,"/tissueAtlas.rds"))) {
    tissueAtlas <- readRDS(paste0(dataFolder,"/tissueAtlas.rds"))[,-1]
    tissueAtlas$total <- rowSums(tissueAtlas) > 0
  } else tissueAtlas <- TRUE
  predicted_ribo$total <- rowSums(predicted_ribo) > 0

  cageTissuesPrediction <- tissueAtlas & predicted_ribo

  insertTable(cageTissuesPrediction, p(preName, "tissueAtlasByCageAndPred"), rmOld = TRUE)
  finalCagePred <- rowSums(cageTissuesPrediction) > 0
  insertTable(finalCagePred, p(preName, "finalCAGEuORFPrediction"), rmOld = TRUE)

  finalCagePred.dt <- data.table(prediction = finalCagePred, prediction[, 2:3])
  insertTable(finalCagePred.dt, p(preName, "finalPredWithProb"), rmOld = TRUE)

  startCodonMetrics(finalCagePred)
  export.bed12(grl, file = p("candidate_and_predicted_uORFs_", preName, "total", ".bed"), rgb = 255*finalCagePred)
  message("Prediction finished, now do the analysis you want")
  return(finalCagePred.dt)
}

#' Find active CDS
#'
#' A filter to say if a CDS has strong indication of translation,
#' used as part of positives training set
#' @param coverage (counts over region)
#' @param fpkm FPKM value of region
#' @param startCodonCoverage (Coverage of startcodon +/- 1 base region)
#' @param startRegionRelative Relative coverage of start region to upstream short region
#' @param ORFScore Periodicity score of triplets, > 0 if frame 0 has most reads.
#' @return logical TRUE/FALSE per row
strongCDS <- function(coverage, fpkm, startCodonCoverage, startRegionRelative, ORFScore) {
  filter <- (coverage > min(quantile(coverage, 0.25), 10)) & (fpkm > max(1, quantile(fpkm, 0.15))) &
    (startCodonCoverage > quantile(startCodonCoverage, 0.75)) &
    (startRegionRelative > 0.1) & (ORFScore > 1.0)
  return(filter)
}

#' Find non-active trailers
#'
#' A filter to say if a trailer has strong indication of no translation,
#' used as part of negative training set
#' @param coverage (counts over region)
#' @param startRegionRelative Relative coverage of start region to upstream short region
#' @param fpkm fragments per kilobase transcript per million
#' @param ORFScore Periodicity score of triplets, > 0 if frame 0 has most reads.
#' @return logical TRUE/FALSE per row
weakTrailer <- function(coverage, fpkm, startRegionRelative,
                         ORFScore) {
  filter <- (((coverage < quantile(coverage, 0.95)) |
   (startRegionRelative < max(1, quantile(startRegionRelative, 0.97)))) |
     fpkm < 1.1 | ORFScore < 0.5)
  return(filter)
}
#' Get specific feature for specific tissues
#'
#' Grouped by rowMeans
#' @param tableName name of table in sql database
#' @param tissue character, single tissue
#' @importFrom data.table data.table
#' @return data.table of row mean values by tissue
getTissueFromFeatureTable <- function(tableName, tissue) {
  if (tableNotExists(tableName)) stop(paste("table does not exist:",
                                            tableName))
  if (length(tissue) > 1) stop("Length of tissue must be exactly 1!")
  uORFomePipe:::createDataBase(p(dataBaseFolder, "/uorfCatalogue.sqlite"))
  riboTable <- readTable(tableName)
  tissues <- readTable("experiment_groups")[[1]]
  if (is.null(tissue) | (tissue == "combined")) {
  } else if (tissue %in% tissues){
    indices <- tissues == tissue
    riboTable <- riboTable[, indices, with = FALSE]
  } else stop("tissue does not exist in db")
  return(data.table(rowMeans(riboTable)))
}

#' Get all data for combined dataset
#' @param prediction Which prediction to use, default:
#' prediction = readTable("finalPredWithProb"
#' @param stringsAsFactors TRUE
#' @return data.table of prediction, NGS features and sequence features
getAllUORFData <- function(prediction = readTable("finalPredWithProb"),
                           stringsAsFactors = TRUE) {
  uorfTable <- makeORFPredictionData(tissue)
  uorfData <- getAllSequenceFeaturesTable()

  return(data.table(prediction, uorfTable, uorfData, stringsAsFactors = stringsAsFactors))
}

#' #' Not used in prediction
#' #'
#' #' A test I made to check what good features are
#' filterHardOnes <- function(prediction, tissue = "all") {
#'   prediction$filtered <- rep(F, nrow(prediction))
#'   uorfTable <- makeUORFPredicateTable()
#'   uorfData <- getAllSequenceFeaturesTable()
#'   grl <- getUorfsInDb()
#'   getCDS()
#'   getCageTx()
#'   getFasta()
#'   table <- startCodonMetrics(as.logical(prediction[,3] >= 0.50))
#'   badStarts <- table[,1] > 0 & table[,2] < 0
#'   badStartCodons <- rownames(table)[badStarts]
#'   goodStartCodons <- paste(rownames(table)[!badStarts], collapse = "|")
#'   goodUORFs <- grl[uorfData$StartCodons %in%  rownames(table)[!badStarts]]
#'   res = c()
#'   for(codon in badStartCodons) {
#'     # find region
#'     agg <- uorfData$StartCodons == codon & prediction[,3] >= 0.75
#'     starts <- startSites(grl[agg], T, T, T)
#'     # make string
#'     hits <- grep(x = ORFik:::startRegionString(grl[agg], tx, fa, 6, 9), pattern = goodStartCodons)
#'     hitsUp <- grep(x = ORFik:::startRegionString(grl[agg], tx, fa, 5, 0), pattern = stopDefinition(1))
#'     hitsOverlapsBetter <- starts %over% goodUORFs
#'     hitsCDS <- to(findOverlaps(startSites(cds, T, T, T), starts, maxgap = 3))
#'     valid <- (!(seq.int(1,sum(agg)) %in% c(hits, hitsUp, hitsOverlapsBetter, hitsCDS)))
#'
#'     # find indices
#'     notBestStart <- (uorfTable$startCodonPerGroupBest == F)[agg]
#'     index <- which((!notBestStart & valid))
#'     toKeep <- which(agg)[index]
#'     res <- c(res, toKeep)
#'   }
#'   if(length(res) != length(unique(res))) stop("error in res creation!")
#'   hits <- (grep(x = uorfData$StartCodons, pattern = paste(badStartCodons, collapse = "|")))
#'
#'   prediction$predict[hits] <- 0
#'   prediction$p0[hits] <- 1
#'   prediction$p1[hits] <- 0
#'   prediction$filtered[hits] <- T
#'
#'   prediction$predict[res] <- 1
#'   prediction$p0[res] <- 0
#'   prediction$p1[res] <- 1
#'   startCodonMetrics(as.logical(prediction[,3] >= 0.50))
#'   save(prediction, file = paste0("forests/finalPrediction_filtered",tissue, ".rdata"))
#'   return(prediction)
#' }
