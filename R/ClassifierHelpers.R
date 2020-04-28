#' Train h2o rf model.
#' negDT if you want own samples for that
forest <- function(dt, cv = 10, ntrees = 64, nthreads = 40,  max_mem_size = "200G"){
  library(h2o)
  h2o.init(nthreads = nthreads, max_mem_size = max_mem_size, port = 20050)
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
makeCombinedPrediction <- function(tissue, dataFolder = get("dataFolder", envir = .GlobalEnv)) {
  # load data
  prediction <- readRDS(paste0("prediction_model/prediction_",tissue,".rds"))
  load(paste0(dataFolder,"/tissueAtlas.rdata"))

  some <- prediction$predict == 1 # set value

  cageTissuesPrediction <- copy(tissueAtlas)
  for(i in colnames(cageTissuesPrediction)[-1]) {
    cageTissuesPrediction[, paste(i) := (tissueAtlas[,i, with=F] & some)]
  }
  insertTable(cageTissuesPrediction, "tissueAtlasByCageAndPred", rmOld = T)
  finalCagePred <- rowSums(cageTissuesPrediction[,-1]) > 0
  insertTable(finalCagePred, "finalCAGEuORFPrediction", rmOld = T)

  startCodonMetrics(finalCagePred)
  return(data.table(prediction = finalCagePred, prediction[, 2:3]))
}

#' strongCDS
#'
#' A filter to say if a CDS has strong indication of translation, used as strong positives
strongCDS <- function(coverage, fpkm, startCodonCoverage, fiveRegionRelative) {
  filter <- (coverage > min(quantile(coverage, 0.25), 10)) & (fpkm > quantile(fpkm, 0.15)) &
    (startCodonCoverage > quantile(startCodonCoverage, 0.75)) & fiveRegionRelative > 0.95
  return(filter)
}

# Get only specific tissues, or all.
# Grouped by rowMeans
getTissueFromFeatureTable <- function(tableName, tissue) {
  if (tableNotExists(tableName)) stop(paste("table does not exist:",
                                            tableName))
  uORFomePipe:::createDataBase(p(dataBaseFolder, "/uorfCatalogue.sqlite"))
  riboTable <- readTable(tableName, with.IDs = FALSE)
  if (is.null(tissue) | (tissue == "all")) {
    print("Grouping all together")
  } else if ((tissue %in% as.character(unique(rpfSamples$tissue)))){
    indices <- rpfSamples$tissue == tissue
    riboTable <- riboTable[,indices, with = F]
  } else stop("tissue does not exist in db")
  return(data.table(rowMeans(riboTable)))
}

#' Not used in prediction
#'
#' A test I made to check what good features are
filterHardOnes <- function(prediction, tissue = "all") {
  prediction$filtered <- rep(F, nrow(prediction))
  uorfTable <- makeUORFPredicateTable()
  uorfData <- getAllSequenceFeaturesTable()
  grl <- getUorfsInDb()
  getCDS()
  getCageTx()
  getFasta()
  table <- startCodonMetrics(as.logical(prediction[,3] >= 0.50))
  badStarts <- table[,1] > 0 & table[,2] < 0
  badStartCodons <- rownames(table)[badStarts]
  goodStartCodons <- paste(rownames(table)[!badStarts], collapse = "|")
  goodUORFs <- grl[uorfData$StartCodons %in%  rownames(table)[!badStarts]]
  res = c()
  for(codon in badStartCodons) {
    # find region
    agg <- uorfData$StartCodons == codon & prediction[,3] >= 0.75
    starts <- startSites(grl[agg], T, T, T)
    # make string
    hits <- grep(x = ORFik:::startRegionString(grl[agg], tx, fa, 6, 9), pattern = goodStartCodons)
    hitsUp <- grep(x = ORFik:::startRegionString(grl[agg], tx, fa, 5, 0), pattern = stopDefinition(1))
    hitsOverlapsBetter <- starts %over% goodUORFs
    hitsCDS <- to(findOverlaps(startSites(cds, T, T, T), starts, maxgap = 3))
    valid <- (!(seq.int(1,sum(agg)) %in% c(hits, hitsUp, hitsOverlapsBetter, hitsCDS)))

    # find indices
    notBestStart <- (uorfTable$startCodonPerGroupBest == F)[agg]
    index <- which((!notBestStart & valid))
    toKeep <- which(agg)[index]
    res <- c(res, toKeep)
  }
  if(length(res) != length(unique(res))) stop("error in res creation!")
  hits <- (grep(x = uorfData$StartCodons, pattern = paste(badStartCodons, collapse = "|")))

  prediction$predict[hits] <- 0
  prediction$p0[hits] <- 1
  prediction$p1[hits] <- 0
  prediction$filtered[hits] <- T

  prediction$predict[res] <- 1
  prediction$p0[res] <- 0
  prediction$p1[res] <- 1
  startCodonMetrics(as.logical(prediction[,3] >= 0.50))
  save(prediction, file = paste0("forests/finalPrediction_filtered",tissue, ".rdata"))
  return(prediction)
}