#' Get start codon bias usage
#'
#' Usually something like this
#' good: ATG, CTG, ACG and GTG
#' bad: AAG AND AGG
#' @param hits a logical vector of TRUE, FALSE (TRUE is predicted)
#' @param tissue default NULL, or a character name of tissue to print
#' @param StartCodons character, sequences of start codons, default:
#' getAllSequenceFeaturesTable()$StartCodons
#' @export
startCodonMetrics <- function(hits, tissue = NULL,
                              StartCodons = getAllSequenceFeaturesTable()$StartCodons) {
  if (!is.null(tissue)) message(p("Start codon distribution for tissue: ", tissue))
  if (sum(hits) == 0) {
    warning("No ORFs predicted as active, check your input data!")
    return(invisible(NULL))
  }
  message("0 is downregulated, 1 is upregulated by chi.sq test")
  ySeq <- rep(0, length(hits))
  ySeq[hits] <- 1
  StartResultsSequences <- chisq.test(table(data.frame(StartCodons, prediction = as.factor(ySeq))))
  res <- round(StartResultsSequences$residuals,1)
  res <- res[order(res[,2],decreasing = T),]
  count <- table(StartCodons[hits])[rownames(res)]
  relativeCount <- round(table(StartCodons[hits])/sum(hits), 2)[rownames(res)]
  res <- cbind(res, relative = round(res[,2]/max(res[,2]), 2), count, relativeCount)
  print(paste("number of uORFs predicted translated:", sum(hits)))
  print(res)
  return(invisible(NULL))
}

#' General feature analysis
#' @param prediction a data.table, 3 columns: (prediction, p0 and p1).
#' Made from uORFomePipe
#' @export
featureAnalysis <- function(prediction, tissue){
  uorfTable <- makeORFPredictionData(tissue)
  uts <- uorfTable

  uorfData <- getAllSequenceFeaturesTable()
  uds <- data.table(uorfData, stringsAsFactors = TRUE)

  pred <- prediction$predict
  print("Output feature summaries")
  print("For ribo seq prediction")
  for(i in colnames(uorfTable)) {
    message(i)
    print(data.table(pos = data.table(summary(uts[pred == 1, i, with = F]))[,3],
                     neg = data.table(summary(uts[pred == 0, i, with = F]))[,3]))
  }
  print("On sequence feature table")
  for(i in colnames(uorfData)) {
    message(i)
    print(data.table(pos = data.table(summary(uds[pred == 1, i, with = F]))[,3],
                     neg = data.table(summary(uds[pred == 0, i, with = F]))[,3]))
  }

  message("Ribo-seq features (mean value) seperated by start codon")
  dt <- data.table(); uniqueStarts <- unique(uorfData$StartCodons)
  for (i in uniqueStarts) {
    dt <- cbind(dt, i = round(colMeans(uorfTable[uorfData$StartCodons == i], na.rm = TRUE), 3))
  }
  dt <- t(dt); colnames(dt) <- colnames(uorfTable)
  rownames(dt) <- uniqueStarts
  print(dt)

  message("Start codon distribution at > 75% certain predication")
  print(uORFomePipe:::startCodonMetrics(prediction$p1 > 0.75))

  # more checks on orf score
  message("Feature mean values at different predication thresholds, 0.95 is strong positive, 0.05 is strong negative")
  x <- seq(0,1, 0.05)[2:20]
  dt <- data.table()
  for(i in x) {
    dt <- cbind(dt , round(colMeans(uorfTable[prediction$p1 > i]), 3))
  }
  dt <- t(dt); colnames(dt) <- colnames(uorfTable)
  rownames(dt) <- x
  print(dt)
  return(invisible(NULL))
}

checkBadPred <- function(){
  # region check
  uorfData <- getAllSequenceFeaturesTable()
  agg <- uorfData$StartCodons == "AGG" & prediction$predict == 1

  starts <- startSites(grl[agg], T, T, T)
  startRegion <- windowPerGroup(starts, tx, 6, 9)

  seqs <- getSequencesFromFasta(startRegion, T)
  allOver <- c("CTG|ATG|ACG|GTG|TTG")
  notBestStart <- (uorfTable$startCodonPerGroupBest == F)[agg]
  hits <- grep(x = seqs, pattern = allOver)
  valid <- (!(seq.int(1,sum(agg)) %in% hits))
  index <- which((!notBestStart & valid))
  t <- 81
  for(i in index[t:(t+2)]){
    print(i)
    print(prediction$predict[agg][i])
    print(uorfData[agg,][i])
    print(uorfTable[agg,][i])
    print(grl[agg][i])
  }
  toKeep <- which(agg)[index]
}

#' Classification bounaries
#'
#' @param tissues the tissues wanted to check
findClassificationBoundary <- function(tissues){
  x <- seq(0, 1, 0.1)[2:10]
  # for riboseq prediction
  for(tissue in tissues) {
    load(paste0("prediction_model/prediction_",tissue, ".rdata"))
    hits <- unlist(lapply(x, function(y) sum(as.logical(prediction$p1 > y))),
                   use.names = F)
    plot(x, hits, main = tissue)
  }
  # we pick 0.5 from here
  # for sequence prediction
  for(tissue in tissues) {
    load(paste0("prediction_model/finalPrediction_",tissue, ".rdata"))
    hits <- unlist(lapply(x, function(y) sum(as.logical(uorfPrediction$p1 > y))),
                   use.names = F)
    plot(x, hits, main = tissue,
         xlab = "prediction cutoff", ylab = "# of predicted uORFs")
  }
  # we pick 0.75 from here

  # validate boundaries
  for(j in seq(0.5, 1, 0.05)[2:10]) {
    print(paste("pred:", j))
    for(i in 1:5) {
      d <- getBestIsoformStartCodonCoverage()[readHits >= i & prediction$predict == 1 & prediction$p1 > j,]
      d <- d[, .SD[which.max(readHits)], by = group]
      print(i)
      print(round((table(uorfData$StartCodons[d$index])[6]/table(uorfData$StartCodons[d$index])[1]), 2))
    }
  }
  # validate read count

  start <- 0.5
  stop <- 1.5

  for(i in start:stop) {
    d <- getBestIsoformStartCodonCoverage()[readHits >= i & prediction$p1 >= 0.65,]
    d <- d[, .SD[which.max(readHits)], by = group]
    print(i)
    len[(i-start+1)] <- length(d$index)
  }
  plot(start:stop, len, main = paste("read count # change in", tissue),
       xlab = "# of reads as cutoff", ylab = "# of predicted uORF groups")

  x <- seq(0, 1, 0.05)[2:20]
  for(i in x){
    print(paste("x:   ",i))
    hits <- which(as.logical(uorfPrediction$p1 > i))
    # good: ATG, CTG, procaryote: GTG, TTG
    # bad: AAG AND AGG
    ySeq <- rep(0, nrow(uorfPrediction))
    ySeq[hits] <- 1
    StartResultsSequences <- chisq.test(table(data.frame(uorfData$StartCodons, prediction = as.factor(ySeq))))
    print(paste("number of uORFs predicted translated:", length(hits)))
    print(round(StartResultsSequences$residuals,1))
    print(round(table(uorfData$StartCodons[hits])/length(hits), 2))
    print("\n\n")
  }
  # We choose 0.75 here
}
