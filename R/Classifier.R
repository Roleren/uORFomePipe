#' Predict uORFs using random forrest classification
#'
#' Training data
#' pos set is by CDS (coding sequences)
#' Neg seg is by trailers (3'UTRs)
#'
#' @param tissues Groups to train on, default (readTable("experiment_groups")[[1]]),
#'  use "combined" if you want mean of all groups
#' @param ip h2o cluster ip, default: "localhost".
#' @param port h2o cluster port, default: 54321
#' @param nthreads_h2o number of cores for H20. Default:
#' \code{max(min(45, as.integer(BiocParallel::bpparam()$workers/2)), 1)}
#' @param max_mem_size max allowed memory for H20: default ("200G")
#' @param mode character, default "uORF"
#' @export
#' @import h2o
predictUorfs <- function(tissues = readTable("experiment_groups")[[1]],
                         ip = "localhost",
                         port = 54321,
                         nthreads_h2o = max(min(45, as.integer(BiocParallel::bpparam()$workers/2)), 1),
                         max_mem_size = "200G",  mode = "uORF") {
  h2o.started <- FALSE
  preName <- ifelse(mode == "CDS", "verify_", "")
  nTissues <- length(tissues)
  for (tissue in tissues) {
    message(p("(",which(tissues == tissue), "/", nTissues,") Prediction for tissue: ", tissue))

    # make uORFTable
    if(file.exists(p("prediction_model/prediction_", preName, tissue, ".rds"))) next
    h2o.started <- TRUE
    training_model <- trainClassifier(tissue, ip, port, nthreads_h2o, max_mem_size)
    # Use training model from cds and trailers to predict on ORF prediction data (uORFs)
    prediction <- as.data.table(h2o.predict(training_model, as.h2o(makeORFPredictionData(tissue, mode))))
    uORFomePipe:::startCodonMetrics(prediction$predict == 1, tissue)
    saveRDS(prediction, file = p("prediction_model/prediction_", preName, tissue, ".rds"))
  }
  if (h2o.started) h2o::h2o.shutdown(prompt = FALSE)
  return(makeCombinedPrediction(tissues, mode = mode))
}

#' The training model with cds and 3' UTRs as random forest
#' @param tissue Tissues to train on, use "combined" if you want all in one
#' @inheritParams predictUorfs
trainClassifier <- function(tissue,
                            ip = "localhost",
                            port = 54321,
                            nthreads_h2o = max(min(45, as.integer(BiocParallel::bpparam()$workers/2)), 1),
                            max_mem_size = "200G") {

  if(file.exists(paste0("prediction_model/randomForrest_",tissue))) {
    #TODO! Double check this one
    h2o.init(ip = ip, port = port, nthreads = nthreads_h2o, max_mem_size = max_mem_size)
    forrest <- h2o.loadModel(path = p("prediction_model/randomForrest_",tissue,"/",
                                          list.files(p("prediction_model/randomForrest_",tissue)[1])))
    return(forrest)
  }
  # define training control
  training <- makeTrainingData(tissue)
  forrest <- forest(training, ntrees = 100, ip = ip, port = port,
                    nthreads_h2o = nthreads_h2o, max_mem_size = max_mem_size)
  if (!is.null(tissue)) {
    h2o.saveModel(forrest, path = p("prediction_model/randomForrest_",tissue))
  }
  return(forrest)
}

#' Alternative ORF sequence classifier
#' @param prediction data.table of predictions
#' @param tissue Tissues to train on, use all if you want all in one
sequenceClassifier <- function(prediction, tissue){
  print("predicting sequence classifier")

  uorfData <- getAllSequenceFeaturesTable()

  # dt <- uorfData[!overCDS()]
  dt <- uorfData
  dt[,y := as.factor(prediction$p1 > 0.55)]

  # table(dt$StartCodons)
  # make classification
  forestH2o <- forest(dt, cv = 5, ntrees = 150)
  print(forestH2o@model)
  h2o.saveModel(forestH2o, path = p("prediction_model/finalForest_",tissue))
  # prediction
  uorfPrediction <- as.data.table(h2o.predict(forestH2o, newdata = as.h2o(uorfData)))

  saveRDS(uorfPrediction, file = p("prediction_model/seqPrediction_",tissue, ".rds"))

  # checking
  hits <- as.logical(uorfPrediction[,3] > 0.50)
  startCodonMetrics(hits)
  return(NULL)
}
