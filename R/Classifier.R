#' Predict uORFs using random forrest classification
#'
#' Training data
#' pos set is by cds
#' Neg seg is by 3'utrs
#'
#' @param tissues Tissues to train on, use all if you want all in one
#' @param nthreads_h2o number of cores for H20: default max(45, detectCores()/3)
#' @param max_mem_size max allowed memory for H20: default ("200G")
#' @export
#' @import h2o
predictUorfs <- function(tissues, nthreads_h2o = max(45, detectCores()/3),
                         max_mem_size = "200G") {
  for (tissue in tissues) {
    message(p("Prediction for tissue: ", tissue))
    # make uORFTable
    if(file.exists(paste0("prediction_model/prediction_", tissue, ".rds"))) {
      prediction <- readRDS(paste0("prediction_model/prediction_", tissue, ".rds"))
    } else {
      forestRibo <- uORFomePipe:::trainClassifier(tissue, nthreads_h2o, max_mem_size)
      prediction <- as.data.table(h2o.predict(forestRibo,
                                              as.h2o(uORFomePipe:::makeORFPredictionData())))
      hits <- as.logical(prediction[,3] > 0.50)
      uORFomePipe:::startCodonMetrics(hits)
      saveRDS(prediction, file = p("prediction_model/prediction_", tissue, ".rds"))
    }
  }
  # TODO: output bed12 with colors
  return(makeCombinedPrediction(tissues))
}

#' The training model with cds and 3' UTRs as random forest
#' @param tissue Tissues to train on, use all if you want all in one
trainClassifier <- function(tissue = NULL, nthreads, max_mem_size) {

  if(file.exists(paste0("prediction_model/randomForrest_",tissue))) {
    forrest <- h2o.loadModel(path = p("prediction_model/randomForrest_",tissue,"/",
                                          list.files(p("prediction_model/randomForrest_",tissue)[1])))
    return(forrest)
  }
  # define training control
  training <- makeTrainingData(tissue)
  forrest <- forest(training, ntrees = 100, nthreads = nthreads, max_mem_size = max_mem_size)
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
