#' shortcut for paste0()
#'
p <- paste0

#' Remove infinite and NAs in matrix
#' @param predicate a matrix
fixNAandINFTable <- function(predicate){
  if(!all((c("RSS", "entropyRFP") %in% colnames(predicate)))) {
    message("predicate did not contain names RSS and entropyRFP")
    return(as.data.table(predicate))
  }

  isNA <- is.na(predicate[,"RSS"])[, 1] #RSS
  predicate <- as.matrix(predicate)
  isINF <- is.infinite(predicate[, "entropyRFP"])
  predicate[isNA, "RSS"] <- 0
  predicate[isINF, "entropyRFP"] <- 0
  return(as.data.table(predicate))
}
