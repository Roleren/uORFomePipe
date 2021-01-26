#' shortcut for paste0()
#' @inheritDotParams base::paste
p <- paste0

#' Remove infinite and NAs in matrix
#' @param predicate a matrix
fixNAandINFTable <- function(predicate){
  if(!all((c("RSS", "RRS","entropyRFP") %in% colnames(predicate)))) {
    message("predicate did not contain names RSS, RRS and entropyRFP")
    return(as.data.table(predicate))
  }

  isNA <- is.na(predicate[,"RSS"])[, 1] #RSS
  isNA.RRS <- is.na(predicate[,"RRS"])[, 1] #RRS
  predicate <- as.matrix(predicate)
  isINF <- is.infinite(predicate[, "entropyRFP"])
  predicate[isNA, "RSS"] <- 0
  predicate[isNA.RRS, "RSS"] <- 0
  predicate[isINF, "entropyRFP"] <- 0
  return(as.data.table(predicate))
}

#' Artificial ORFs from CDS
#'
#' A wrapper to ORFik:::artificial.orfs, that create artificial ORFs from CDS,
#' with given start and stop codons of size size * 6
#' @param cds GRangesList of CDS
#' @param size integer, default: 100, 1/6 of maximum size of ORFs (max size 600 if 100)
#' @param startCodons, character vector, default: c("ATG", "CTG", "TTG", "AAG", "AGG")
#' @param stopCodons character vector, default: c("TAA", "TGA", "TAG")
#' @param fa a FaFile or path to fasta index file
#' @param mode character, default: "aCDS", must be either CDS or aCDS (artificial CDS)
#' @return GRangesList of new ORFs
artificial.from.cds <- function(cds, size = 100,
                                startCodons = c("ATG", "CTG", "TTG", "AAG", "AGG"),
                                stopCodons = c("TAA", "TGA", "TAG"),
                                fa,
                                mode = "aCDS") {
  if (!(mode %in% c("CDS", "aCDS"))) stop("mode must be either CDS or aCDS (artificial CDS)")
  orfs <- subset.ORF(cds, startCodons = startCodons,
                    stopCodons = stopCodons, fa = fa)
  names(orfs) <- paste0(names(orfs), "_1")
  if (mode == "aCDS") {
    message("Using artificial CDS' as sequences")
    orfs <- ORFik:::artificial.orfs(orfs, end5 = size, start3 = -size)
  } else {
    message("Using CDS' as sequences")
  }
  return(orfs)
}
#' Subset orfs by filters given
#' @inheritParams artificial.from.cds
#' @param orfs GRangesList of ORFs
#' @param minimum.length numeric, default 6
#' @param must.be.mod.3.0 logical, default TRUE
#' @return GRangesList of subset that is valid by filter
subset.ORF <- function(orfs, minimum.length = 6,
                       must.be.mod.3.0 = TRUE,
                       startCodons = c("ATG", "CTG", "TTG"),
                       stopCodons = c("TAA", "TGA", "TAG"),
                       fa) {
  if (!is(minimum.length, "numeric")) stop("minimum.length must be numeric!")
  # If regex notation, switch it to character vector style
  if (length(startCodons) == 1 & length(grep("|", startCodons)))
    startCodons <- unlist(strsplit(startCodons, split = "\\|"))
  if (length(stopCodons) == 1 & length(grep("|", stopCodons)))
    stopCodons <- unlist(strsplit(stopCodons, split = "\\|"))
  # Length & modulus filter
  orfs <- orfs[widthPerGroup(orfs) >= minimum.length]
  if (must.be.mod.3.0)
    orfs <- orfs[(widthPerGroup(orfs, F) %% 3) == 0]

  orfs <- ORFik:::removeMetaCols(orfs)
  starts <- ORFik:::txSeqsFromFa(startCodons(orfs), fa, TRUE, FALSE)
  stops <- ORFik:::txSeqsFromFa(stopCodons(orfs), fa, TRUE, FALSE)
  return(orfs[(starts %in% startCodons) & (stops %in% stopCodons)])
}
