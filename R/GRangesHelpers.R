#' Create GRangesList of ORF IDs
#' @param uniqueIDs a character vector of ORF ids, or single column data.table
#' @return GRangesList
#' @importFrom data.table fread
#' @importFrom data.table tstrsplit
#' @importFrom data.table rbindlist
toGRFromUniqueID <- function(uniqueIDs) {

  if (class(uniqueIDs)[1] != "data.table"){
    if (class(uniqueIDs)[1] == "character"){
      uniqueIDs <- data.table(uorfID = uniqueIDs)
    } else stop("uniqueIDs must either be data.table or character!")
  }else if (ncol(uniqueIDs) != 1) stop("must only be 1 column")
  colnames(uniqueIDs) = "uorfID"

  splitList <- fread(text = uniqueIDs$uorfID, sep = ",", header = FALSE)
  if (ncol(splitList) != 3) stop("not correct ncols in uniqueIDs")
  # Set strand and seqnames (col 1 and 2)
  seqnamesUsed <- as.character(splitList[[1]])
  strands <- as.character(splitList[[2]])

  # Split on exons
  colnames(splitList) <- paste0("uorfID_", 1:3)
  numExons <- nchar(gsub("[^;]*", "", splitList$uorfID_3))
  whichMax <- which.max(numExons)

  splitList <- rbindlist(list(splitList[whichMax,],splitList))
  a <- fread(text = splitList$uorfID_3, header = FALSE, fill = TRUE, sep = ";")
  a <- a[-1, ]

  exonsList <- as.data.table(matrix(ncol = ncol(a)*2, nrow = nrow(a)))
  j = 1
  for(i in colnames(a)) {
    exons <- data.table:::tstrsplit(a[, get(i)], split = " ")
    exonsList[,j] <- exons[[1]]
    exonsList[,ncol(a)+j] <- exons[[2]]
    j <- j + 1
  }

  starts <- exonsList[, 1:ncol(a)]
  widths <- exonsList[, (ncol(a)+1):(ncol(a)*2)]
  counts <- rowSums(!is.na(starts))
  grouping <- unlist(lapply(seq.int(length(counts)), function(x) {
    rep.int(x, counts[x])
  }))

  #IntegerLists
  starts <- c(t(starts))
  starts <- as.integer(starts[!is.na(starts)])
  widths <- c(t(widths))
  widths <- as.integer(widths[!is.na(widths)])

  gr = GRanges(seqnames = seqnamesUsed[grouping],
               ranges = IRanges(start = starts, width = widths),
               strand = strands[grouping])

  grl = groupGRangesBy(gr, grouping)
  grl <- sortPerGroup(grl)

  #test the orfs
  widths <- ORFik:::widthPerGroup(grl, F)
  if(!all((widths %% 3) == 0)) stop("widths of uorfs are not %3 = 0!")
  if (length(grl) != nrow(uniqueIDs)) stop("not all ranges was reconstructed properly, check data!")
  return(grl)
}
