#' grl A GRangesList
#' bedName is bed name
#' uscs bed 12 format
#' @param grl a grl
#' @param bedName output name
#' @param fixChromoNaming default FALSE
#' @export
bed12 <- function(grl, bedName, fixChromoNaming = FALSE){
  # TODO, check seqlevels style should be forced to UCSC ?
  if(!ORFik:::is.grl(class(grl))) stop("grl, must be of class GRangesList")
  if (fixChromoNaming) print(seqlevels(grl))
  grl <- sortPerGroup(grl,ignore.strand = T) # <- sort bed way!

  dt.grl <- data.table(seqnames = ORFik:::seqnamesPerGroup(grl, F))
  dt.grl$start <- as.integer(ORFik:::firstStartPerGroup(grl,keep.names = F) -1)
  dt.grl$end <- ORFik:::lastExonEndPerGroup(grl, keep.names = F) #non inclusive end
  dt.grl$name <- names(grl)
  dt.grl$score <- widthPerGroup(grl, keep.names = F)
  dt.grl$strand <- strandPerGroup(grl, F)
  dt.grl$thickStart <- dt.grl$start
  dt.grl$thickEnd <- dt.grl$end
  dt.grl$rgb <- rep(0, length(grl))
  dt.grl$blockCount <- ORFik:::numExonsPerGroup(grl)
  blockSizes <- paste(width(grl), collapse = ",")
  names(blockSizes) <- NULL
  dt.grl$blockSizes <- blockSizes
  relativeStarts <- (start(grl) -1) - dt.grl$start
  blockStarts <- paste(relativeStarts, collapse = ",")
  names(blockStarts) <- NULL
  dt.grl$blockStarts <- blockStarts

  #chromStart + chromStarts[last] + blockSizes[last])
  #must equal chromEnd.
  data.table::fwrite(x = dt.grl, file = bedName,
                     sep = "\t", col.names = F, row.names = F, quote = F)
  return(NULL)
}
