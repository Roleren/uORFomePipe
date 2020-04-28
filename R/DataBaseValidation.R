#' Plot top tissues
#'
#' Shows max 20 tissues
#' @import ggplot2
#' @import gridExtra
#' @export
predictionVsCageHits <- function(saveName){
  cageTissues <- readTable("tissueAtlasByCage", with.IDs = FALSE)

  cageRed <- cageTissues[, colSums(cageTissues) > 0, with = F]
  inAll <- sum(rowSums(cageRed) == ncol(cageRed))
  bestNames <- names(sort(-colSums(cageRed)))[1:min(20, ncol(cageRed))]
  cageRed <- cageRed[,bestNames, with = F]
  values <- c(colSums(cageRed) - inAll, rep(inAll, ncol(cageRed)))
  variable <- c(rep(colnames(cageRed), 2))
  type <- c(rep("Unique uORFs tissue", ncol(cageRed)), rep("uORFs in all tissues", ncol(cageRed)))
  df <- data.table(value = values, variable, type)
  df <- df[order(value),]
  df$variable <- factor(df$variable, levels = unique(df$variable), ordered = T)
  cageAll <- ggplot(df, aes(x=variable,y=value,fill=type)) +
    geom_bar(stat="identity", position="stack") +
    xlab("Tissue")+ylab("# uORFs found by CAGE") +
    theme(axis.text.y = element_text(size = 10)) +
    guides(fill=FALSE) +
    coord_flip()

  # for prediction
  cageTissuesPrediction <- readTable("tissueAtlasByCageAndPred", with.IDs = FALSE)
  inAll <- sum(rowSums(cageTissuesPrediction) == ncol(cageTissuesPrediction))
  bestNames <- names(sort(-colSums(cageTissuesPrediction)))[1:min(20, ncol(cageRed))]

  cageRed <- cageTissuesPrediction[,bestNames, with = F]
  values <- c(colSums(cageRed) - inAll, rep(inAll, ncol(cageRed)))
  variable <- c(rep(colnames(cageRed), 2))
  type <- c(rep("Unique uORFs tissue", ncol(cageRed)), rep("uORFs in all tissues", ncol(cageRed)))
  df <- data.table(value = values, variable, type)
  df <- df[order(value),]
  df$variable <- factor(df$variable, levels = unique(df$variable), ordered = T)
  predAll <- ggplot(df, aes(x=variable, y=value, fill=type)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_discrete(name="uORF counts in final prediction") +
    xlab("Tissue")+ylab("# predicted active uORFs") +
    guides(fill=FALSE) +
    theme(axis.text.y = element_text(size = 10)) +
    coord_flip()
  pred <- gridExtra::grid.arrange(top = "uORF prediction", cageAll, predAll, ncol = 2)
  codons <- startAndStopCodonPlots()
  gridExtra::grid.arrange(pred, codons, ncol = 1)

  return(invisible(NULL))
}

#' @import gridExtra
startAndStopCodonPlots <- function() {
  cageTissuesPrediction <- readTable("finalCAGEuORFPrediction")
  uorfData <- getAllSequenceFeaturesTable()

  startAndStop <- data.table(StartCodons = factor(uorfData$StartCodons),
                             StopCodons = factor(uorfData$StopCodons),
                             prediction = cageTissuesPrediction$Matrix == 1)
  startCandidates <- ggplot(data = startAndStop, aes(StartCodons)) +
    geom_bar(width = 0.3)
  startPredicted <- ggplot(data = startAndStop[prediction == TRUE,], aes(StartCodons)) +
    geom_bar(width = 0.3)
  stopCandidates <- ggplot(data = startAndStop, aes(StopCodons)) +
    geom_bar(width = 0.3)
  stopPredicted <- ggplot(data = startAndStop[prediction == TRUE,], aes(StopCodons)) +
    geom_bar(width = 0.3)
  return(gridExtra::grid.arrange(startCandidates, startPredicted, stopCandidates,
                          stopPredicted, ncol = 2))
}

#' Venn diagram between two tissues
#' @import VennDiagram
varianceTissueUsage <- function(tissue1, tissue2, pred = readTable("tissueAtlasByCageAndPred", with.IDs = FALSE)) {
  stop("not working yet")
  tab <- table(pred[,get(tissue1)], pred[,get(tissue2)])
  chi <- chisq.test(tab)
  chi$residuals

  grid.newpage()
  boOver <- draw.pairwise.venn(15021, 14213,
                               11523, category = c("Ovary", "Brain"),
                               lty = rep("blank", 2), fill = c("red", "light blue"),
                               alpha = rep(0.5, 2), cat.pos = c(0, 0),
                               cat.dist = rep(0.025, 2), title = "abc")
  boOver <- grid.arrange(gTree(children=boOver), top=textGrob("uORF overlaps", gp=gpar(fontsize=20,font=8)),
                         bottom="")

  library(cowplot)
  plot_grid(predAll,boOver, align='hv',nrow=1,labels=c('A','B'))
}



#' Validate uorfs of data-base
#'
#' This is a check to see that pipeline have done everything correctly
#' if redoing the findOverlaps does not find all orfs within fiveUTRs
#' it means that some orfs are outside the mapping area
#' this should not happen!
validateExperiments <- function(grl){

  fiveUTRs <- leaderAllSpanning()
  a <- findOverlaps(query = unlist(grl, use.names = F), fiveUTRs)
  a <- a[!duplicated(from(a))]
  if(length(a) != length(unlist(grl))){
    stop("Not all orfs where within the FiveUTRs used
         to make them, something is wrong!")
  } else print("experiments look valid")
}

#' check that all orfs acctually have a valid start codon
#' A good check for minus strand errors.
validateStartCodons <- function(uniqueIDs, startCodons){
  gr <- toGRFromUniqueID(uniqueIDs)
  names(gr) <- seq(length(gr))
  g <- unlist(gr, use.names = TRUE)
  gr <- groupGRangesBy(gr)

  st <- startCodons(gr)
  getSequencesFromFasta(st, T)

  starts <- unlist(str_split(startCodons, pattern = "\\|"))
  temp <- 0
  for( codon in starts) {
    temp <- temp + sum(seqs == codon)
  }

  return(temp == length(gr))
}

#' Get variance between different leader versions
#'
#' For validation
#' @import ggplot2
getAllLeaderChanges <- function(){
  getLeaders()
  widths <- widthPerGroup(fiveUTRs)

  leadersList = list.files(leadersFolder, full.names = TRUE)
  output <- bplapply(leadersList, function(x, widths) {
    widthsCage <- ORFik:::widthPerGroup(readRDS(i))
    diffWidths <- widths - widthsCage
    same <- sum(diffWidths == 0)
    bigger <- sum(diffWidths < 0)
    smaller <- sum(diffWidths > 0)
    meanDifBigger <- mean(diffWidths[diffWidths < 0])
    meanDifSmaller <- mean(diffWidths[diffWidths > 0])
    return(c(same,bigger,smaller,meanDifBigger,meanDifSmaller))
  }, widths = widths)
  output <- setDT(unlist(output, recursive = FALSE))

  dt <- as.data.table(matrix(output, ncol = 5))
  colnames(dt) <- c("same", "bigger", "smaller", "meanBigger", "meanSmaller")
  save(dt, file = "leaderWidthChanges.rdata")

  # as box plot per tissue, 1 box per, width change
  meanFiveWidth <- mean(widths)

  changes <- dt
  changes$widthChange <- -((changes$same * meanFiveWidth) + (changes$bigger * changes$meanBigger) + (changes$smaller * changes$meanSmaller))

  cageTable <- getCageInfoTable()
  cageWeHave <- getAllUsableCage(cageTable)


  changes <- changes[cageWeHave$cage_index,]
  changes$tissue <- cageWeHave$Characteristics.Tissue.
  changes$widthChangePer <- changes$widthChange/length(fiveUTRs)
  boxplot(changes$widthChange, fill = changes$tissue)

  ggplot(changes, aes(x=as.factor(tissue),y=widthChangePer-mean(widthChangePer)), fill = as.factor(tissue))+
    geom_boxplot() +
    coord_flip() +
    geom_hline(aes(yintercept = 0, colour = "red")) +
    theme(axis.text=element_text(size=5)) +
    xlab("Tissue") + ylab("change of widths leaders")

  # new test, check per leader change
  changes <- changes[, cageWeHave$cage_index, with = F]

  setwd(dataBaseFolder)
  save(changes, file = "leaderWidthChangesPerLeader.rdata")


  tissue <- gsub(" ", ".", cageWeHave$Characteristics.Tissue.)
  uniqueTissues <- unique(tissue)
  getCDS()
  cdsLength <- widthPerGroup(firstExonPerGroup(cds[names(fiveUTRs)]), keep.names = F)
  c <- as.data.table(matrix(nrow = nrow(changes), ncol = 1))
  for(i in 1:(length(uniqueTissues))){
    cageFilestoCheck <- cageWeHave[which(tissue == uniqueTissues[i]),]$cage_index
    c[, uniqueTissues[i] := rowMeans(changes[,paste0("result.",cageFilestoCheck), with = F]) - cdsLength]
  }
  c$V1 <- NULL

  cageTissues <- readTable("tissueAtlasByCage")

  whichColNames <- which(colSums(cageTissues) > 0)
  c <- c[, whichColNames, with = F]

  uniqueTissues <- colnames(c)
  cc <- matrix(unlist(c, use.names = F), ncol = 1, nrow = nrow(c)*ncol(c))


  # test
  tissueExpand <- unlist(lapply(uniqueTissues, function(x) rep(x, nrow(c))), use.names = F)

  df <- data.table(value = cc, variable = tissueExpand)
  colnames(df) <- c("value", "variable")
  df$value <- -df$value

  res <- ggplot(df, aes(x=variable,y=value))+
    geom_boxplot(alpha = 0.7) +
    coord_flip() +
    geom_hline(aes(yintercept = 0), colour = "red") +
    theme(axis.text.y = element_text(size=6)) +
    xlab("Tissue") + ylab("change in leader widths")
  plot(res)
  # for cds te variance number of uORFs in changed te value tx
  # take brain adoult vs glioblastoma
  # fold change categories ( 4 types (both, 1st or 2nd, none))
  # type will geom_ecdf (4 distribution)
  # data.frame with foldchange, id, type(4))

  c <- as.matrix(c)
  # log10(new leader length / old leader length )
  a <- rowMeans(c, na.rm = T)
  df <- data.frame(var = a)
  res <- ggplot(df, aes(x = var))+
    geom_density(fill = "blue") +
    xlim(-500, 500) +
    theme(axis.text.x = element_text(angle=0, size = 10), axis.text.y = element_text(size = 12),
          plot.title = element_text(size = 12), plot.subtitle = element_text(hjust = 0.5)) +
    labs(title="Mean change of leader lengths",
         subtitle="Over all tissues by CAGE", x = "Average change (bp length)", y = "proportion")
  plot(res)


  #library(cowplot)
  #plot_grid(,align='hv',nrow=1,labels=c('A','B'))
  return(gridExtra::grid.arrange(cageAll, res, nrow = 1))
}

#' validate features using
#' seperate them using pca and quantiles
pcaCAGEValidation <- function(uids1, uids2) {
  stop("not ready yet!")
  # goal: what is different between groups and within group
  # variance / aov
  # clustering, all features hclust(dists) , kmean, jaccard index
  # scale to normalize
  # which are important, pca, svd etc.
  # get uorf names , do for both, change standardCage

  # now merge, either all, or cross set

  dtUid1 <- data.table(uid = uids1, index = 1:length(uids1))
  dtUid2 <- data.table(uid = uids2, index = 1:length(uids2))
  # cross valid set
  dtMerged <- merge(x = dtUid1, y = dtUid2, by = "uid")
  # all
  #dtMerged <- merge(x = dtUid1, y = dtUid2, by = "uid", all = T)

  # now get data
  experiments <- grep(x = list.files(getwd()), pattern = "dt_", value = T)
  experiments <- paste0(getwd(),"/", experiments)
  for(i in experiments){
    assign(x = paste0("dt", i), data.table::fread(i), envir = .GlobalEnv)
  }

  # fix na, inf, nan values, set to 0
  dtExper <- paste0("dt", experiments)
  for(i in dtExper){
    invisible(lapply(names(get(i)),function(.name) set(get(i), which(!is.finite(get(i)[[.name]])), j = .name,value =0)))
  }
  # now use merge to make nrows equal for both lists
  for(i in dtExper[1:5]){
    assign(i, get(i)[dtMerged$index.x,], .GlobalEnv)
    invisible(lapply(names(get(i)),function(.name) set(get(i), which(!is.finite(get(i)[[.name]])), j = .name,value =0)))
  }
  for(i in dtExper[6:10]){
    assign(i, get(i)[dtMerged$index.y,], .GlobalEnv)
    invisible(lapply(names(get(i)),function(.name) set(get(i), which(!is.finite(get(i)[[.name]])), j = .name,value =0)))
  }


  # get a data.table for each feature
  ncols <- ncol(get(dtExper[1]))
  nrows <- nrow(get(dtExper[1]))
  names <- colnames(get(dtExper[1]))

  for(i in names){
    assign(i, data.table(matrix(nrow = nrows, ncol = 10)), .GlobalEnv)
  }
  colnames <- colnames(get(i))
  # for each column data.table(ig. floss data.table), update each column
  x <- 1L
  for(i in names){

    currentCol <- get(i)
    y <- 1L
    for(j in dtExper){
      currentCol[, colnames[y] :=  get(j)[, x, with = F]]
      y <- y + 1L
    }
    assign(names[x], currentCol, .GlobalEnv)
    x <- x + 1L
  }
  # scale them ?

  # now do pca:
  for(i in names[2:length(names)]){
    assign(paste0("pca_",i), prcomp(get(i), scale = T), .GlobalEnv)
  }

  # make the pca, do the plots, make quantiles,
  pcaNames <- paste0("pca_",names[2:length(names)])
  for(i in pcaNames){
    pca <- get(i)
    barplot(pca$sdev/sum(pca$sdev), main = i) # variance in percentage

    # split tissue groups
    plot(pca$rotation, col = c(rep(1, 5), rep(2, 5)), main = i)
    abline(h = mean(pca$rotation[,2]))
    abline(v = mean(pca$rotation[,1]))

    seperator = pca$x[,1]
    nfquant <- quantile(seperator, 0.95)
    zfquant <- quantile(seperator, 0.05)

    fivepercentSet <- seperator <= zfquant
    ninetyFivePercentSet <- seperator >= nfquant

    seperator = pca$x[, 2]
    nfquant <- quantile(seperator, 0.95)
    zfquant <- quantile(seperator, 0.05)

    fivepercentSet <- fivepercentSet | (seperator <= zfquant)
    ninetyFivePercentSet <- ninetyFivePercentSet | (seperator >= nfquant)

    fivepercentWhich <- which(fivepercentSet)
    ninetyFiveWhich <- which(ninetyFivePercentSet)
    assign(paste0("fq", i), fivepercentWhich, .GlobalEnv)
    assign(paste0("nfq", i), ninetyFiveWhich, .GlobalEnv)
  }
  # find the quntile overlaps of uorfs, all duplicated
  # these are the potentialy interesting uorfs
  # look for overlaps between the quantile sets
  # remove booleans
  # which features are good, which uorfs are regulated
  pcaNamesNonBool <- pcaNames[-c(9, 11, 12, 13, 14)]

  nfqNames <- paste0("nfq", pcaNamesNonBool)

  fqNames <- paste0("fq", pcaNamesNonBool)
  tempMatches <- NULL
  for(i in 1:length(fqNames)){
    tempMatches <- c(tempMatches, get(fqNames[i]))
  }
  oriMatches <- tempMatches
  oriLengthMatches <- length(oriMatches)
  fquniqueMatches <- unique(tempMatches)
  lengthUniqueMatches <- length(fquniqueMatches)
  fqdiffMatches <-lengthUniqueMatches / oriLengthMatches

  tempMatches <- NULL
  for(i in 1:length(nfqNames)){
    tempMatches <- c(tempMatches, get(nfqNames[i]))
  }

  # cluster
  for(i in names){
    clusterUorfFeature(get(i), inBoth, paste0("./clustering/","cluster_",i,".pdf"))
  }

}

bestTeUorfsOnTxEffect <- function(){
  uorfTEs <- readTable("teFiltered", with.IDs = T)

  #maxTEs <- uorfTEs[, lapply(.SD, max), by = txNames]

  rowMeansUorfs <- rowMeans(uorfTEs[,3:ncol(uorfTEs)])
  q90 <- quantile(rowMeansUorfs, 0.90)
  best <- which(rowMeansUorfs > q90)
  bestNames <- uorfTEs[best, txNames]

  numberOfUorfsPerTx <- readTable("numberOfUorfsPerTx")

  whichNames <- numberOfUorfsPerTx$txNames %in% bestNames
  meanNumberBest <- mean(numberOfUorfsPerTx$nUorfs[whichNames])

  # worst
  q10 <- quantile(rowMeansUorfs, 0.10)
  worst <- which(rowMeansUorfs < q10)
  worstNames <- uorfTEs[worst, txNames]
  whichNames <- numberOfUorfsPerTx$txNames %in% worstNames

  meanNumberWorst <- mean(numberOfUorfsPerTx$nUorfs[whichNames])

  # so this is not a good feature by itself
  cdsTEs <- readTable("cdsTeFiltered", with.IDs = T)
  rowMeansCDS <- rowMeans(cdsTEs[,2:ncol(cdsTEs)])

  q90 <- quantile(rowMeansCDS, 0.90)
  best <- which(rowMeansCDS > q90)
  bestNames <- cdsTEs[best, txNames]

  whichNames <- numberOfUorfsPerTx$txNames %in% bestNames
  meanNumberBest <- mean(numberOfUorfsPerTx$nUorfs[whichNames])

  # worst
  q10 <- quantile(rowMeansCDS, 0.10)
  worst <- which(rowMeansCDS < q10)
  worstNames <- cdsTEs[worst, txNames]
  whichNames <- numberOfUorfsPerTx$txNames %in% worstNames

  meanNumberWorst <- mean(numberOfUorfsPerTx$nUorfs[whichNames])

  # number of uORFs seem to correlate with TE of cds
}
