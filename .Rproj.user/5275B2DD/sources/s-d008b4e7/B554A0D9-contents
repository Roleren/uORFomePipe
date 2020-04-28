#' Comparison with uORFome of McGilivray et al.
#'
#' Uses the supplemental data from article
#'
#' Remap tool:
#' https://www.ncbi.nlm.nih.gov/genome/tools/remap
#' @importFrom openxlsx read.xlsx
verifyOthersWorks <- function(work = "/export/valenfs/projects/uORFome/Supplemental_Data_Tables_.xlsx") {

  # old grch37 annotation
  d <- read.xlsx(work, sheet = 6, colNames = T, startRow = 3)

  strand <- d$strand
  startSite <- as.integer(d$start_coordinate)
  stopSite <- as.integer(d$end_coordinate)
  transcript <- sub("\\..*","",d$uORF_ID)
  chromosome <- d$chromosome

  sta <- startSite
  sta[strand == "-"] <- stopSite[strand == "-"]
  sto <- stopSite
  sto[strand == "-"] <- startSite[strand == "-"]
  getLeaders()
  goodHits <- which(transcript %in% names(fiveUTRs))
  bed <- GRanges(seqnames = chromosome, ranges = IRanges(sta, sto), strand = strand, name = transcript, score = rep(0, length(chromosome)))

  rtracklayer::export.bed(object = bed, con = "mcgilivray_uorfs.bed") # now do conversion
  print("@ NCBI Genome Remapping Service GrCH37 -> 38")
  print("File: mcgilivray_uorfs.bed")
  print("Go to: https://www.ncbi.nlm.nih.gov/genome/tools/remap")
  my.name <- readline(prompt="Press enter when you finished remapping and downloaded annotation data: ")
  g <- rtracklayer::import.bed(con = "remapped_mcgilivray_uorfs.bed")
  g$score <- NULL
  g <- groupGRangesBy(g, seq(length(g)))
  if( length(strand) != length(startSite) | length(strand) != length(transcript)){
    stop("wrong input readings for others")
  }
  # Test their start codons:
  getFasta()
  codons <- ORFik:::startCodons(g, is.sorted = T)
  startcods <- ORFik:::txSeqsFromFa(codons[seqnamesPerGroup(codons, F) %in% seqlevels(fa)], fa, is.sorted = TRUE)
  table(as.character(startcods, use.names = FALSE))


  grl <- getUorfsInDb()

  #how many did we find ?
  startsOur <- startSites(grl, is.sorted = T, keep.names = FALSE)
  stopsOur <- stopSites(grl, is.sorted = T, keep.names = FALSE)
  strandOur <- strandPerGroup(grl, keep.names = F)
  chromosomeOurs <- seqnamesPerGroup(grl, F)

  ours <- paste(startsOur, stopsOur, chromosomeOurs, strandOur)
  theirs <- paste(startSites(g, is.sorted = T, keep.names = FALSE),
                  stopSites(g, is.sorted = T, keep.names = FALSE),
                  seqnamesPerGroup(g, F), strandPerGroup(g, keep.names = F))
  hitsOurs <- which(ours %in% theirs)
  hitsTheirs <- which(theirs %in% ours)

  print(paste("Mapping step found from theirs, total uORFs of:", length(hitsOurs), "that is:", length(hitsOurs) /length(g), "of full set"))
  #how many match our prediction

  # Using only predicted uORFs
  predicted <- which(readTable("finalCAGEuORFPrediction")$Matrix == 1)
  ours <- ours[predicted]

  hitsOurs <- which(ours %in% theirs)
  hitsTheirs <- which(theirs %in% ours )
  print(paste("Mapping step found from theirs, total uORFs of:", length(hitsOurs), "that is:", length(hitsOurs) /length(g), "of full set"))

  #ratio from ours to theirs

  ratio <- 1675/645 # ours by theirs to get real number of hits
  ratio <- length(hitsOurs)
  theirSearchSpace <- length(theirs)
  ourSearchSpace <- length(ours)
  library(grid)
  grid.newpage()
  library(VennDiagram)
  ggplot <- draw.pairwise.venn(ourSearchSpace,
                               theirSearchSpace,
                               ratio, category = c("uORFome prediction", "McGillivray et al. prediction"),
                               lty = rep("blank", 2), fill = c("cyan", "red"),
                               alpha = rep(0.5, 2), cat.pos = c(0, 0),
                               cat.dist = rep(0.025, 2), cat.cex = c(1,1), cex = c(1,1,1))
  library(gridExtra)
  vennPred <- grid.arrange(gTree(children=ggplot), top=textGrob("Overlap between prediction pipelines", gp=gpar(fontsize=20,font=8)),
                           bottom="")

  # Start codon metrics
  uorfData <- getAllSequenceFeaturesTable()
  StartCodons <- uorfData$StartCodons
  table(StartCodons[hitsOurs])
  #table(StartCodons[uniqueOrder][finalCagePred])
  table(StartCodons[predicted])
  tab1 <- table(StartCodons[hitsOurs])/sum(table(StartCodons[hitsOurs]))
  #tab2 <- table(StartCodons$startCodon[uniqueOrder][finalCagePred])/sum(table(StartCodons$startCodon[uniqueOrder][finalCagePred]))
  tab2 <- table(StartCodons[predicted])/sum(table(StartCodons[predicted]))

  df <- data.frame(value = c(tab1, tab2), variable =c(names(tab1), names(tab2)),
                   pred  = c(rep("McGillivray et al. prediction", length(tab1)), rep("uORFome prediction", length(tab2))))

  cstarts <- ggplot(df, aes(x=variable,y=value,fill=factor(pred)))+
    geom_bar(stat="identity",position="dodge")+
    scale_fill_discrete(name="Prediction pipeline")+
    xlab("Start codon")+ylab("percentage")

  library(cowplot)
  plot_grid(vennPred,cstarts, align='hv',nrow=2,labels=c('A','B'))
}

validateRiboRNAPairs <- function(){
  matching_rna_ribo <- getRiboRNAInfoTable()

  validate <- c()
  for(SRR in matching_rna_ribo$ribo){
    validate <- c(validate, grep(pattern = SRR ,x = list.files(rfpFolder), value = T))
  }
  if(nrow(matching_rna_ribo) != length(validate)) stop(paste0("did not find all paths for riboseq, only found ", length(validate)," of ", nrow(matching_rna_ribo)))

  validate <- c()
  for(SRR in matching_rna_ribo$rna){
    RNAPath <- grep(SRR,list.files(path = rnaFolder, all.files = TRUE, full.names = TRUE, recursive = TRUE), value=TRUE)
    if (length(RNAPath) != 1) stop(p("did not find unique RNA path for SRR: ", SRR))
    validate <- c(validate, RNAPath)
  }
  if(nrow(matching_rna_ribo) != length(validate)) stop(paste0("did not find all paths for rnaseq, only found ", length(validate)," of ", nrow(matching_rna_ribo)))

  print("Experiments validated, found all pairs of ribo-seq and rna-seq, you can continue!")
  return(NULL)
}
#' @import ggplot2
teVariance <- function(){
  if(tableNotExists("biggestTEVariance")){


    uorfTEs <- readTable("teFiltered", with.IDs = T)
    dt <- teAtlasTissueNew(uorfTEs)
    # take the different tissues, te, make ratio matrix

    ratioMatrix <- data.table(ratio = 1)
    # total ratio
    for( i in 2:ncol(dt)){
      ratioMatrix <- data.table(ratioMatrix, sum(dt[,i, with = F] > 1.1)/nrow(dt))
    }
    colnames(ratioMatrix) <- c("ratio", uniques )

    # comparison matrices

    comparisonMatrix <- data.table(txNames = dt$txNames)
    # total ratio
    for( i in 2:(ncol(dt)-1)){
      for(j in (i+1):(ncol(dt))){
        comparisonMatrix <- data.table(comparisonMatrix, abs(log10(dt[,i, with = F] / dt[,j, with = F])))
      }
    }

    comparisonMeanMatrix <- rowMeans(comparisonMatrix[,2:ncol(comparisonMatrix)])

    q90 <- quantile(comparisonMeanMatrix, 0.90)
    comparisonMeanMatrix <- data.table(uorfIDs = uorfTEs$uorfID,
                                       txNames = dt$txNames,
                                       comparisonMeanMatrix)
    colnames(comparisonMeanMatrix)[3] <- c("meanDiffTE")



    biggestTEVariance <- comparisonMeanMatrix[ comparisonMeanMatrix$meanDiffTE > q90,]
    biggestTEVariance$which <- which(comparisonMeanMatrix$meanDiffTE > q90)
    insertTable(biggestTEVariance, "biggestTEVariance", rmOld = T)

    q10 <- quantile(comparisonMeanMatrix$meanDiffTE, 0.10)
    smallestTEVariance <- comparisonMeanMatrix[ comparisonMeanMatrix$meanDiffTE < q10,]
    smallestTEVariance$which <- which(comparisonMeanMatrix$meanDiffTE < q10)
    insertTable(smallestTEVariance, "smallestTEVariance", rmOld = T)
    # presence
    # difference tissue

    dtPlot <- dt[,2:ncol(dt)]
    dtPlot$x <- as.factor(1:nrow(dtPlot))

    dt.m <- melt(dtPlot, id="x")
    plotTitle <- "Translational efficiency variance in Tissues"
    ggplot(dt.m) + geom_boxplot(aes(x = variable, log10(value))) +
      xlab("Tissue") + ylab("Te") +
      ggtitle(plotTitle)
    ggsave(plotTitle)
  }

  # cds effect of top bottom uorf variance, remove duplicated top bottom
  cdsTETissue <- readTable("cdsTETissueMean", with.IDs = T)


  cdsTETissueBig <- cdsTETissue[txNames %in% biggestTEVariance$txNames,]

  cdsTETissueSmall <- cdsTETissue[txNames %in% smallestTEVariance$txNames,]
  # filter out equal tx
  cdsTETissueBig <- cdsTETissueBig[!(txNames %in% cdsTETissueSmall$txNames),]
  cdsTETissueSmall <- cdsTETissueSmall[!(txNames %in% cdsTETissueBig$txNames),]


  colnames(cdsTETissueBig) <- paste0(colnames(cdsTETissueBig), "_high")
  cdsTETissueBig$high <- rep(T,nrow(cdsTETissueBig))
  cdsTETissueSmall$high <- rep(F,nrow(cdsTETissueSmall))
  merged <- rbindlist(list(cdsTETissueBig, cdsTETissueSmall))
  merged$x <- as.factor(1:nrow(merged))

  merged.m <- melt(merged[,2:ncol(merged)], id.vars=c("x", "high"))
  ggplot(merged.m) + geom_boxplot(aes(x = variable, log10(value), fill = high)) +
    xlab("Tissue") + ylab("Te") +
    ggtitle(plotTitle)
  # how many overlap the cds of its tx
  # check with row selection per tissue not mean, is this important ?

  # te cds vs te uorf, color the extremes


  plot(cdsTETissueBig[,2:(ncol(cdsTETissueBig)-1)], cdsTETissueSmall[,2:(ncol(cdsTETissueSmall)-1)])
  ggplot(merged.m, aes(x = variable, log10(value), fill = high)) +
    geom_dotplot()

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
  } else { print("experiments look valid")}
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
getAllLeaderChanges <- function(){
  if(!file.exists(p(dataFolder,"/leaderOriginalWidths.rdata"))){
    getLeaders()
    widths <- ORFik:::widthPerGroup(fiveUTRs)
    save(widths, file = p(dataFolder,"/leaderOriginalWidths.rdata"))
    rm(fiveUTRs)
  }

  setwd("/export/valenfs/projects/uORFome/RCode1/")
  pipelineCluster(75)
  nLeadersList = length(leadersList)
  rm(fiveUTRs)
  output <- foreach(i=1:nLeadersList, .combine = 'rbind') %dopar% {
    source("./uorfomeGeneratorHelperFunctions.R")
    leadersList = list.files(leadersFolder)

    load(p(dataFolder,"/leaderOriginalWidths.rdata"))
    load(p(leadersFolder,leadersList[i]))
    widthsCage <- ORFik:::widthPerGroup(fiveUTRs)

    diffWidths <- widths - widthsCage
    same <- sum(diffWidths == 0)
    bigger <- sum(diffWidths < 0)
    smaller <- sum(diffWidths > 0)
    meanDifBigger <- mean(diffWidths[diffWidths < 0])
    meanDifSmaller <- mean(diffWidths[diffWidths > 0])
    return(c(same,bigger,smaller,meanDifBigger,meanDifSmaller))
  }
  dt <- as.data.table(matrix(output, ncol = 5))
  colnames(dt) <- c("same", "bigger", "smaller", "meanBigger", "meanSmaller")
  stopCluster(cl)
  setwd(dataBaseFolder)
  save(dt,file = "leaderWidthChanges.rdata")

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
  setwd(codeFolder)
  output <- foreach(i=1:length(list.files(leadersFolder)), .combine = 'cbind') %dopar% {
    source("./HelperVariables.R")
    leadersList = list.files(leadersFolder)

    load(p(dataFolder,"/leaderOriginalWidths.rdata"))
    load(p(leadersFolder,leadersList[i]))

    return(ORFik:::widthPerGroup(fiveUTRs) - widths )
  }


  changes <- as.data.table(as.matrix(output))

  changes <- changes[, cageWeHave$cage_index, with = F]

  setwd(dataBaseFolder)
  save(changes, file = "leaderWidthChangesPerLeader.rdata")

  load("leaderWidthChangesPerLeader.rdata")


  tissue <- gsub(" ", ".", cageWeHave$Characteristics.Tissue.)
  uniqueTissues <- unique(tissue)
  library(ORFik)
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


  library(cowplot)
  plot_grid(cageAll, res ,align='hv',nrow=1,labels=c('A','B'))
}

uidFromCage <- function(cage = standardCage, asUID = TRUE,
                        with.transcriptNames = TRUE){

  rm(cageFiveUTRs)
  rm(fiveUTRs)
  getCDS()
  getThreeUTRs()
  getLeaders()
  cageFiveUTRs <- ORFik:::reassignTSSbyCage(fiveUTRs, standardCage, 1000, 1, cds)
  originalUorfsByTx <- getUnfilteredUORFs(cageFiveUTRs, assignRanges = F)
  gr <- unlist(originalUorfsByTx, use.names = F)
  grl <- groupGRangesBy(gr, gr$names)
  grl <- removeORFsWithinCDS(grl)

  if (!asUID) {
    return(grl)
  }
  uids <- toUniqueIDFromGR(grl)
  if (with.transcriptNames) {
    return(paste(uids, ORFik:::OrfToTxNames(grl)))
  }
  return(uids)
}

#' validate features using brain and hela to see that features
#'  seperate them
#'  using pca and quantiles
validateAllFeatures <- function(){

  #rm(list=ls())
  #setwd("/export/valenfs/projects/uORFome/RCode1/")
  #source("./DataBaseCreator.R")

  # goal: what is different between groups and within group
  # variance / aov
  # clustering, all features hclust(dists) , kmean, jaccard index
  # scale to normalize
  # which are important, pca, svd etc.
  # get uorf names , do for both, change standardCage
  uids1 <- uidFromCage()

  uids2 <- uidFromCage("./../DATA/CAGE/human/kidney%2c%20adult%2c%20pool1.CNhs10622.10017-101C8.hg38.nobarcode.ctss.bed.gz")


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
  # find overlapping uorfs between quantile sets
  # a way to find interesting uorfs
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
  oriMatches <- tempMatches
  oriLengthMatches <- length(oriMatches)
  nfquniqueMatches <- unique(tempMatches)
  lengthUniqueMatches <- length(nfquniqueMatches)
  nfqdiffMatches <-lengthUniqueMatches / oriLengthMatches

  inBoth <- fquniqueMatches[fquniqueMatches %in% nfquniqueMatches]

  onlyFirst <- fquniqueMatches[!(fquniqueMatches %in% nfquniqueMatches)]
  onlySecond <- nfquniqueMatches[!(nfquniqueMatches %in% fquniqueMatches)]


  # for FPKM RFP
  View(fpkmRFP[onlyFirst])
  View(fpkmRFP[onlySecond])
  meanFirst <- mean(colMeans(fpkmRFP[onlyFirst]))
  meanSecond <- mean(colMeans(fpkmRFP[onlySecond]))
  # for FPKM RNA
  meanFirstRNA <- mean(colMeans(fpkmRNA[onlyFirst]))
  meanSecondRNA <- mean(colMeans(fpkmRNA[onlySecond]))
  # the ones that overlap, check

  #for ioscore
  meanFirstIO <- mean(colMeans(ioScore[onlyFirst]))
  meanSecondIO <- mean(colMeans(ioScore[onlySecond]))

  #for orfsScore
  meanFirstORFScore <- mean(colMeans(ORFScores[onlyFirst]))
  meanSecondORFScore <- mean(colMeans(ORFScores[onlySecond]))

  #for fractionLengths
  meanFirstfractionLengths <- mean(colMeans(fractionLengths[onlyFirst]))
  meanSecondfractionLengths <- mean(colMeans(fractionLengths[onlySecond]))

  #for te
  meanFirstTe <- mean(colMeans(te[onlyFirst]))
  meanSecondTe <- mean(colMeans(te[onlySecond]))
  nInBoth <- sum(fquniqueMatches %in% nfquniqueMatches)
  print(paste(round(nInBoth/max(length(fquniqueMatches),
                                length(fquniqueMatches)), 2)*100,
              "% of uorfs are in both quantile sets"))

  # cluster
  for(i in names){
    clusterUorfFeature(get(i), inBoth, paste0("./clustering/","cluster_",i,".pdf"))
  }

}

# look like a good features test, the cds does not seperate well
kozakVsRiboseq <- function(){
  kozak <- as.numeric(unlist(readTable("kozak", with.IDs = FALSE), use.names = FALSE))
  q99 <- quantile(kozak, 0.99)
  best <- which(kozak > q99)
  ribo <- readTable("Ribofpkm", with.IDs = FALSE)
  bestMean <- mean(colMeans(ribo[best,]))
  worst <- which(kozak < quantile(kozak, 0.01))
  worstMean <- mean(colMeans(ribo[worst,]))

  getCDS()

  cdsFPKM <- readTable("cdsRfpFPKMs", with.IDs = F)
  kozakCDS <- readTable("cdsKozak", with.IDs = F)$kozakCDS
  bestCDS <-  which(kozakCDS > quantile(kozakCDS, 0.90))
  bestMeanCDS <- mean(colMeans(cdsFPKM[bestCDS,]))
  worstCDS <-  which(kozakCDS < quantile(kozakCDS, 0.10))
  worstMeanCDS <- mean(colMeans(cdsFPKM[worstCDS,]))
}

# Looks like the distance increases with better kozak sequence
kozakVsdistanceToCDS <- function(){
  kozak <- as.numeric(unlist(readTable("kozak", with.IDs = FALSE), use.names = FALSE))
  q99 <- quantile(kozak, 0.99)
  best <- which(kozak > q99)
  distCDS <- readTable("distORFCDS", with.IDs = FALSE)
  bestMean <- mean(colMeans(distCDS[best,]))
  worst <- which(kozak < quantile(kozak, 0.01))
  worstMean <- mean(colMeans(distCDS[worst,]))
}

# Looks like a good feature test, must check cds
kozakVsORFScores <- function(){
  kozak <- as.numeric(unlist(readTable("kozak", with.IDs = FALSE), use.names = FALSE))
  q99 <- quantile(kozak, 0.99)
  best <- which(kozak > q99)
  ORFScores <- readTable("ORFScores", with.IDs = FALSE)
  bestMean <- mean(colMeans(ORFScores[best,]))
  worst <- which(kozak < quantile(kozak, 0.01))
  worstMean <- mean(colMeans(ORFScores[worst,]))

  # now cds
  kozak <- as.numeric(unlist(readTable("cdsKozak", with.IDs = FALSE), use.names = FALSE))
  q99 <- quantile(kozak, 0.99)
  best <- which(kozak > q99)
  ORFScores <- readTable("cdsORFScore", with.IDs = FALSE)
  bestMean <- mean(rowMeans(ORFScores[best,]))
  worst <- which(kozak < quantile(kozak, 0.01))
  worstMean <- mean(rowMeans(ORFScores[worst,]))

  # Looks like a good seperator
}

#' How much does the TE go down for CDS with uorfs in tx
uorfTeVsCDSTe <- function(){
  cdsTEs <- readTable("cdsTeFiltered", with.IDs = T)
  grl <- getUorfsInDb()
  uorfNames <- unique(txNames(grl))
  uorfTXCDS <- cdsTEs$txNames %in% uorfNames
  cdsTEUORFs <- cdsTEs[uorfTXCDS, 2:ncol(cdsTEs)]
  withUorfCDS <- mean(colMeans(cdsTEUORFs))
  withoutUorfCDS <- mean(colMeans(cdsTEs[!uorfTXCDS, 2:ncol(cdsTEs)]))
  hist(withUorfCDS, withoutUorfCDS)
  boxplot(colMeans(cdsTEs[uorfTXCDS, 2:ncol(cdsTEs)]),colMeans(cdsTEs[!uorfTXCDS, 2:ncol(cdsTEs)]))
  # quantile
  cdsTxUorfs <- txNames(grl) %in% cdsTEs$txNames
  uorfTEs <- readTable("teFiltered", with.IDs = F)


  # uorfTEs <- data.table(riboFPKM[,1:2], uorfTEs)
  uorfTEs <- uorfTEs[, lapply(.SD, mean), by = txNames(grl)]
  colnames(uorfTEs)[1] <- "txNames"
  # best
  rowMeansUorfs <- rowMeans(uorfTEs[,2:ncol(uorfTEs)])
  q90 <- quantile(rowMeansUorfs, 0.90)
  best <- which(rowMeansUorfs > q90)
  bestCDSTEs <- mean(colMeans(cdsTEUORFs[best, ]))
  #worst
  q10 <- quantile(rowMeansUorfs, 0.10)
  worst <- which(rowMeansUorfs < q10)
  worstCDSTEs <- mean(colMeans(cdsTEUORFs[worst,]))
  # conclusion, no q90/q10 correlation it looks like

  rowMeansCDS <- rowMeans(cdsTEs[uorfTXCDS, 2:ncol(cdsTEs)])

  corResult <- cor.test(rowMeansUorfs, rowMeansCDS)


  return((1 - withoutUorfCDS/withUorfCDS)*100)
}

distVSCDSTe <- function(){

  cdsTEs <- readTable("cdsTeFiltered", with.IDs = T)
  dists <- readTable("distORFCDS")
  riboFPKM <- readTable("Ribofpkm")

  q90 <- quantile(dists$distORFCDS, 0.90)
  best <- which(dists$distORFCDS > q90)
  bestDists <- dists[best,]

  longestTXNames <- cdsTEs$txNames %in% riboFPKM$txNames[riboFPKM$uorfID %in% bestDists$uorfID]
  longestCDSTe <- mean(colMeans(cdsTEs[longestTXNames, 2:ncol(cdsTEs)]))

  q10 <- quantile(dists$distORFCDS, 0.10)
  worst <- which(dists$distORFCDS < q10)
  worstDists <- dists[worst,]

  shortestTXNames <- cdsTEs$txNames %in% riboFPKM$txNames[riboFPKM$uorfID %in% worstDists$uorfID]
  shortestCDSTe <- mean(colMeans(cdsTEs[shortestTXNames, 2:ncol(cdsTEs)]))
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

# Looks like the No correlation
kozakVsCDSTe <- function(){

  # mean grouped by transcript
  kozak <- readTable("kozak", with.IDs = T)
  mappingTxUorfs <- readTable("linkORFsToTx")
  kozak$txNames <- mappingTxUorfs$txNames
  kozak$uorfID <- NULL

  kozakMean <- kozak[, lapply(.SD, mean), by = txNames]
  kozakNames <- kozakMean$txNames
  kozak <- as.numeric(unlist(kozakMean$kozak, use.names = FALSE))
  cdsTEs <- readTable("cdsTeFiltered", with.IDs = T)
  rowMeansCDS <- rowMeans(cdsTEs[,2:ncol(cdsTEs)])
  q90 <- quantile(rowMeansCDS, 0.99)
  best <- which(rowMeansCDS > q90)
  bestNames <- cdsTEs[best, txNames]

  whichNames <- kozakNames %in% bestNames
  meanNumberBest <- mean(kozak[whichNames])

  #worst
  q10 <- quantile(rowMeansCDS, 0.01)
  worst <- which(rowMeansCDS < q10)
  worstNames <- cdsTEs[worst, txNames]
  whichNames <- kozakNames %in% worstNames

  meanNumberWorst <- mean(kozak[whichNames])

  # max grouped by transcript

  kozak <- readTable("kozak", with.IDs = T)
  mappingTxUorfs <- readTable("linkORFsToTx")
  kozak$txNames <- mappingTxUorfs$txNames
  kozak$uorfID <- NULL

  kozakMean <- kozak[, lapply(.SD, max), by = txNames]
  kozakNames <- kozakMean$txNames
  kozak <- as.numeric(unlist(kozakMean$kozak, use.names = FALSE))
  cdsTEs <- readTable("cdsTeFiltered", with.IDs = T)
  rowMeansCDS <- rowMeans(cdsTEs[,2:ncol(cdsTEs)])
  q90 <- quantile(rowMeansCDS, 0.99)
  best <- which(rowMeansCDS > q90)
  bestNames <- cdsTEs[best, txNames]

  whichNames <- kozakNames %in% bestNames
  meanNumberBest <- mean(kozak[whichNames])

  #worst
  q10 <- quantile(rowMeansCDS, 0.01)
  worst <- which(rowMeansCDS < q10)
  worstNames <- cdsTEs[worst, txNames]
  whichNames <- kozakNames %in% worstNames

  meanNumberWorst <- mean(kozak[whichNames])

}

# te of next uORF downstream vs distance to it
teDownstreamVSDistance <- function(){
  uorfTEs <- readTable("teFiltered", with.IDs = T)

  numberOfUorfsPerTx <- readTable("numberOfUorfsPerTx")

  cdsTEs <- readTable("cdsTeFiltered", with.IDs = T)

  #maxTEs <- uorfTEs[, lapply(.SD, max), by = txNames]

  rowMeansUorfs <- rowMeans(uorfTEs[,3:ncol(uorfTEs)])
  q90 <- quantile(rowMeansUorfs, 0.90)
  best <- which(rowMeansUorfs > q90)
  bestNames <- uorfTEs[best, txNames]

  # group by tx
  #

}

distanceVsUorfTE <- function(){
  uorfTEs <- readTable("teFiltered", with.IDs = T)

  dists <- readTable("distORFCDS")

  rowMeansUorfs <- rowMeans(uorfTEs[,3:ncol(uorfTEs)])
  which <- rowMeansUorfs > 1
  plot(dists$distORFCDS[which], rowMeansUorfs[which],
       ylim = c(0.7,100), xlim = c(-1000,1000))
  cor.test(dists$distORFCDS[which], rowMeansUorfs[which])
  numberOfUorfsPerTx <- readTable("numberOfUorfsPerTx")
  whichOne <- numberOfUorfsPerTx$nUorfs == 1
  txToUse <- numberOfUorfsPerTx$txNames[whichOne]

  matchedNames <- uorfTEs$txNames %in% txToUse

  mean(rowMeansUorfs[matchedNames])


  ### for cds te
  cdsTEs <- readTable("cdsTeFiltered", with.IDs = T)
  rowMeansCDS <- rowMeans(cdsTEs[,2:ncol(cdsTEs)])

  mappingTxUorfs <- readTable("linkORFsToTx")$txNames
  cdsDis <- dists$distORFCDS[]
  names(rowMeansCDS) <- cdsTEs$txNames
  rowMeansCDS <- rowMeansCDS[mappingTxUorfs]

  which <- rowMeansCDS > 1
  whichPos <- dists$distORFCDS > 0
  usedDists <- dists$distORFCDS[which]
  group = (usedDists<=0) + (usedDists<=200) + (usedDists <= 400) +(usedDists <= 600)
  groups = factor(group,levels=0:4,labels=c("<0","<200","<400","<600", ">600"))
  boxplot(rowMeansCDS[which]~groups)

  boxplot(dists$distORFCDS[whichPos], rowMeansCDS[whichPos])

  plot(dists$distORFCDS[whichPos], rowMeansCDS[whichPos],
       ylim = c(0.7,100), xlim = c(-1000,1000),
       lines(x = c(-185,-185), y = c(1,100)))


}



predictionVsCageHits <- function(){
  cageTissuesPrediction <- readTable("tissueAtlasByCageAndPred")

  cageTissues <- readTable("tissueAtlasByCage")

  cageRed <- cageTissues[, colSums(cageTissues) > 0, with = F]
  inAll <- sum(rowSums(cageRed) == ncol(cageRed))
  # finalCagePred from table
  bestNames <- names(sort(-colSums(cageRed)))[1:20]
  cageRed <- cageRed[,bestNames, with = F]


  values <- c(colSums(cageRed) - inAll, rep(inAll, ncol(cageRed)))
  variable <- c(rep(colnames(cageRed), 2))
  type <- c(rep("Unique uORFs tissue", ncol(cageRed)), rep("uORFs in all tissues", ncol(cageRed)))
  df <- data.table(value = values, variable, type)
  df <- df[order(value),]
  df$variable <- factor(df$variable, levels = unique(df$variable), ordered = T)
  cageAll <- ggplot(df, aes(x=variable,y=value,fill=type)) +
    geom_bar(stat="identity", position="stack") +
    xlab("Tissue")+ylab("Number of uORFs found by CAGE") +
    theme(axis.text.y = element_text(size = 10)) +
    guides(fill=FALSE) +
    coord_flip() +
    labs(title="Number of uORFs per tissue")

  # for prediction
  inAll <- sum(rowSums(cageTissuesPrediction) == ncol(cageTissuesPrediction))
  bestNames <- names(sort(-colSums(cageTissuesPrediction)))[1:20]

  cageRed <- cageTissuesPrediction[,bestNames, with = F]
  values <- c(colSums(cageRed) - inAll, rep(inAll, ncol(cageRed)))
  variable <- c(rep(colnames(cageRed), 2))
  type <- c(rep("Unique uORFs tissue", ncol(cageRed)), rep("uORFs in all tissues", ncol(cageRed)))
  df <- data.table(value = values, variable, type)
  df <- df[order(value),]
  df$variable <- factor(df$variable, levels = unique(df$variable), ordered = T)
  predAll <- ggplot(df, aes(x=variable,y=value,fill=type)) +
    geom_bar(stat="identity", position="stack") +
    scale_fill_discrete(name="uORF counts in final prediction") +
    xlab("Tissue")+ylab("Number of predicted active uORFs") +
    guides(fill=FALSE) +
    theme(axis.text.y = element_text(size = 10)) +
    coord_flip()

}

varianceTissueUsage <- function(){
  tab <- table(cageTissuesPrediction$ovary, cageTissuesPrediction$brain)
  chi <- chisq.test(tab)
  chi$residuals
  # plan:
  # do pairwise tests, see that things are ok.
  # venn diagram:

  library(VennDiagram)
  grid.newpage()
  boOver <- draw.pairwise.venn(15021, 14213,
                               11523, category = c("Ovary", "Brain"),
                               lty = rep("blank", 2), fill = c("light blue", "yellow"),
                               alpha = rep(0.5, 2), cat.pos = c(0, 0),
                               cat.dist = rep(0.025, 2), title = "abc")
  boOver <- grid.arrange(gTree(children=boOver), top=textGrob("uORF overlaps", gp=gpar(fontsize=20,font=8)),
                           bottom="")

  library(cowplot)
  plot_grid(predAll,boOver, align='hv',nrow=1,labels=c('A','B'))
}


#' How good are the features
#'
#' Te will be used as translation
#' Use cds as positive set and randomized columns of
#' cds as negative.
#'
#' Can a random forrest seperate these two sets ?
correlateFeatures <- function(){
  cdsTEs <- readTable("cdsTeFiltered", with.IDs = F)
  cdsKozak <- readTable("cdsKozak", with.IDs = F)
  cdsORFScore <- readTable("cdsORFScore", with.IDs = F)
  cdsRfpFPKMs <- readTable("cdsRfpFPKMs", with.IDs = F)

  dt <- data.table(cdsTEs[,1], cdsKozak[,1], cdsORFScore[,1], cdsRfpFPKMs[,1])

  # add these to cds
  # ribosomeReleaseScore
  # insideOutsideORF
  # ribosomeStalingScore

  dtRandom <- dt
  for (i in 1:ncol(dt)) {
    ran <- sample(nrow(dtRandom))
    dtRandom[, i] <- dtRandom[ran, i, with = F]
  }

  merged <- rbindlist(list(dt,dtRandom))
  l <- nrow(dtRandom)
  y <- c(rep(1,l),rep(0,l))

  max <- nrow(merged)
  testMerged <- rbindlist(list(dt[1:max,],dtRandom[1:max,]))
  testY <- as.factor(c(rep(1,max),rep(0,max)))


  rf <- randomForest(x = merged, y = y)

}
