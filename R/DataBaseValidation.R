#' Plot top tissues
#'
#' Shows max 20 tissues candidate and predicted uORF counts
#' And start and stop codon usage for combined prediction
#' @param saveName full name of location to save to, default
#' "differential_uORF_usage.png"
#' @inheritParams feature.boxplots
#' @inheritParams startAndStopCodonPlots
#' @import ggplot2
#' @import gridExtra
#' @import cowplot
#' @export
predictionVsCageHits <- function(saveName = p("differential_", mode, "_usage.png"),
                                 mode = "uORF",
                                 prediction = predictUorfs(mode = mode),
                                 stop.codons = FALSE, sort.by = "counts") {
  if (!(mode %in% c("uORF", "CDS", "aCDS")))
    stop("mode must be uORF or CDS or aCDS (artificial CDS)")
  message("--------------------------------------")
  message("Creating result plots:")
  pred <- uORFomePipe:::perSampleCountPlot(mode, sort.by)
  codons <- uORFomePipe:::startAndStopCodonPlots(stop.codons)
  features <- uORFomePipe:::feature.boxplots(prediction)
  lay <- rbind(1, 1, 2, 3, 3)
  grid <- cowplot::plot_grid(pred, codons, features, ncol = 1,
                     rel_heights = c(2, 1, 2),
                     labels = c("A", "B", "C"), label_y = c(1, 1.1, 1.1))
  # grid <- gridExtra::grid.arrange(top = paste0(mode," prediction"),
  #                                 pred, codons, features,
  #                                 layout_matrix = lay)
  ggsave(saveName, grid, width = 220, height = 180, units = "mm", dpi = 300)
  # Venn diagram
  venn.diagram.uORFs()
  return(invisible(NULL))
}

#' Per sample library plot overall counts
#'
#' Divide into active and candidate uORFs
#' @param mode character, default: "uORF". alternative "ORF" or what ever
#' region you are checking.
#' @param sort.by character, default "count". Sort rows by counts,
#' or by alphabetical "alpha"
perSampleCountPlot <- function(mode = "uORF", sort.by = "count") {
  themes.no.y <- theme(axis.text.y = element_text(size = 8),
                       panel.grid.major.y = element_blank(),
                       panel.grid.minor.y = element_blank(),
                       panel.grid.minor.x = element_blank())

  if (file.exists(paste0(dataFolder,"/tissueAtlas.rds"))) {
    cageTissues <- readRDS(paste0(dataFolder,"/tissueAtlas.rds"))[,-1]
    cageTissues$union <- rowSums(cageTissues) > 0

    ncolHits <- rowSums(cageTissues)
    inAll <- sum(ncolHits == ncol(cageTissues))
    inMoreThanOneNotAll <- (ncolHits > 2) & (ncolHits != ncol(cageTissues))
    inUnique <- (ncolHits <= 2 & ncolHits > 0) & (ncol(cageTissues) > 2)
    if ((inAll + sum(inMoreThanOneNotAll) + sum(inUnique)) != sum(ncolHits > 0))
      stop("A bug in split of CAGE Prediction hits, report bug on github!")
    top20 <- names(sort(-colSums(cageTissues)))[1:min(20, ncol(cageTissues))]

    cageRed <- cageTissues[, top20, with = FALSE]
    values <- c(rep(inAll, ncol(cageRed)), colSums(cageRed[inMoreThanOneNotAll]),
                colSums(cageRed[inUnique]))
    variable <- rep(colnames(cageRed), 3)
    type <- c(rep("All", ncol(cageRed)), rep("More than one", ncol(cageRed)),
              rep("Unique", ncol(cageRed)))
    df <- data.table(value = values, variable, type)
    df$type <- factor(df$type, levels = c("All", "More than one", "Unique"), ordered = TRUE)
    if (sort.by == "count") {
      df$variable <- factor(df$variable, levels = df[, sum(value),
                                                     by = variable][order(V1),]$variable, ordered = TRUE)
    } else {
      df$variable <- factor(df$variable, levels = rev(c("union",
                                                        as.character(unique(sort(df$variable[df$variable != "union"]))))),
                                                    ordered = TRUE)
    }

    df.cage.levels <- levels(df$variable)
    cageAll <- ggplot(df, aes(x=variable,y=value,fill=type)) +
      geom_bar(stat="identity", position =  position_stack(reverse = TRUE)) +
      scale_fill_manual(values = c("#009E73", "#0072B2", "#D55E00")) +
      xlab("Stage")+ylab("Candidate uORFs by sequence") +
      scale_y_continuous(labels = scales::scientific) +
      theme_bw() +
      themes.no.y +
      guides(fill=FALSE) +
      coord_flip()
    cageAll
  }

  # for prediction
  addit <- ifelse(mode == "CDS", "verify_", "")
  cageTissuesPrediction <- readTable(paste0(addit, "tissueAtlasByCageAndPred"), with.IDs = FALSE)
  colnames(cageTissuesPrediction)[colnames(cageTissuesPrediction) == "total"] <- "union"
  ncolHits <- rowSums(cageTissuesPrediction)
  inAll <- sum(ncolHits == ncol(cageTissuesPrediction))
  inMoreThanOneNotAll <- (ncolHits > 2) & (ncolHits != ncol(cageTissuesPrediction))
  inUnique <- ncolHits <= 2 & ncolHits > 0 & ncol(cageTissuesPrediction) > 2
  if ((inAll + sum(inMoreThanOneNotAll) + sum(inUnique)) != sum(ncolHits > 0))
    stop("A bug in split of Ribo-seq Prediction hits, report bug on github!")
  top20 <- names(sort(-colSums(cageTissuesPrediction)))[1:min(20, ncol(cageTissuesPrediction))]

  cageRed <- cageTissuesPrediction[, top20, with = FALSE]
  values <- c(rep(inAll, ncol(cageRed)), colSums(cageRed[inMoreThanOneNotAll]), colSums(cageRed[inUnique]))
  variable <- rep(colnames(cageRed), 3)
  # type <- c(rep(p(mode,"s in all"), ncol(cageRed)),
  #           rep(p(mode,"s > 1"), ncol(cageRed)), rep(p(mode,"s unique"), ncol(cageRed)))
  type <- c(rep("All", ncol(cageRed)),
            rep("More than one", ncol(cageRed)), rep("Unique", ncol(cageRed)))
  df <- data.table(value = values, variable, type)
  # df$type <- factor(df$type, levels = c(p(mode,"s in all"), p(mode,"s > 1"), p(mode,"s unique")), ordered = TRUE)
  df$type <- factor(df$type, levels = c("All", "More than one", "Unique"), ordered = TRUE)
  if (exists("df.cage.levels")) {
    df$variable <- factor(df$variable, levels = df.cage.levels, ordered = TRUE)
  } else {
    df$variable <- factor(df$variable, levels = df[, sum(value),
                                                   by = variable][order(V1),]$variable, ordered = TRUE)
  }

  predAll <- ggplot(df, aes(x=variable, y=value, fill=type)) +
    geom_bar(stat="identity", position =  position_stack(reverse = TRUE)) +
    scale_fill_manual(values = c("#009E73", "#0072B2", "#D55E00")) +
    scale_y_continuous(labels = scales::scientific) +
    xlab("")+ylab(paste0("Predicted translated ", mode, "s")) +
    guides() +
    theme_bw() +
    themes.no.y +
    coord_flip()
  if (file.exists(paste0(dataFolder,"/tissueAtlas.rds"))) {
    pred <- gridExtra::grid.arrange(cageAll, predAll, ncol = 2)
  } else pred <- predAll
  return(pred)
}

#' Distribution of start and stop codon according to total prediction
#' @param stop.codons = FALSE, if TRUE add a row of stop codon plots
#' @importFrom ggplot2 theme_bw
#' @import gridExtra
#' @return grid object of plot
startAndStopCodonPlots <- function(stop.codons = FALSE) {
  cageTissuesPrediction <- readTable("finalCAGEuORFPrediction")
  uorfData <- getAllSequenceFeaturesTable()

  startAndStop <- data.table(StartCodons = factor(uorfData$StartCodons),
                             StopCodons = factor(uorfData$StopCodons),
                             prediction = cageTissuesPrediction$Matrix == 1)
  x_size <- min(11, (14 / max(8, max(length(levels(startAndStop$StartCodons)), length(levels(startAndStop$StartCodons)))))*5)
  theme <- ggplot2::theme_bw() + theme(axis.text.x = element_text(size = x_size),
                              panel.grid.minor = element_blank())
  candidates <- copy(startAndStop)
  candidates$group <- "candidate"
  translated <- startAndStop[prediction == TRUE,]
  translated$group <- "translated"
  startAndStops <- rbindlist(list(candidates, translated))
  starts <- startAndStops[, (.N), by = .(StartCodons, group)]
  starts[, percent := (V1 / sum(V1))*100, by = group]

  startCandidates <- ggplot(data = starts, aes(x = StartCodons, y = percent, fill = group)) +
    geom_bar(stat="identity", position = "dodge", width = 0.3,) +
    theme + xlab("Start codon") +
    ylab("Percent") + scale_fill_manual(values = c("gray", "#390e70"))

  if (stop.codons) {
    stops <- startAndStops[, (.N), by = .(StopCodons, group)]
    stops[, percent := (V1 / sum(V1))*100, by = group]

    stopCandidates <- ggplot(data = stops, aes(x = StopCodons, y = percent, fill = group)) +
      geom_bar(stat="identity", position = "dodge", width = 0.3,) +
      theme + xlab("Stop codon") +
      ylab("Percent") + scale_fill_manual(values = c("gray", "#390e70"))
    grid <- gridExtra::grid.arrange(startCandidates, stopCandidates,
                                    ncol = 1, top = "Start and stop codon by total prediction")
  } else {
    grid <- gridExtra::grid.arrange(startCandidates, top = "Start codon by total prediction")
  }
  return(grid)
}

#' Box plot of feature difference between active and inactive uORFs
#'
#' @param prediction data.table of prediction values in c(0, 1).
#' @param tissue, character, default readTable("experiment_groups")[[1]][1].
#' @importFrom data.table melt
#' @importFrom ggplot2 ggplot
#' @return ggplot object of boxplots
feature.boxplots <- function(prediction = predictUorfs(mode = mode),
                             tissue = readTable("experiment_groups")[[1]][1]) {
  uorfTable <- makeORFPredictionData(tissue)
  uts <- uorfTable; uts$startCodonPerGroupBest <- NULL; uts$RSS <- NULL; uts$startRegionCoverage <- NULL

  uorfData <- getAllSequenceFeaturesTable()
  uds <- data.table(uorfData, stringsAsFactors = TRUE)
  # Add CDS
  cdsData <- readRDS(file = paste0("features/TrainingData_",tissue,".rds"))[y == 1,]
  cdsData <- cdsData[, which(colnames(cdsData) %in% colnames(uts)), with = FALSE]
  uts <- rbindlist(list(uts, cdsData))

  #if (!is.null(uds$lengths)) uts[,ORF_length := uds$lengths]
  #if (!is.null(uds$distORFCDS)) uts[,distORFCDS := uds$distORFCDS]
  if (!is.null(uts$disengagementScores)) colnames(uts)[which(colnames(uts) == "disengagementScores")] <- "disengagement"
  if (!is.null(uts$startCodonCoverage)) colnames(uts)[which(colnames(uts) == "startCodonCoverage")] <- "startCoverage"
  if (!is.null(uts$startRegionRelative)) colnames(uts)[which(colnames(uts) == "startRegionRelative")] <- "startRelative"
  pred <- c(prediction$predict, rep(2, nrow(cdsData)))


  uts.melt <- data.table::melt(uts)
  uts.melt[, translated := factor(rep(pred, ncol(uts)), levels = c(1, 0, 2), labels = c("yes", "no", "CDS"),
                                  ordered = TRUE)]
  uts.melt[abs(value) < 0.01, value := 0]
  uts.melt[value < -0.1, value := -log2(-value)]
  uts.melt[value > 0.1, value := abs(log2(value))]


  uts.melt.sample <- uts.melt[sample(nrow(uts.melt), size = min(1e6, nrow(uts.melt))),]
  plot <- ggplot(uts.melt.sample, aes(x=variable, y=value)) +
    geom_boxplot(aes(fill = translated), outlier.size = 0.7, outlier.alpha = 0.5) +
    xlab("Feature")+ylab("Score (pseudo-log2)") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    coord_cartesian(ylim = c(-5, 12))
  return(plot)
}

expression.check <- function(prediction = predictUorfs(mode = mode),
                             tissue = readTable("experiment_groups")[[1]][1]) {
  uorfTable <- makeORFPredictionData(tissue)[, fpkmRFP]
  uorfTable <- data.table(txNames = readTable("uORFTxToGene")$txNames, uorfTable)
  uorfTable <- uorfTable[prediction$predict == 1,]
  # Add CDS
  cdsData <- readTable("cdsfpkmRFP")[,unlist(readTable("experiment_groups_all")) %in% tissue, with = FALSE]
  cdsData <- data.table(txNames = names(uORFomePipe:::getCDSFiltered()), cdsData)
  dt <- data.table::merge.data.table(uorfTable, cdsData, by = "txNames", all.x = TRUE)
  colnames(dt) <- c("txNames", "uORF", "CDS")


  plot <- ggplot(dt, aes(x=log2(uORF), y=log2(CDS))) +
    geom_point() +
    xlab("Ribo-seq FPKM log2 (CDS)")+ylab("Ribo-seq FPKM log2 (uORF)") +
    theme_bw()
}

#' venn diagram of uORFs
#'
#' Shows overlap between all stages / tissues
#' @param predictions a data.table of predictions, default: readTable("tissueAtlasByCageAndPred")
#' @param filename 'venn_diagramm_uORFs.png'
#' @param width numeric, default: 5 (unit is inch)
#' @param height numeric, default: 5 (unit is inch)
#' @return invisible(NULL), object saved to disc
#' @importFrom VennDiagram venn.diagram
venn.diagram.uORFs <- function(predictions = readTable("tissueAtlasByCageAndPred"),
                               filename = 'venn_diagramm_uORFs.png',
                               width = 5, height = 5) {
  preds <- copy(predictions)
  preds$total <- NULL
  preds <- preds[rowSums(preds) > 0, ]

  my_list <- list()
  for(i in seq_along(preds)) {
    my_list[[length(my_list) + 1]] <- which(preds[, i, with = F] == 1)
  }
  # Chart
  VennDiagram::venn.diagram(imagetype = "png",
                             x = my_list,
                             category.names = colnames(preds),
                             filename = filename,
                             output=TRUE, units = "in",
                            width = width, height = height)
  return(invisible(NULL))
}


#' Test artificial vs original cds
#'
#' Will see how well the model handles smaller versions of it self
#' @param artificial path to uORFomePipe run of artificial
#' @param original path to uORFomePipe run of original
#' @param output path to save plot
#' @return data.table of hits per length of artificial
#' @importFrom ggpubr ggscatter
#' @export
test.artificial <- function(artificial,
                            output,
                            original = artificial,
                            minimum.group.size = 30) {
  dt <- verification.conf.matrix(artificial, original, minimum.group.size)
  dt2 <- copy(dt)
  dt2 <- dt2[, .(true.positive.rate = sum(true.positive) / (sum(true.positive) + sum(false.negative)),
                 true.negative.rate = sum(true.negative) / (sum(false.positive) + sum(true.negative)),
                 group.size = .N), by = lengths]
  dt2 <- dt2[order(lengths),]
  dt2 <- dt2[group.size >= minimum.group.size,]
  if (!any(is.finite(dt2$true.positive.rate))) {
    warning("No finite group of true positivies to create comparison, returning!")
    return(dt2)
  }
  if (!all(is.finite(dt2$true.positive.rate))) warning("Found verify groups without true positivies, check data!")

  dt2$lengths <- dt2$lengths / 3
  # "TPR & TNR vs length of Artificial CDS"
  gg <- ggpubr::ggscatter(dt2, "lengths", c("true.positive.rate", "true.negative.rate"),
                          numeric.x.axis = TRUE, merge = TRUE,
                          ylab = "rate", xlab = "artificial CDS length (# codons)", add = "loess"); plot(gg)
  ggsave(output, gg, width = 4, height = 3)
  if (sum(dt2$true.positive, na.rm = T) == 0) {
    warning("Not enough true positives (predicted CDS) to make comparison model!")
  } else
    print(cor.test(dt2$lengths, dt2$true.positive))
  return(dt2)
}

#' Get comparison of validation model to full CDS length
#' @inheritParams test.artificial
#' @export
verification.conf.matrix <- function(artificial,
                                     original = artificial,
                                     minimum.group.size = 30) {
  # artificial <- mainPath_aCDS; original = artificial;minimum.group.size = 30
  # test
  uorfdb_full <- dbConnect(RSQLite::SQLite(), p(original,"/dataBase/uorfCatalogue.sqlite"))
  uorfdb_art <-  dbConnect(RSQLite::SQLite(), p(artificial,"/dataBase/uorfCatalogue.sqlite"))
  if (artificial == original) {
    pred_full <- readTable("verify_finalPredWithProb", uorfDB = uorfdb_full)
  } else pred_full <- readTable("finalPredWithProb", uorfDB = uorfdb_full)

  pred_art <- readTable("finalPredWithProb",  uorfDB = uorfdb_art)
  message("Prediction model for full CDS")
  print(summary(pred_full))
  if (sum(pred_full$prediction) == 0) warning("No CDS predicted, check input data!")
  message("Prediction model for truncated CDS (artificial)")
  print(summary(pred_art))
  if (sum(pred_art$prediction) == 0) warning("No artificial CDS predicted, check input data!")

  lengths <- readTable("lengths", uorfDB = uorfdb_art)
  dt <- data.table(pred_full = pred_full$prediction == 1, pred_art = pred_art$prediction == 1, lengths)
  dt[, true.positive := (pred_art & pred_full)]
  dt[, false.positive := (pred_art & !pred_full)]
  dt[, true.negative := (!pred_art & !pred_full)]
  dt[, false.negative := (!pred_art & pred_full)]
  print(summary(dt))
  return(dt)
}

#' Venn diagram between two tissues
#' @param tissue1 logical TRUE / FALSE
#' @param tissue2 logical TRUE / FALSE
#' @param pred default
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

  # library(cowplot)
  # plot_grid(predAll,boOver, align='hv',nrow=1,labels=c('A','B'))
  gridExtra::grid.arrange(predAll, boOver, nrow = 1)
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
#' @param uniqueIDs unique ORF ids from ORFik
#' @param startCodons character, sequence of start codons
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
