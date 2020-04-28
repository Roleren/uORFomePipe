#' method can also be binary
clusterHelper <- function(table, saveLocation, method = NULL, threads = 25,
                          saveDistVectorName = NULL, mainText = NULL){
  library(parallelDist)
  if (!is.null(method)) {
    binDist <- parallelDist(x = t(table), method = method, threads = threads)
  } else {
    binDist <- parallelDist(x = t(table), threads = threads)
  }
  if (!is.null(saveDistVectorName)) {
    save(binDist, file = saveDistVectorName)
  } 

  library(fastcluster)
  hclustResult <-  hclust(binDist)
  pdfPath <- p(dataBaseFolder, saveLocation)
  if (ncol(table) > 50) {
    if (ncol(table) > 100) {
      if (ncol(table) > 1000) {
        if(ncol(table) > 1500) {
          pdf(pdfPath, width = 55, height = 32)
        } else {
          pdf(pdfPath, width = 40, height = 24)
        }
      } else {
        pdf(pdfPath, width = 25, height = 16)
      }
    } else {
      pdf(pdfPath, width = 20, height = 12)
    }
  } else {
    pdf(pdfPath)
  }
  if (!is.null(mainText)) {
    plot(hclustResult, main = "Clustering of cage experiments")
  } else {
    plot(hclustResult)
  }
  
  dev.off()
  print(paste("ok clustering", saveLocation))
  return(NULL)
}

clusterUorfFeature <- function(uorfFeature, indices = NULL, saveLocation, folder = "/export/valenfs/projects/uORFome/dataBase/clustering/"){
  # cluster
  if (is.null(indices)) {
    uorfFeature <- t(uorfFeature) # transpose for tissue
  } else {
    uorfFeature <- t(uorfFeature[indices]) # transpose for tissue
  }
  mainName <- sub(x = saveLocation, pattern = "*.pdf", replacement = " clustering")
  clusterHelper(uorfFeature, saveLocation =  paste0(folder, saveLocation), mainText = mainName)
  pdfPath <- paste0(folder, saveLocation)
  library(fastcluster)
  clustering <- hclust(dists)
  if (ncol(uorfFeature) > 50) {
    pdf(pdfPath, width = 20, height = 12)
  } else {
    pdf(pdfPath)
  }
  
  plot(clustering)
  dev.off()
}

clusterAllUorfFeatures <- function(){
  # first filter out non-feature tables
  list <- listTables()
  list <- list[-grep(x = list, pattern = "cds", ignore.case = F)]
  list <- list[-grep(x = list, pattern = "uorf", ignore.case = T)]
  list <- list[-grep(x = list, pattern = "information", ignore.case = T)]
  list <- list[-grep(x = list, pattern = "Variance", ignore.case = F)]
  list <- list[-grep(x = list, pattern = "link", ignore.case = F)]
  list <- list[-grep(x = list, pattern = "unique", ignore.case = F)]
  list <- list[-grep(x = list, pattern = "Info", ignore.case = F)]
  list <- list[-grep(x = list, pattern = "Tissue", ignore.case = T)]
  list <- list[-grep(x = list, pattern = "teUnfiltered", ignore.case = T)]
  
  riboInfo <- readTable("RiboSeqInfo")
  rnaInfo <- readTable("RNASeqInfo")
  for (feature in list) {
    tab <- readTable(feature, with.IDs = F)
    if(ncol(tab) == 103) { #rfp seq
      colnames(tab) <- riboInfo$Tissue.Cell_line
    } else if (ncol(tab) == 43) {
      colnames(tab) <- rnaInfo$Tissue.Cell_line
    } else if(ncol(tab) == 35) {
      indices <- as.integer(gsub(pattern = "fpkmRFP_", replacement = "", x = colnames(tab)))
      colnames(tab) <- riboInfo$Tissue.Cell_line[indices]
    } else {
      if (ncol(tab) != 1) {
        stop(paste("could not find ncol for feature: ", feature))
      } 
    }
    
    if (ncol(tab) != 1) {
      clusterUorfFeature(uorfFeature = tab, saveLocation = paste0(feature, ".pdf"))
    } 
    
  }
}

clusterCDSTEs <- function(){
  inputDT <- readTable("cdsTeFiltered", with.IDs = T)
  pattern <- "result."
  inputDTNon <- removeIDColumns(inputDT)
  indices <- as.integer(gsub(pattern = pattern, replacement = "", x = colnames(inputDTNon)))
  if(length(indices) == 0) stop("could not find te indices from colExclusion")
  tissues <- info$Tissue.Cell_line[indices] 
  
  colnames(inputDTNon) <- tissues
  clusterUorfFeature(inputDTNon,
                     saveLocation = "cdsTEsCluster.pdf")
  
  #jaccard index
}

clusterByCageTissue <- function(){
  load(p(dataBaseFolder, "/tissueAtlas.rdata"))
  tissueAtlas$uorfID <- NULL
  binDist <- dist(x = t(tissueAtlas), method = "binary")
  
  library(fastcluster)
  hclustResult <-  hclust(binDist)
  pdfPath <- p(dataBaseFolder, "/clustering/cageByTissue.pdf")
  if (ncol(tissueAtlas) > 50) {
    if (ncol(tissueAtlas) > 100) {
      pdf(pdfPath, width = 25, height = 16)
    } else {
      pdf(pdfPath, width = 20, height = 12)
    }
  } else {
    pdf(pdfPath)
  }
  
  plot(hclustResult, main = "Clustering of cage tissues")
  dev.off()
}

clusterByCageExperiments <- function() {
  load(p(dataBaseFolder, "/UORFAtlas.rdata"))
  uorfAtlas[, uorfID := NULL]
  cageTable <- getCageInfoTable()
  
  cageWeHave <- getAllUsableCage(cageTable)
  uorfAtlasReord <- uorfAtlas
  
  uorfAtlasReord <- uorfAtlasReord[, cageWeHave$cage_index, with = F]
  colnames(uorfAtlasReord) <- cageWeHave$Characteristics.Tissue.
  rm(cageTable)
  rm(uorfAtlas)
  library(parallelDist)
  binDist <- parallelDist(x = t(uorfAtlasReord), method = "binary", threads = 70)
  save(binDist, file = "binDist.rdata")
  library(fastcluster)
  hclustResult <-  hclust(binDist)
  pdfPath <- p(dataBaseFolder, "/clustering/cageByExperiments.pdf")
  if (ncol(uorfAtlasReord) > 50) {
    if (ncol(uorfAtlasReord) > 100) {
      if (ncol(uorfAtlasReord) > 1000) {
        pdf(pdfPath, width = 40, height = 24)
      } else {
        pdf(pdfPath, width = 25, height = 16)
      }
    } else {
      pdf(pdfPath, width = 20, height = 12)
    }
  } else {
    pdf(pdfPath)
  }
  library(dendextend)
  dend <- as.dendrogram(hclustResult)
  colors <- as.numeric(as.factor(colnames(uorfAtlasReord)))
  labels_colors(dend) <- colors
  plot(dend, main = "Clustering of cage experiments")
  rect.hclust(hclustResult, k = 30, border = "red")
  dev.off()
  rm(binDist)
}

#' Ribo seq clustering of 4 groups
clustering <- function(){
  
  # gonzales 1:5, kidney, 18:25, remember +1 for uorfid!
  riboAtlas <- readTable("riboAll")
  
  usedSamples <- riboAtlas[, .(V2,V3,V4,V5,V6, V19,V20,V21,V22,V23,
                               V49,V50,V51,V52,V53, V7,V8,V9,V10,V11 )]
  names1 = paste("Brain",1:5)
  names2 = paste("Kidney",1:5)
  names3 = paste("Blood",1:5)
  names4 = paste("FibroBlast",1:5)
  colnames(usedSamples)[1:5] = names1
  colnames(usedSamples)[6:10] = names2
  colnames(usedSamples)[11:15] = names3
  colnames(usedSamples)[16:20] = names4
  binDist <- dist(x = t(usedSamples), method = "euclidean")
  hclustResult <-  hclust(binDist)
  plot(hclustResult)
  
  # cage seq of clustering of 4 groups
  load("./UORFAtlas.rdata")
  #blood: 284-288, brain: 562-566, kidney: 857-861, fibro(gum): 811,812,814, 815,817
  usedCage <- uorfAtlas[, .(`562`,`563`,`564`,`565`,`566` ,`857`,`858`,`859`,`860`,`861`
                            ,`284`,`285`,`286`,`287`,`288` ,`811`,`812`,`814`,`815`,`817`)]
  colnames(usedCage)[1:5] = names1
  colnames(usedCage)[6:10] = names2
  colnames(usedCage)[11:15] = names3
  colnames(usedCage)[16:20] = names4
  binDist <- dist(x = t(usedCage), method = "binary")
  hclustResult <-  hclust(binDist)
  plot(hclustResult)
}