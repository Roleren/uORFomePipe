#' Get table of cage data
#' @importFrom xlsx read.xlsx
getCageInfoTable <- function(file = "/export/valenfs/projects/uORFome/HumanSamples2.0.sdrf.xlsx"){
  if (tableNotExists("cageInformation")){
    cageTable <- read.xlsx(file, sheetName = "Sheet1")
    if(length(cageFiles) != nrow(cageTable)) warning(paste("Not equal size of cageFiles and cageTable",
                                                           length(cageFiles), "vs", nrow(cageTable)))
    ############################### UPDATE THIS IF USING SELF DEFINED TABLE!!!!!!!################################################
    #filter bad tissues
    # change name to tissue
    colnames(cageTable)[10] <- "tissue"
    # all to lower case
    cageTable$tissue <- tolower(cageTable$tissue)
    # bone <- osteosarcoma
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "osteosarcoma", value = F)] <- "bone"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "omental", value = F)] <- "omentum"
    cageTable$tissue[grep(x = cageTable$Characteristics..Cell.type., pattern = "monocyte", value = F)] <- "blood"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "Hep-2", value = F)] <- "cervix"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "langerhans cell", value = F)] <- "skin"
    cageTable$tissue[ grep(x = cageTable$Charateristics..description., pattern = ", adult", value = F)] <-
      sub("\\,.*", "", grep(x = cageTable$Charateristics..description., pattern = ", adult", value = T))
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "perirenal", value = F)] <- "kidney"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "HES3-GFP Embryonic Stem cells, cardiomyocytic", value = F)] <- "heart"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "mesenchymal stem cells", value = F)] <- "adipose"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "iPS differentiation to neuron", value = F)] <- "brain"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "H1 embryonic stem cells", value = F)] <- "stem cell H1"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "hIPS", value = F)] <- "stem cell iPS"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "melanocytic", value = F)] <- "stem cell H9"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "293SLAM", value = F)] <- "kidney"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "COBL-a rinderpest", value = F)] <- "mesoderm"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "ARPE-19", value = F)] <- "eye"
    cageTable$tissue[grep(x = cageTable$Charateristics..description., pattern = "Renal Glomerular", value = F)] <- "kidney"

    if(any(is.na(cageTable$tissue))) warning("Tissue column contains NAs")

    cageTable <- as.data.table(cageTable)
    cageTable[is.na(tissue)] <- "unclassifiable"
    cageTable[tissue == "adipose tissue"]$tissue <-"adipose"

    #############Fixer region#########################
    #### Good place to debug if you are missing some CAGE libs
    indices <- vector()
    for(i in seq.int(nrow(cageTable))) {
      index <- grep(x = cageFiles, pattern = cageTable$Source.Name[i])
      if (length(index) != 1){
        print(paste("warning for index", i," ",unlist(index)))
      }
      indices <- c(indices, index[1])
    }
    if(length(unique(indices)) != length(indices)) stop("using duplicated CAGE files, check how you group!")
    cageTable$cage_index <- indices

    insertTable(Matrix = cageTable,tableName = "cageInformation", rmOld = TRUE)
  } else{
    cageTable <- readTable("cageInformation")
  }
  return(cageTable)
}

#' get ribo-seq/rna-seq info table
#'
#' If not already made, make it and save it
#' @param name character vector names of .rdata file
getRiboRNAInfoTable <- function(tableName = "/matching_rna_ribo.rdata"){
  load(file = p(mainFolder,tableName))
  return(matching_rna_ribo)
}

#' Fix this
getRiboMatchedToAll <- function() {
  grep(paste(getRiboRNAInfoTable()$ribo, collapse="|"), grep(pattern = "merged",
                                                             x = list.files(rfpFolder), value = T))
}
