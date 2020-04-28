#' Reassign leaders by CAGE
#'
#' Step 1 of uORFome pipeline
#' @param cageFiles a ORFik experiment with CAGE files and annotation
#' @return uORFs search region (CAGE leaders + CDS)
#' @export
getLeadersFromCage <- function(cageFiles, BPPARAM = bpparam()){
  message("Starting to find uORF search spaces")
  if (is(cageFiles, "experiment")) cageFiles <- filepath(cageFiles, "bedo")
  uORFomePipe:::getLeaders()
  bplapply(cageFiles, FUN = function(cageName, dataFolder, leadersFolder, regionUORFsFolder) {
    fiveUTRsCage <- reassignTSSbyCage(fiveUTRs, cageName, filterValue = 3)
    exportNamerdataLeader = paste0(leadersFolder, basename(p(cageName, ".leader.rds")))
    saveRDS(fiveUTRsCage, file = exportNamerdataLeader)

    # Extend cage leaders with CDS
    getCDS()
    uORFSearchRegion <- ORFik:::addCdsOnLeaderEnds(fiveUTRsCage, cds)
    exportNamerdata = paste0(regionUORFsFolder, basename(p(cageName, ".regionUORF.rds")))
    saveRDS(uORFSearchRegion, file = exportNamerdata)

  }, dataFolder = get("dataFolder", envir = .GlobalEnv),
  leadersFolder = get("leadersFolder", envir = .GlobalEnv),
  regionUORFsFolder = get("regionUORFsFolder", envir = .GlobalEnv), BPPARAM = BPPARAM)
  message("finished new 5' UTRs and uORF search regions")
}

#' Find uORFs from new leader regions
#'
#' Step 2 of uORFome pipeline
#' @param startCodons default "ATG|CTG|TTG|GTG|AAG|AGG|ACG|ATC|ATA|ATT"
#' @export
getCandidateuORFs <- function(folder = regionUORFsFolder,
                              startCodons = "ATG|CTG|TTG|GTG|AAG|AGG|ACG|ATC|ATA|ATT",
                              BPPARAM = bpparam()) {
  message("Searching for candidate uORFs")
  leadersList = list.files(folder, full.names = TRUE)
  uORFomePipe:::getCDS()
  bplapply(leadersList, function(i, cds) {
    saveName = p(uorfFolder, gsub(pattern = "regionUORF.rds", replacement = "uorf.rds",
                                  x = basename(i)))
    if (!file.exists(saveName)) {
      uORFomePipe:::getFasta()
      rangesOfuORFs <- findUORFs(readRDS(i), fa, startCodon = startCodons,
                                 minimumLength = 0, longestORF = FALSE)
      rangesOfuORFs <- ORFik:::filterUORFs(rangesOfuORFs, get("cds", mode = "S4"))
      saveRDS(rangesOfuORFs, file = saveName)
      return(i)
    }
    return(0)
  }, cds = get("cds", mode = "S4"), BPPARAM = BPPARAM)
}

#' Find unique uORF ID's from uORFs
#'
#' Step 3 of uORFome pipeline
#' @export
getIDsFromUorfs <- function(folder = uorfFolder, BPPARAM = bpparam()){
  uorfFiles = list.files(folder, full.names = TRUE)

  message("Creating uORF ID's")
  bplapply(uorfFiles, function(i) {
    saveName = paste0(idFolder, gsub("uorf.rds","", basename(i)), "uorfID.rds")
    saveRDS(unique(ORFik:::orfID(readRDS(i))), file = saveName)
    print("ok")
  }, BPPARAM = BPPARAM)
}

#' Main function to fill uORF database
#'
#' Step 4 of uORFome pipeline
#' The rows of tables will always be uORFs in order as the candidate uORF file
#' For transcripts it will be the transcript order in cageTx
#' @export
createCatalogueDB <- function(df.cage,
                              dataBaseFolder = get("dataBaseFolder", envir = .GlobalEnv),
                              idFolder = get("idFolder", envir = .GlobalEnv),
                              dataFolder = get("dataFolder", envir = .GlobalEnv),
                              leadersFolder = get("leadersFolder", envir = .GlobalEnv)) {
  uORFomePipe:::createUniqueIDs(idFolder) # IDs for uORFs as matrix
  uORFomePipe:::createGRObjects(dataFolder, leadersFolder) # GRanges objects for all uORFs
  uORFomePipe:::createUORFAtlas(idFolder, dataFolder) # Per CAGE reassigned tx annotation, does uORF exist ?
  uORFomePipe:::getTissueTable(df.cage, dataFolder)  # Per CAGE tissue, does uORF exist ?
}

#' All features from sequence, Riboseq and RNAseq
#'
#' Step 5 of uORFome pipeline
#' @export
makeTrainingAndPredictionData <- function(df.rfp, df.rna, tissue, organism, biomart) {
  # first sequence features
  getSequenceFeatures(organism, biomart)
  # Ribo-seq features for ORFs
  getGeneralRiboFeatures(df.rfp, df.rna,
                         grl = uORFomePipe:::getUorfsInDb())
  # Ribo-seq features for cds and 3'
  getCDS()
  getGeneralRiboFeatures(df.rfp, grl = cds[widthPerGroup(cds) > 5], preName = "cds")
  getThreeUTRs()
  getGeneralRiboFeatures(df.rfp, grl = threeUTRs[widthPerGroup(threeUTRs) > 5],
                         preName = "three", threeUTRsSpecial = getSpecialThreeUTRs())

  uORFomePipe:::makeTrainingData(tissue)
  uORFomePipe:::makeORFPredictionData(tissue)
  message("Training complete")
  return(invisible(NULL))
}
