#' Run whole uORFomePipe prediction
#'
#' Steps:\cr
#' 1: Make directory structure for orf finding, create database
#' assign variables and validate input data.\cr
#' 2. Find CAGE transcripts\cr
#' 3. Find uORFs\cr
#' 4. Create database\cr
#' 5. Fill database with NGS and sequence features\cr
#' 6. Train the random forrest model\cr
#' 7. Predict on uORFs\cr
#' 8. Get analysis plots\cr
#' \cr
#' NOTE: IF it crashes it will continue from the point you quit, so delete the mainPath
#' folder if you want fresh rerun.\cr Also do not change working directory after you
#' started running, as this might make the program crash
#' @inheritParams checkAndInitPipe
#' @inheritParams getCandidateuORFs
#' @inheritParams BiocParallel::register
#' @param max.artificial.length integer, default: 100, only applies if mode = "aCDS",
#' so ignore this for most people,
#' when creating artificial ORFs from CDS, how large should maximum ORFs be,
#' this number is 1/6 of maximum size of ORFs (max size 600 if artificialLength is 100)
#' Will sample random size from 6 to that number, if max.artificial.length is
#' 2, you can get artificial ORFs of size (6, 9 or 12) (6, + 6 + (3x1), 6 + (3x2))
#' @param features features to train model on, any of the features created
#' during ORFik::computeFeatures, default:
#' \code{c("countRFP", "disengagementScores", "entropyRFP", "floss",
#' "fpkmRFP","ioScore", "ORFScores", "RRS", "RSS", "startCodonCoverage",
#' "startRegionCoverage","startRegionRelative")}
#' @param requiredActiveCds numeric, default 30. How many CDSs are required to be
#' detected active. Size of minimum positive training set. Will abort if not
#' bigger than this number.
#' @return the prediction as data.table with 3 columns. Prediction (0 or 1),
#' p0 (probability of a negtive prediction), p1 (probability of positive prediction).
#' Only one of p0 and p1 can be > 0.5, and that value will decide if
#' prediction is 0 or 1.
#' @importFrom BiocParallel register
#' @importFrom BiocParallel bpparam
#' @import ORFik
#' @export
#' @examples
#' mainPath <- "~/bio/results/uORFome_Zebrafish"
#' # df.rfp <- read.experiment("path/to/rfp.csv")
#' # df.rna <- read.experiment("path/to/rna.csv") # Not required
#' # df.cage <- read.experiment("path/to/CAGE.csv") # Not required
#' # find_uORFome(mainPath, df.rfp, df.rna, df.cage)
find_uORFome <- function(mainPath, organism = organism(df.rfp), df.rfp, df.rna, df.cage,
                         startCodons = "ATG|CTG|TTG|GTG|AAG|AGG|ACG|ATC|ATA|ATT",
                         stopCodons = "TAA|TAG|TGA", mode = "uORF",
                         requiredActiveCds = 30, max.artificial.length = 100,
                         startCodons.cds.allowed = startCodons,
                         stopCodons.cds.allowed = stopCodons,
                         biomart = "ensembl",
                         features = c("countRFP", "disengagementScores", "entropyRFP", "floss",
                                      "fpkmRFP","ioScore", "ORFScores", "RRS", "RSS", "startCodonCoverage",
                                      "startRegionCoverage","startRegionRelative"),
                         BPPARAM = bpparam()) {
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Create folders, variables and validate input
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  checkAndInitPipe(mainPath = mainPath,
            df.rfp, df.rna, df.cage,
            organism = organism, biomart = biomart,
            mode = mode,
            startCodons.cds.allowed = startCodons.cds.allowed,
            stopCodons.cds.allowed = stopCodons.cds.allowed,
            features = features)
  if (mode == "uORF") {
    #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
    # 2. Find uORF search region per CAGE
    #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
    getLeadersFromCage(df.cage, BPPARAM = BPPARAM)

    #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
    # 3. Find candidate uORFs per CAGE
    #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
    getCandidateuORFs(startCodons = startCodons,
                      stopCodons = stopCodons,
                      BPPARAM = BPPARAM)

    #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
    # 4. make uorf IDs (to get unique identifier per uORF)
    #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
    getIDsFromUorfs(BPPARAM = BPPARAM)

    #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
    # 5. CAGE atlas per tissue and uORF / cage leader objects
    #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
    createCatalogueDB(df.cage)
    #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
    # 6. Find sequence, Ribo-seq and RNA-seq features for training model
    #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  }

  makeTrainingAndPredictionData(df.rfp, df.rna, organism = organism, mode = mode,
                                max.artificial.length = max.artificial.length,
                                features = features,
                                requiredActiveCds = requiredActiveCds,
                                BPPARAM = BPPARAM)

  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # 7. Predict uORFs
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # Do all prediction, the returned object will be the prediction for the "combined" predicted
  # over all stages / tissues
  prediction <- predictUorfs(mode = mode)

  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # 8. Analysis
  #¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤¤#
  # CAGE usage analysis (all tissues and total)
  predictionVsCageHits(prediction = prediction, mode = mode)

  # Feature analysis (Just on first tissue)
  featureAnalysis(prediction, tissue = readTable("experiment_groups")[[1]][1])

  if (mode == "aCDS") { # If artificial, now run for whole to verify.
    find_uORFome(mainPath, organism, df.rfp, df.rna, df.cage,
                 startCodons = startCodons,
                 stopCodons = stopCodons, mode = "CDS",
                 max.artificial.length = max.artificial.length,
                 startCodons.cds.allowed = startCodons.cds.allowed,
                 stopCodons.cds.allowed = stopCodons.cds.allowed)
  } else if (mode == "CDS") {
    dt2 <- test.artificial(artificial = mainPath,
                           output = p(mainPath, "/validation_comparison_",
                                      max.artificial.length, "_rates.png"))
  }

  return(prediction)
}


#' Reassign leaders by CAGE
#'
#' Step 1 of uORFome pipeline: Save leaders and uORF search region
#' (leader + cds)
#' @param cageFiles a ORFik experiment with CAGE files and annotation
#' @inheritParams ORFik::reassignTSSbyCage
#' @inheritParams find_uORFome
#' @return invisible(NULL), files saved to disc
#' @export
getLeadersFromCage <- function(cageFiles, filterValue = 3,
                               BPPARAM = bpparam()) {
  message("Starting to find uORF search spaces")
  uORFomePipe:::getLeaders()
  if (is.null(cageFiles)) {
    message("Running pipeline without CAGE data")
    groups <- readTable("experiment_groups")[[1]]
    if (file.exists(paste0(regionUORFsFolder, "_", groups[length(groups)], ".regionUORF.rds"))) {
      message("finished new 5' UTRs and uORF search regions (already exist)")
      return(invisible(NULL))
    }
    getCDS()
    uORFSearchRegion <- ORFik:::addCdsOnLeaderEnds(fiveUTRs, cds)
    for(g in groups) {
      saveRDS(fiveUTRs, file = paste0(leadersFolder, "_", g, ".leader.rds"))
      saveRDS(uORFSearchRegion, file = paste0(regionUORFsFolder, "_", g, ".regionUORF.rds"))
    }
    return(invisible(NULL))
  }

  if (is(cageFiles, "experiment")) cageFiles <- filepath(cageFiles, "ofst")

  bplapply(cageFiles, FUN = function(cageName, dataFolder, leadersFolder,
                                     regionUORFsFolder, filterValue) {
    exportNamerdata <- paste0(regionUORFsFolder, basename(p(cageName, ".regionUORF.rds")))
    if (file.exists(exportNamerdata)) return(0)
    fiveUTRsCage <- reassignTSSbyCage(fiveUTRs, cageName, filterValue = filterValue)
    exportNamerdataLeader <- paste0(leadersFolder, basename(p(cageName, ".leader.rds")))
    saveRDS(fiveUTRsCage, file = exportNamerdataLeader)

    # Extend cage leaders with CDS
    getCDS()
    uORFSearchRegion <- ORFik:::addCdsOnLeaderEnds(fiveUTRsCage, cds)
    saveRDS(uORFSearchRegion, file = exportNamerdata)
    return(0)
  }, dataFolder = get("dataFolder", envir = .GlobalEnv),
  leadersFolder = get("leadersFolder", envir = .GlobalEnv),
  regionUORFsFolder = get("regionUORFsFolder", envir = .GlobalEnv),
  filterValue = filterValue, BPPARAM = BPPARAM)
  message("finished new 5' UTRs and uORF search regions")
}

#' Find uORFs from new leader regions
#'
#' Step 2 of uORFome pipeline
#' @param folder folder to save, default .GlobalEnv arugment: regionUORFsFolder
#' @param startCodons default "ATG|CTG|TTG|GTG|AAG|AGG|ACG|ATC|ATA|ATT", set to
#' "ATG|CTG|TTG|GTG" for a more certain set.
#' @param stopCodons default "TAA|TAG|TGA"
#' @inheritParams find_uORFome
#' @export
getCandidateuORFs <- function(folder = regionUORFsFolder,
                              startCodons = "ATG|CTG|TTG|GTG|AAG|AGG|ACG|ATC|ATA|ATT",
                              stopCodons = "TAA|TAG|TGA",
                              BPPARAM = bpparam()) {
  message("Searching for candidate uORFs")
  leadersList <- list.files(folder, full.names = TRUE)
  uORFomePipe:::getCDS()
  convertCodonStyle <- function(codons) {
    if (length(codons) > 1) {
      if (all(nchar(codons) == 3)) {
        return(paste(codons, collapse = "|"))
      } else stop("Wrong input of start or stop codons!")
    } else if (length(codons) == 0) stop("Either no start or stop codon given!")
    return(codons)
  }
  startCodons <- convertCodonStyle(startCodons)
  stopCodons <- convertCodonStyle(stopCodons)
  bplapply(leadersList, function(i, cds) {
    saveName = p(uorfFolder, gsub(pattern = "regionUORF.rds", replacement = "uorf.rds",
                                  x = basename(i)))
    if (!file.exists(saveName)) {
      getFasta()
      rangesOfuORFs <- findUORFs(readRDS(i), fa, startCodon = startCodons,
                                 stopCodon = stopCodons,
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
#' @param folder folder to save, default .GlobalEnv arugment: uorfFolder
#' @inheritParams find_uORFome
#' @export
getIDsFromUorfs <- function(folder = uorfFolder, BPPARAM = bpparam()){
  uorfFiles = list.files(folder, full.names = TRUE)

  message("Creating uORF ID's")
  bplapply(uorfFiles, function(i) {
    saveName = paste0(idFolder, gsub("uorf.rds","", basename(i)), "uorfID.rds")
    if (!file.exists(saveName)) {
      saveRDS(unique(ORFik:::orfID(readRDS(i))), file = saveName)
      return(i)
    }
    return(0)
  }, BPPARAM = BPPARAM)
}

#' Main function to fill uORF database
#'
#' Step 4 of uORFome pipeline
#' The rows of tables will always be uORFs in order as the candidate uORF file
#' For transcripts it will be the transcript order in cageTx
#' @inheritParams find_uORFome
#' @export
createCatalogueDB <- function(df.cage,
                              dataBaseFolder = get("dataBaseFolder", envir = .GlobalEnv),
                              idFolder = get("idFolder", envir = .GlobalEnv),
                              dataFolder = get("dataFolder", envir = .GlobalEnv),
                              leadersFolder = get("leadersFolder", envir = .GlobalEnv)) {
  createUniqueIDs(idFolder) # IDs for uORFs as matrix
  createGRObjects(dataFolder, leadersFolder) # GRanges objects for all uORFs
  createUORFAtlas(idFolder, dataFolder) # Per CAGE reassigned tx annotation, does uORF exist ?
  getTissueTable(df.cage, dataFolder)  # Per CAGE tissue, does uORF exist ?
}

#' All features from sequence, Riboseq and RNAseq
#'
#' Step 5 of uORFome pipeline
#' @inheritParams find_uORFome
#' @param biomart character or NULL, default: get("biomart_dataset", envir = .GlobalEnv)
#' @return invisible(NULL)
#' @export
makeTrainingAndPredictionData <- function(df.rfp, df.rna,
                                          organism = get("organism", mode = "character", envir = .GlobalEnv),
                                          biomart = get("biomart_dataset", envir = .GlobalEnv),
                                          mode = "uORF",
                                          features = c("countRFP", "disengagementScores", "entropyRFP", "floss",
                                                       "fpkmRFP","ioScore", "ORFScores", "RRS", "RSS", "startCodonCoverage",
                                                       "startRegionCoverage","startRegionRelative"),
                                          max.artificial.length, requiredActiveCds = 30,
                                          BPPARAM = bpparam()) {
  message("--------------------------------------")
  if (!(mode %in% c("uORF", "CDS", "aCDS")))
    stop("mode must be uORF or CDS or aCDS (artificial CDS)")

  ## TRAINING data: Ribo-seq features for cds and 3'
  message("Making training data from CDS and trailers (3' UTRs):")
  getGeneralRiboFeatures(df.rfp, grl = getCDSFiltered(), preName = "cds", BPPARAM = BPPARAM)
  getThreeUTRs()
  getGeneralRiboFeatures(df.rfp, grl = threeUTRs[widthPerGroup(threeUTRs) > 5],
                         preName = "three", threeUTRsSpecial = getSpecialThreeUTRs(),
                         BPPARAM = BPPARAM)
  uORFomePipe:::makeTrainingData(df.rfp, max.artificial.length = max.artificial.length,
                                 mode = mode, features = features,
                                 requiredActiveCds = requiredActiveCds)

  ## Prediction data: Sequence and Ribo-seq features for uORFs / artificial CDS
  # first sequence features
  message("--------------------------------------")
  orfs <- uORFomePipe:::getUorfsInDb(mode = mode)
  getSequenceFeatures(organism, biomart, orfs, mode = mode)
  # Ribo-seq features for ORFs

  getGeneralRiboFeatures(df.rfp, df.rna, orfs,
                         preName = ifelse(mode == "CDS", "verify", ""),
                         BPPARAM = BPPARAM)
  uORFomePipe:::makeORFPredictionData(df.rfp, mode = mode, features = features)

  message("Training complete")
  return(invisible(NULL))
}
