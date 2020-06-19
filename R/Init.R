#' Init uORFome pipeline
#'
#' Make directory structure for orf finding, create database
#' assign variables and validate input data.
#' @param mainPath folder for uORFome to put results
#' @param df.rfp ORFik experiment of Ribo-seq
#' @param df.rna ORFik experiment of RNA-seq
#' @param df.cage ORFik experiment of CAGE
#' @param dataFolder ORFik experiment that contains Annotation
#' @param makeDatabase default FALSE, set to TRUE if you want
#' to predict translated uORFs
#' @param organism scientific name of organism,
#' like Homo sapiens, Danio rerio, etc.
#' @param biomart default "not_decided" will be automaticly detected by
#' organism name. Set it if you know.
#' @param mode character, default: "uORF". alternative "ORF". Do you want to predict
#' on uORFs or ORFs (CDS etc.)
#' @importFrom tools file_ext
#' @importFrom BiocParallel registered
#' @export
orfikDirs <- function(mainPath, df.rfp, df.rna, df.cage,
                      organism, biomart = "not_decided", mode = "uORF") {
  if (!(mode %in% c("uORF", "ORF")))
    stop("mode must be uORF or ORF")

  dir.create(mainPath, showWarnings = FALSE, recursive = TRUE)
  setwd(mainPath)
  message("Welcome, setting up uORFome folders\n")
  message(p("Registered organism is: ", organism))
  biomart_dataset <- getBiomartFromOrganism(organism)

  message(paste("main path for project will be: ", mainPath))
  # Set directory paths
  if (mode == "uORF") {
    resultsFolder <- p(mainPath, "/results")
    regionUORFsFolder = p(resultsFolder,"/uORF_searchRegion/")
    leadersFolder = p(resultsFolder,"/New_Cage_Leaders/")
    uorfFolder = p(resultsFolder,"/candidate_uORFs/")
    idFolder = p(resultsFolder,"/uorfIDs/")
    plottingFolder = p(resultsFolder,"/Plotting/")
  }

  dataFolder <- p(mainPath, "/helper_files")
  featuresFolder <- p(mainPath, "/features")
  modelFolder <- p(mainPath, "/prediction_model")

  # Create directories
  if (mode == "uORF") {
    dir.create(resultsFolder, showWarnings = FALSE)
    dir.create(leadersFolder, showWarnings = FALSE)
    dir.create(regionUORFsFolder, showWarnings = FALSE)
    dir.create(uorfFolder, showWarnings = FALSE)
    dir.create(idFolder, showWarnings = FALSE)
  }
  dir.create(featuresFolder, showWarnings = FALSE)
  dir.create(dataFolder, showWarnings = FALSE)
  dir.create(modelFolder, showWarnings = FALSE)

  # Create database
  dir.create("dataBase", showWarnings = FALSE)
  dataBaseFolder <- p(mainPath,"/dataBase")
  assign("dataBaseFolder", dataBaseFolder, envir = .GlobalEnv)
  uORFomePipe:::createDataBase(p(dataBaseFolder, "/uorfCatalogue.sqlite"))

  # Assign variables to global environment
  if (mode == "uORF") {
    assign("resultsFolder", resultsFolder, envir = .GlobalEnv)
    assign("regionUORFsFolder", regionUORFsFolder, envir = .GlobalEnv)
    assign("leadersFolder", leadersFolder, envir = .GlobalEnv)
    assign("uorfFolder", uorfFolder, envir = .GlobalEnv)
    assign("idFolder", idFolder, envir = .GlobalEnv)
    assign("plottingFolder", plottingFolder, envir = .GlobalEnv)
  }
  assign("mainFolder", mainPath, envir = .GlobalEnv)
  assign("dataFolder", dataFolder, envir = .GlobalEnv)
  assign("featuresFolder", featuresFolder, envir = .GlobalEnv)
  assign("organism", organism, envir = .GlobalEnv)
  assign("biomart_dataset", biomart_dataset, envir = .GlobalEnv)

  # now validate all that directories exist
  if(mode == "uORF") {
    if (!dir.exists(dataFolder))
      stop(p("Could not find directory: ", c(dataFolder)[!file.exists( dataFolder)]))
  }


  # Set up annotation and save transcript-regions in helper_libraries
  gtfdb = df.rfp@txdb
  txdb <- NULL
  if (!file.exists(p(dataFolder, "/threeUTRs.rds"))) {
    seqstyle <- seqlevelsStyle(ORFik:::findFa(df.rfp@fafile))[1]
    f <- file_ext(gtfdb)
    if (f == "gff" | f == "gff2" | f == "gff3" | f == "gtf") {
      if (!file.exists(p(gtfdb, ".db")) | !file.exists(p(gtfdb, "sqlite"))) {
        txdb <- GenomicFeatures::makeTxDbFromGFF(gtfdb)
      } else {
        if (file.exists(p(gtfdb, ".db"))) {
          gtfdb <- p(gtfdb, ".db")
          txdb <- loadTxdb(gtfdb)
        } else {
          gtfdb <- p(gtfdb, ".sqlite")
          txdb <- loadTxdb(gtfdb)
        }
      }
      txdb <- ORFik:::matchSeqStyle(txdb, seqstyle)
      saveDb(txdb, file = p(dataFolder, "/Gtf.db"))
      txdb <- loadTxdb(p(dataFolder, "/Gtf.db"))

    } else if(f == "db" | f == "sqlite") {
      txdb <- loadTxdb(gtfdb, chrStyle = seqstyle)
      saveDb(txdb, p(dataFolder, "/Gtf.db"))
    } else stop("when txdb is path, must be one of .gff, .gtf and .db")
    loadRegions(txdb, c("tx", "cds", "fiveUTRs", "threeUTRs"))

    saveRDS(tx, file = p(dataFolder, "/tx.rds"))
    saveRDS(cds, file = p(dataFolder, "/cds.rds"))
    saveRDS(fiveUTRs, file = p(dataFolder, "/fiveUTRs.rds"))
    saveRDS(threeUTRs, file = p(dataFolder, "/threeUTRs.rds"))
    if(mode != "uORF") saveRDS(tx, file = p(dataFolder, "/CageFiveUTRs.rds"))

    rm(list = c("tx", "cds", "fiveUTRs", "threeUTRs"), envir = .GlobalEnv)
  }

  assign("faiName", df.rfp@fafile, envir = .GlobalEnv)
  assign("gtfdb", p(dataFolder, "/Gtf.db"), envir = .GlobalEnv)

  validateInputs(df.rfp, df.rna, df.cage)

  message("This is default backend:")
  #print(registered()[1])
  message(p("Using ", registered()[[1]]$workers, " threads from CPU"))
  message(p("Example on how to register default backend to 10 cores:"))
  print("register(MulticoreParam(workers = 10), default=TRUE)")


  message("\nuORFome is now ready")
  return(invisible(NULL))
}

#' Validation of experiments
#' @inheritParams orfikDirs
validateInputs <- function(df.rfp, df.rna, df.cage) {
  samples.rfp <- nrow(df.rfp);samples.rna <- nrow(df.rna)
  if (samples.rfp != samples.rna) stop("Not equal samples in RNA-seq and Ribo-seq")
  if (!all(df.rfp$stage %in% df.rna$stage))
    stop("Not equal tissues/stages in RNA-seq and Ribo-seq")
  if (!all(df.rfp$condition %in% df.rna$condition))
    stop("Not equal conditions in RNA-seq and Ribo-seq")
  if (is.null(df.cage)) {
    message("Running without CAGE")
    message("Cancel if this was wrong, and input correct ORFik experiment!")

  } else if (!all(df.rfp$stage %in% df.cage$stage))
      stop("Not equal tissues/stages in CAGE and Ribo-seq")
  groups <- bamVarName(df.rfp, skip.replicate = TRUE, skip.libtype = TRUE)
  ugroups <- unique(groups)
  if (tableNotExists("experiment_groups")) insertTable(ugroups, "experiment_groups")
  if (tableNotExists("experiment_groups_all")) insertTable(groups, "experiment_groups_all")

  message(p("Tissues validated, will run for ", length(ugroups), " groups:"))
  message("Replicates per group:")
  print(table(groups))
  return(invisible(NULL))
}

# load cds
# remove random size of middle so length resembles ORF, must still have same frame!
# save these CDS and use them as test uORFs
