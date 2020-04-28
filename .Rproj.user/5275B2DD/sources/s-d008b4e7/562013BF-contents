#########Variables used by uorfome data-base############
#Change them if needed on other folders!
#####SLASH IS ALWAYS ADDED IN START#####

# set these 4 directories

mainFolder = "/export/valenfs/projects/uORFome" # main folder is one back from RCode1/ folder

codeFolder = p(mainFolder,"/RCode1") # the rcode location

resultsFolder = p(mainFolder,"/results") #output folder

dataFolder = p(mainFolder,"/Annotations") #location of created gtf, fasta and .fai


# now validate all that directories exist
if(!all(dir.exists(c(codeFolder, resultsFolder, dataFolder)))){
  stop(p("Could not find directory: ", c(codeFolder, resultsFolder, dataFolder)[!file.exists(c(codeFolder, resultsFolder, dataFolder))]))
}

# input folders (cage, ribo and rna) Optional (set to NULL if not needed)
cageFolder = "/export/valenfs/projects/uORFome/DATA/CAGE/human/"

faiName = p(dataFolder,"/Homo_sapiens.GRCh38.dna.primary_assembly.chr.fa")
gtfName = p(dataFolder,"/Homo_sapiens.GRCh38.79.chr.NO_PATCH.gtf")
gtfdb = p(dataFolder,"/Gtf.db")  ### a speed up for Gtf, remove if not used
# now validate all files exist
if(!all(file.exists(c(faiName, gtfName, gtfdb)))){
  stop(p("Could not find file: ", c(faiName, gtfName, gtfdb)[!file.exists(c(faiName, gtfName, gtfdb))]))
}

# if(!is.null(cageFolder)){
#   if(!dir.exists(cageFolder)){
#     stop("cage folder not found")
#   }
# }
#
# if(!is.null(rfpFolder)){
#   if(!dir.exists(cageFolder)){
#     stop("ribo-seq folder not found")
#   }
# }
#
# if(!is.null(rfpFolder)){
#   if(!dir.exists(cageFolder)){
#     stop("rna-seq folder not found")
#   }
# }
