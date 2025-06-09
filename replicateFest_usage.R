# 2025-01-24
# Ludmila Danilova, Leslie Cope
# Usage of replicateFest

library(tools)
library(dplyr)
# load functions
devtools::install_github("OncologyQS/replicateFest")
library(replicateFest)


#================================================================
# usage of replicateFest on experimental data with replicates
#================================================================
# load data
# list files with data
files=list.files("Data/",pattern = ".txt", full.name=T, recursive=T)
# extract file names without extension
filenames = file_path_sans_ext(basename(files))
sampAnnot = splitFileName(filenames)

######################
# usage
#====================
# run all clones in a patient and time point and return the results
res = runExperiment(files,
                    peptides = timeSamples$condition,
                    cont= "DMSO",
                    fdrThr = 0.05,
                    xrCond = NA,
                    outputFile = "test.xlsx",
                    saveToFile = T)



#================================
# running analysis for one clone
#================================
# read in data
mergedData=readMergeSave(files, filenames = NULL)$mergedData

# specify a clone
clone1 = "CASSFGRGAEKLFF"
# fit model for the clone
fitModel(clone1,mergedData,
         peptides = sampAnnot$condition,
         control="DMSO")

