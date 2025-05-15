# 2025-01-24
# Ludmila Danilova, Leslie Cope
# Usage of replicateFest on experimental data with replicates
# as an example, using HIV experiment from UCSF, Rutishauser, Rachel group

#These assays were performed to test for HIV-specific responses
#in HIV+ patients who received a vaccine. There are 3 patients,
#an assay performed before and after vaccination in each patient.
#“DMSO” are the negative controls, CEF are the positive controls.
# “HIV08”, “HIV10”, and “HIV12” are known immunodominant epitopes
#derived from the HIV Gag protein, which means clones that are
# positive for the minimal epitopes SHOULD have considerable expansion
#in “Gag” peptide pool conditions as well. However it is also likely
#that many clones that expand to “Gag” will NOT expand in response
#to the minimal epitopes.

#BL is baseline (uncultured) and Env is just another peptide pool condition.
# The reference should be DMSO. Even though it is technically
#the positive control, it would still be nice to see the expansion and
#magnitude of these clones, if that makes sense. So I would just treat
#it as another condition.

library(tools)
library(dplyr)
# load functions
devtools::install_github("OncologyQS/replicateFest")
library(replicateFest)

source("R/manafest_shiny_functions.r")
source("R/repManFunctions.r")

#================================================================
# usage of replicateFest on experimental data with replicates
#================================================================
# load data
# list files with data
files=list.files("Data/",pattern = ".txt", full.name=T, recursive=T)
# parse file names to exptract patient ID, time point, condition, and replicate
# extract file names without extension
filenames = file_path_sans_ext(basename(files))
sampAnnot = splitFileName(filenames)

# create a data frame with sample annotations
sampAnnot = do.call(rbind,strsplit(filenames,"_"))
colnames(sampAnnot)=c("Patient","Time","Condition","Replicate")
sampAnnot = cbind(sampAnnot, path = files)
# fix baseline sample
sampAnnot[grep("BL",sampAnnot[,"Condition"]), "Replicate"]="BL"
sampAnnot = as.data.frame(sampAnnot, stringsAsFactors = F)

# remove baseline samples
sampAnnot = sampAnnot %>% filter(Replicate != "BL")

######################
# usage
#==================================
# run all data iterating through patients and time points
#==================================
#

## iterate unique patients
for(p in unique(sampAnnot$Patient)){
  # get patient samples
  patSamples = sampAnnot %>% filter(Patient ==p)
  # iterate unique time points
  for(t in unique(patSamples$Time)){
    # get time point samples
    timeSamples = patSamples %>% filter(Time == t)
    # run one experiment
    out = runExperiment(timeSamples$path,
                        peptides = timeSamples$Condition,
                        cont= "DMSO", # control condition
                        FDRthr = 0.05, # FDR threshold
                        percentThr = 0.5, # percentage threshold
                        # a list of cross-reactive conditions
                        xrCond = c("Gag","HIV08","HIV10","HIV12"),
                        # output file
                        outputFile = paste0("results/results_xr_FDR_05_perc05/",
                                            p,"_",t,"_xrGag_FDR_05_perc05.xlsx"))
  }
}

#================================
# debugging, running one clone or one patient and time point
#================================
p = "7735"
t = "02242020"
patSamples = sampAnnot %>% filter(Patient ==p)
timeSamples = patSamples %>% filter(Time == t)

# run one clone in a patient and time point
# read in data
mergedData=readMergeSave(timeSamples$path, filenames = NULL)$aaData

# fit model for one clone
clone1 = "CASSFGRGAEKLFF"
# fit model for one clone
fitModel(clone1,mergedData,
         peptides = timeSamples$Condition,control="DMSO",
         c.corr=1)

# run all clones in a patient and time point and return the results
newOut1 = runExperiment(timeSamples$path,
                        peptides = timeSamples$Condition,
                        cont= "DMSO",
                        FDRthr = 0.05,
                        xrCond = c("Gag","HIV08","HIV10","HIV12"),
                        saveToFile = F)



