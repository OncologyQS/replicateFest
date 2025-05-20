# 2025-02-06
# Ludmila Danilova, Leslie Cope

# This script demonstrates the usage of the replicateFest package
# and the runExperimentFisher function on experimental data without replicates

# load package
devtools::install_github("OncologyQS/replicateFest")
library(replicateFest)
library(WriteXLS)

source("R/manafest_shiny_functions.r")

#================================================================
# usage of replicateFest on experimental data without replicates
#================================================================
# load previously saved data
load("no_replicates_inputData.rda")

# set parameters for the analysis
input = list()
input$excludeSamp = ''
input$refSamp = "FW1525_NoPep"
input$baselineSamp = NULL

input$nReads = 50
input$fdrThr = .05
input$orThr = 5
input$percentThr = 0

# specify samples to analyze
sampForAnalysis = setdiff(names(mergedData),
                          c(input$excludeSamp,input$refSamp,
                            input$baselineSamp))
#============================================
# run the analysis with reference
# create comparing pairs (to refSamp)
compPairs = cbind(sampForAnalysis,rep(input$refSamp,
                                      length(sampForAnalysis)))
# run pair-wise Fisher's test
fisherRes = apply(compPairs,1,runFisher,mergedData,
                  clones = NULL,
                  nReadFilter = c(as.numeric(input$nReads),0))
# add names of compared conditions
names(fisherRes) = apply(compPairs,1,paste,collapse = '_vs_')

# select positive clones with specified thresholds
posClones = getPositiveClones(fisherRes, mergedData, samp = sampForAnalysis,
                              orThr = as.numeric(input$orThr),
                              fdrThr=as.numeric(input$fdrThr),
                              percentThr = as.numeric(input$percentThr))


#============================================
# run the analysis without reference
fisherRes = compareWithOtherTopConditions(mergedData,
                                          sampForAnalysis,
                                          nReads = input$nReads,
                                          clones = NULL)

posClones = getPositiveClonesFromTopConditions(fisherRes,
                                               orThr = as.numeric(input$orThr),
                                         fdrThr=as.numeric(input$fdrThr),
                                         percentThr = as.numeric(input$percentThr),
                                         mergedData, samp = sampForAnalysis)


#==============
# create and save output
tablesToXls1 = createPosClonesOutput(posClones, mergedData,
                                    input$refSamp,
                                    input$baselineSamp, addDiff = F)

WriteXLS('tablesToXls1', 'results_without_ref.xlsx',
         SheetNames = names(tablesToXls), row.names = T)


#========================
# write a function to run the Fisher's test as part of package
#========================
source("R/manafest_shiny_functions.r")
#========================
# function debug
inputDir = "C:/Users/Luda/OneDrive - Johns Hopkins/JHU/Manafest/20241009_v10_1_debug_for_Eleni/Data"
files = list.files(inputDir, full.names = T, pattern = "txt")
exSamp = c("FW1525_HIV","FW1525_FFPE_cancer_3417A",
           "FW1525_CEF","FW1525_Baseline_1",
           "FW1525_RL2807_FFPE_skin")
runExperimentFisher(files[15:27],
                    refSamp = "FW1525_NoPep",
                    baselineSamp = NULL,
                    ignoreBaseline = TRUE,
                    nCells = 100000,
                    prob = .99,
                    nReads = 50,
                    fdrThr = .05,
                    orThr = 5,
                    percentThr = 0,
                    excludeSamp = exSamp,
                    compareToRef = TRUE,
                    outputFile = "output.xlsx",
                    saveToFile = T)

