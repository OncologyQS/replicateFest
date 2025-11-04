# 2025-02-06
# Ludmila Danilova, Leslie Cope

# This script demonstrates the usage of the replicateFest package
# and the runExperimentFisher function on experimental data without replicates

# load package
devtools::install_github("OncologyQS/replicateFest")
library(replicateFest)


#================================================================
# usage of replicateFest on experimental data without replicates
#================================================================#========================
# specify a folder with an input data
inputDir = "./tests/testthat/testdata/no_replicates/"

# list paths to files with data
files = list.files(inputDir, full.names = T,
                   pattern = "tsv", recursive = TRUE)
# run analysis and save results
res = runExperimentFisher(files,
                          refSamp = "sample1_control",
                          nReads = 10,
                          fdrThr = .05,
                          orThr = 1,
                          percentThr = 0,
                          condThr = 20,
                          excludeSamp = "",
                          compareToRef = TRUE,
                          outputFile = "testdata-output.xlsx",
                          saveToFile = T)
#==================
# when there is no comparison to ref
res = runExperimentFisher(files,
                          refSamp = NULL,
                          nReads = 10,
                          fdrThr = .05,
                          orThr = 1,
                          percentThr = 0,
                          condThr = 0,
                          excludeSamp = "sample1_control",
                          compareToRef = FALSE,
                          outputFile = "testdata-output_noRef.xlsx",
                          saveToFile = T)


#==============================
# step by step analysis of previously saved data
#==============================
# load previously saved data
load("no_replicates_inputData.rda")

# set parameters for the analysis
input = list()
input$excludeSamp = ''
input$refSamp = "FW1525_NoPep"

input$nReads = 50
input$fdrThr = .05
input$orThr = 5
input$percentThr = 0

# specify samples to analyze
sampForAnalysis = setdiff(names(mergedData),
                          c(input$excludeSamp,input$refSamp))
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
tablesToXls = createPosClonesOutput(posClones, mergedData,
                                    input$refSamp)

saveResults(tablesToXls,"analysisResults.xlsx")
