# 2025-02-06
# Ludmila Danilova, Leslie Cope

# load package
devtools::install_github("OncologyQS/replicateFest")
library(replicateFest)
library(WriteXLS)

#================================================================
# usage of replicateFest on experimental data without replicates
#================================================================
# load previously saved data
load("C:/Users/Luda/OneDrive - Johns Hopkins/JHU/Manafest/20200429_v13/Ny016-014Input.rda")

# set parameters for the analysis
input = list()
input$excludeSamp = ''
input$refSamp = "NY016-014_No_Pep"
input$baselineSamp = NULL
input$ignoreBaseline = TRUE

input$nCells = 100000
input$prob = .99
input$nReads = 500
input$fdrThr = .05
input$orThr = 5

#============================================
# run the analysis with reference
# specify samples to analyze
sampForAnalysis = setdiff(names(mergedData),
                          c(input$excludeSamp,input$refSamp,
                            input$baselineSamp))
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
                              nReads = as.numeric(input$nReads))
# create all output tables
tablesToXls = createPosClonesOutput(posClones, mergedData,
                                    input$refSamp,
                                    input$baselineSamp,
                                    addDiff = F)
# ad write them into Excel file
WriteXLS('tablesToXls', 'test_v13.xlsx',
         SheetNames = names(tablesToXls), row.names = T)

#============================================
# run the analysis without reference
fisherRes = compareWithOtherTopConditions(mergedData, productiveReadCounts, sampForAnalysis, nReads = 10, clones = NULL)

pos = getPositiveClonesFromTopConditions(fisherRes, orThr = 5, fdrThr = 0.05)

tablesToXls = createPosClonesOutput(pos, mergedData,
                                    productiveReadCounts,
                                    input$refSamp,
                                    input$baselineSamp, addDiff = F)

WriteXLS('tablesToXls', 'test_v13.xlsx',
         SheetNames = names(tablesToXls), row.names = T)

