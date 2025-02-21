# 2025-02-06
# Ludmila Danilova, Leslie Cope

# load package
devtools::install_github("OncologyQS/replicateFest")
library(replicateFest)

#================================================================
# usage of replicateFest on experimental data without replicates
#================================================================
# load previously saved data
load("C:/Users/Luda/OneDrive - Johns Hopkins/JHU/Manafest/20200429_v13/Ny016-014Input.rda")

# set parameters for the analysis
excludeSamp = ''
refSamp = "NY016-014_No_Pep"
baselineSamp = NULL
ignoreBaseline = TRUE

sampForAnalysis = setdiff(names(mergedData), c(excludeSamp,refSamp, baselineSamp))
nCells = 100000
prob = .99
nReads = 500
fdrThr = .05
orThr = 5

#============================================
# run the analysis with reference
# create comparing pairs (to refSamp)
compPairs = cbind(sampForAnalysis,rep(refSamp,
                                      length(sampForAnalysis)))
fisherRes = apply(compPairs,1,runFisher,mergedData,
                  clones = clonesToTest,
                  nReadFilter = c(as.numeric(nReads),0))
names(fisherRes) = apply(compPairs,1,paste,collapse = '_vs_')

posClones = getPositiveClones(fisherRes, mergedData, samp = sampForAnalysis,
                              orThr = as.numeric(input$orThr), fdrThr=as.numeric(input$fdrThr), nReads = as.numeric(input$nReads))

#============================================
# run the analysis without reference
fisherRes = compareWithOtherTopConditions(mergedData, productiveReadCounts, sampForAnalysis, nReads = 10, clones = NULL)

pos = getPositiveClonesFromTopConditions(fisherRes, orThr = 5, fdrThr = 0.05)

tablesToXls = createPosClonesOutput(pos, mergedData, productiveReadCounts, refSamp, baselineSamp, addDiff = F)

WriteXLS('tablesToXls', 'test_v13.xlsx', SheetNames = names(tablesToXls), row.names = T)

