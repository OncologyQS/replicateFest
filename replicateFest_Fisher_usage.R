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


