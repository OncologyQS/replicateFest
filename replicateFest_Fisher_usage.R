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
                              nReads = as.numeric(input$nReads))


#============================================
# run the analysis without reference
fisherRes = compareWithOtherTopConditions(mergedData,
                                          sampForAnalysis,
                                          nReads = input$nReads,
                                          clones = NULL)

posClones = getPositiveClonesFromTopConditions(fisherRes,
                                         orThr = as.numeric(input$orThr),
                                         fdrThr=as.numeric(input$fdrThr))


#==============
# create and save output
tablesToXls1 = createPosClonesOutput(pos, mergedData,
                                    input$refSamp,
                                    input$baselineSamp, addDiff = F)

WriteXLS('tablesToXls1', 'results_without_ref.xlsx',
         SheetNames = names(tablesToXls), row.names = T)


#========================
# write a function to run the Fisher's test as part of package
#========================
runExperimentFisher=function(files,
                             refSamp,
                             baselineSamp = NULL,
                             ignoreBaseline = TRUE,
                             nCells = 100000,
                             prob = .99,
                             nReads = 500,
                             fdrThr = .05,
                             orThr = 5,
                             percentThr = 0,
                             excludeSamp = '',
                             compareToRef = TRUE,
                             outputFile = "output.xlsx",
                       saveToFile = T)
{
  #### start algorithm read data
  mergedData = readMergeSave(files, filenames = NULL)$mergedData
  # specify samples to analyze
  sampForAnalysis = setdiff(names(mergedData),
                            c(excludeSamp,refSamp,
                              baselineSamp))


  #======================
  # run the analysis
  posClones = NULL
  # with reference
  if(compareToRef) #if there is comparison to the ref sample
  {
    #============================================
    # run the analysis with reference
    #===================
    # select clones to test
    # if the Ignore baseline flag is on,
    #then all clones will be tested
    clonesToTest = NULL
    fisherRes = NULL
    if(!ignoreBaseline)
    {
      baselineFreq = getFreq(clones = names(mergedData[[baselineSamp]]),obj,samp=baselineSamp)
      clonesToTest = rownames(baselineFreq)[which(baselineFreq[,1] > getFreqThreshold(as.numeric(input$nCells),input$prob)*100)] #
    }
    # create comparing pairs (to refSamp)
    compPairs = cbind(sampForAnalysis,rep(input$refSamp,
                                          length(sampForAnalysis)))
    # run pair-wise Fisher's test
    fisherRes = apply(compPairs,1,runFisher,mergedData,
                      clones = clonesToTest,
                      nReadFilter = c(nReads,0))
    # add names of compared conditions
    names(fisherRes) = apply(compPairs,1,paste,collapse = '_vs_')


    # select positive clones with specified thresholds
    posClones = getPositiveClones(fisherRes, mergedData,
                                  samp = sampForAnalysis,
                                  orThr = orThr,
                                  fdrThr=fdrThr,
                                  nReads = nReads)
  }else{ # if there is no comparison to ref sample

    fisherRes = compareWithOtherTopConditions(mergedData,
                                              sampForAnalysis,
                                              nReads = nReads,
                                              clones = NULL)

    # select positive clones with specified thresholds
    posClones = getPositiveClonesFromTopConditions(fisherRes,
                                             orThr = orThr,
                                             fdrThr = fdrThr)
  }


  # create  output
  tablesToXls = createPosClonesOutput(posClones, mergedData,
                                       refSamp,
                                       baselineSamp, addDiff = F)
  #===================
  # add the ref_comparison_only sheet
  #===================
  if(compareToRef) #if there is comparison to the ref sample
  {
    # create a table with results
    resTable = createResTable(fisherRes,mergedData,
                              orThr = orThr,
                              FDR_threshold = fdrThr, saveCI = F)
    if (!is.null(resTable))
    {
      # make a numeric table to save in an Excel spreadsheet
      refCompRes = resTable[,2:ncol(resTable)]
      refCompRes = t(do.call('rbind',refCompRes))
      rownames(refCompRes) = rownames(resTable)

      # table with results of comparison to the reference sample only
      tablesToXls$ref_comparison_only = data.frame(refCompRes,check.names = F)
    }else{
      tablesToXls$ref_comparison_only = data.frame(res = 'There are no significant clones')
    }
  }

  #============
  # add a sheet with parameters
  #============
  s = c(refSamp, baselineSamp,sampForAnalysis)
  # add the baseline threshold percentage and
  # the corresponding number of templates in baseline sample
  baselineThrNames = baselineThrVal = NULL
  if(!input$ignoreBaseline)
  {
    baselineThrNames =c('confidence','nCells',
                        'baseline_threshold_percent',
                        'baseline_threshold_templates')
    freq = getFreqThreshold(as.numeric(input$nCells),input$prob)
    if(input$baselineSamp == 'None')
      baselineSamp = input$refSamp
    baselineThrVal = c(input$prob,input$nCells,freq*100,
                       floor(round(freq*productiveReadCounts[baselineSamp])))
  }
  param = c('Reference_samp','Baseline_sample',
            'Excluded samples','Compare to reference',
            'nTemplates_threshold','FDR_threshold',
            'OR_threshold','Ignore_baseline_threshold',
            'Nucleotide_level', baselineThrNames,
            'nAnalyzedSamples',
            paste(s, 'nTemplates',sep = '_'))
  value = c(toString(refSamp), toString(baselineSamp),
            paste(input$excludeSamp, collapse = ', '),
            input$compareToRef, input$nReads,input$fdrThr,
            input$orThr, input$ignoreBaseline,
            input$nuctleotideFlag, baselineThrVal,
            length(sampForAnalysis), productiveReadCounts[s])

  tablesToXls$parameters = data.frame(param, value)

  # save into Excel file
  if(saveToFile)
  {
    saveResults(tablesToXls, outputFile = outputFile)
  } else return(tablesToXls)


}
