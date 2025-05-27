# create a function that returns positive clones for replicateFest

getPositiveClonesReplicates = function(analysisRes,
                                       mergedData,
                                       control,
                                       samp = names(mergedData),
                                       excludeCond = NA,
                                       orThr = 1,
                                       fdrThr = 0.05,
                                       percentThr = 0)
{
  # get all expanded clones relative to the control
  # find colunms with control to get clones expanded relative to reference
  contCol = grep(control,colnames(analysisRes), value = T)
  res_exp = getExpanded(analysisRes[,c("clone",contCol)], mergedData,
                        ORthr = orThr, FDRthr = fdrThr)

  #=================
  # keep clones with maximum percentage across
  # all analyzed samples higher that a specified threshold
  #=================
  # grep columns with percentage
  percCol = grep("percent",colnames(res_exp), value = T)
  #browser()
  # exclude columns with excludeCond
  percCol = percCol[!grepl(paste(excludeCond,collapse = "|"),percCol)]
  # get clones with maximum percentage higher than specified threshold
  res_exp_filtered = res_exp[apply(res_exp[,percCol],1,max) >percentThr,]

  #======
  # find uniquely expanded clones by checking the second best clone
  # save the second best comparison results
  # these are the rest of the columns that are not comparison to reference
  screen_scndBest = analysisRes[rownames(res_exp_filtered),
                               setdiff(colnames(analysisRes),contCol)]
  # check for uniqueness. it should be expanded and significant in comparison to the second best as well
  unique_exp = (screen_scndBest[,grep("OR", colnames(screen_scndBest))]>1 &
                  screen_scndBest[,grep("FDR", colnames(screen_scndBest))]<fdrThr)
  # merge the results
  # add columns that indicates how many and what comparisons were significant
  # results for the second best comparison
  # the rest of info
  sigComCol = c("clone","n_significant_comparisons","significant_comparison")
  res = cbind(res_exp_filtered[unique_exp,sigComCol],
                        res_exp_filtered[unique_exp,setdiff(colnames(res_exp),sigComCol)],
                          screen_scndBest[unique_exp,])
  #=====
  return(res)
}

load("7099_03182029_inputData.rda")

res = fitModelSet(clonesToTest, obj,
                  sampAnnot$condition,
                  excludeCond = input$excludeSamp,
                  control=input$refSamp,
                  c.corr=1)
rownames(res) = res$clone

posCloneRep = getPositiveClonesReplicates(analysisRes,
                                          mergedData,
                                          "DMSO")

# made output of positive clones as a data.frame
posClonesdf = getPositiveClones(fisherRes, mergedData,
                        samp = sampForAnalysis,
                        orThr = as.numeric(input$orThr),
                        fdrThr = as.numeric(input$fdrThr),
                        percentThr = as.numeric(input$percentThr))

tablesToXls = createPosClonesOutput(posClonesdf,
                                    mergedData,
                                    input$refSamp,
                                    NULL,
                                    addDiff = F)

