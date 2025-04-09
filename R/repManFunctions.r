#' @import tools
#' @import immunarch
#' @import lme4
#' @import contrast
#' @import multcomp
#' @import MASS
#' @import dplyr
#' @import openxlsx
#' @import stats
#'
#import description end
0


# Functions for analysis of FEST data in shiny app
# using input with replicates
# based on Leslie's development for the Cervical SPORE project


#' readMergeSave
#' Reads files, removes non-productive sequences, extracts counts,
#' creates all necessary objects for further analysis
#' @export
#' @param files a list of filenames with full paths
#' @param filenames a vector of filenames
#' @return a list of AA level counts, nucleotide level counts, and total number of reads
#' @examples
#' readMergeSave(files = c("data/0.csv","data/1.csv","data/2.csv"))
#'
# input: a path to folder with input files
# output:
# countData is AA level counts (merged by AA sequences)
# ntData is nucleotide level counts
# productiveReadCounts is the total number of reads of productive sequencies
readMergeSave = function(files, filenames = NULL)
{
	if (length(files) == 0)
	{
		print('There is no files to read')
			return(NULL);
	}
		require(tools)
		require(immunarch)

  # output objects
		mergedData = ntData = list()

		# read all files with immunarch functionality
		repertoire = repLoad(.path = files)

		# the names of files that were actually read in
		readFiles = c()
		for(i in names(repertoire$data))
		{
#			print(i)
		  dat = repertoire$data[[i]]
				#count reads of productive sequences only
				mergedData[[i]] = tapply(dat$Clones, dat$CDR3.aa, sum, na.rm = T)
				# nucleotide level data
				ntData[[i]] =tapply(dat$Clones, dat$CDR3.aa, sum, na.rm = T)
				readFiles = c(readFiles, i)
		}
		if (length(mergedData) == 0)
		{
			print(paste('There are no data to read'))
			return(NULL)
		}
		# if file names are not supplied, use internal shiny server file names (0,1,2,..)
		if (is.null(filenames))
		{
		  filenames = readFiles
		} else {
		  # it files names are supplied (actual file names that were loaded)
		  # take file names that were actually read and create an index
		  # shiny server saves files with 0,1,2,.. names
		  ind = as.numeric(names(repertoire$data))+1
		  filenames = sapply(unlist(filenames[ind]),file_path_sans_ext)
		}
		# assign file names as names to objects
#		browser()
		# if not all loaded files
		names(mergedData) = names(ntData) = filenames
		return(list(mergedData = mergedData,ntData = ntData))
}

#' cTabPR
#' function to make the count matrix for 1 clone, from merged data
#' @param clone a clone to get counts for
#' @param countData a list of per clone counts for all samples from readMergeSave
#' @param correct a parameter to add to all counts to avoid 0s
#' @return a matrix with counts for a clone across all samples
#'
#'
#### requires 1 input:
###### 1) countData=merged data from readMergeSave
#### correct is a parameter to add to all counts to avoid 0s
cTabPR=function(clone,countData,correct=.5){
    # replace missing values with 0
     minna=function(x){  ### function to deal with missing values
        if(is.na(x)) x=0
        return(x)}
    # get counts for a clone across all samples
     cts=sapply(countData,function(x) return(x[clone]))
     # replace missing values with 0
     cts=sapply(cts,minna)
     # get total read count for each sample minus the count for the clone
     sms=sapply(countData,sum)-cts
     ans=cbind(cts,sms)+correct
     return(ans)
}

#### function to make the data frame for regression
#### requires 3 inputs:
###### 1) cTab = counts object created by cTabPR,
###### a matrix with two columns: counts of a clone and sums minus counts of the clone in all samples
###### 2) peps = vector of peptides corresponding to columns in merged data
###### 3) control=name of control peptide
cDfPR=function(cTab,peps,control="NPA"){
  # if(!(control%in%peps)){ break("control peptide not found")}else{
  datPR=data.frame(as.vector(cTab), ## convert matrix to vector: counts of a clone (c+) come first, then sums minus counts of the clone (c-)
                                    ## then sums minus counts of the clone (c-)
                   factor(rep(c("c+","c-"),rep(nrow(cTab),2)),levels=c("c-","c+")),
                   factor(rep(peps,2),
                          levels=c(control,sort(setdiff(peps,control)))))
  colnames(datPR)=c("cts","clone","pep")
  return(datPR)
}

#### function to convert regression output to OR scale and format
#### requires 1 input:
###### 1) dt = row from mhcMod$coef matrix
######

prepStatsPR=function(dt){  ## dt=row from coefficent matrix

    names(dt)=c("Est","SE","Z","P")
    OR=round(exp(dt["Est"]),3)
    LCB=round(exp(dt["Est"]-1.96*dt["SE"]),3)
    UCB=round(exp(dt["Est"]+1.96*dt["SE"]),3)
    pval=dt["P"]
    ans=c(OR,LCB,UCB,pval)
    names(ans)=c("OR","LCB","UCB","pval")
    return(ans)
}

#### function to organize counts and regression results into 1 data.frame for export
#### messy, surely can be improved
#### requires 2 inputs:
###### 1) cts = counts object created by cTabPR
###### 2) cfs = matrix extracted from model results as
###### coefs[interact,] where
##### coefs=summary(mhcMod)$coef; interact=grep(":",rownames(coefs),fixed=T)
##### probably better to either operate straight on mhcMod or add a function like cTabPR to do so
dfResultPR=function(cts,cfs){
    #rownames(cfs)=gsub(":","|",rownames(cfs),fixed=T)
    ctsM=rbind(c("","counts","sums","",""),cbind(rownames(cts),cts,matrix(rep("",nrow(cts)*2),ncol=2)))
    cfsM=rbind(c("",colnames(cfs)),cbind(rownames(cfs),cfs))
    data=data.frame(rbind(ctsM,matrix(rep("",10),nrow=2),cfsM))
    return(data)
}


#####################################
#' fitModel
#' Perform analysis for one clone
#' @description  The function fits negative binomial regression model for a clone
#' @export
#' @param clone a clone to get counts for
#' @param countData a list of per clone counts for all samples from readMergeSave
#' @param peptides a vector of peptides corresponding to columns in merged data
#' @param control name of control/reference condition
#' @param c.corr a parameter to add to all counts to avoid 0s
#' @return odds ratios and p-values of regression for comparisons to the control condition
#' and comparison of the best to the second best condition
#'


fitModel = function(clone,countData,peptides,control,
                    c.corr=1){

  #### peptides is a vector, equal in length to merged data indicating
  ####  which peptide is represented in each rep.
  #### could extend with a matrix of covariates, rows = length merged data
  require(lme4)
  require(contrast)
  require(multcomp)
  require(MASS)
  require(dplyr)

  #  print(clone)
  ### make the count matrix for a clone across all samples
  ctsPR=cTabPR(clone,countData,correct=c.corr)
  ### make the regression data
  datPR=cDfPR(ctsPR,peps=peptides, control=control)

  ### perform the regression
  #mhcMod <- glm(cts ~ clone*pep, data = datPR, family = poisson(link = "log"))

  # run negative binomial regression
  # if the model fails, return NULL
  tryCatch({
    mhcMod <- glm.nb(cts ~ clone*pep, data = datPR)
  }, error = function(e) {
    print(clone)
    print(e)
    return(NULL)
  })
  # if the model doesn't fail, extract coefficients
  #browser()
  if(!exists("mhcMod")){return(NULL)}

  # if the model fits fine, extract coefficients
  coefs=summary(mhcMod)$coef
  # find interactions to include in the output
  interact=grep(":",rownames(coefs),fixed=T)
  # subset coefficients for interactions that represents results of comparison to control
  coefs_int = coefs[interact,]
  # add condition from which coefficients were extracted
  # and convert coefficients to OR
  condition = gsub("clonec+:pep","",rownames(coefs_int), fixed = T)
  OR = setNames(round(exp(coefs_int[,"Estimate"]),3),
                paste0("OR: ", condition, "_vs_",control))
  pval = setNames(coefs_int[,"Pr(>|z|)"],
                paste0("pval: ", condition, "_vs_",control))

#browser()
  # compare with the second best clone and output OR and p-value
  ### order conditions and find the best and second best interaction coefficients
  pepOrd=coefs_int[order(coefs_int[, "z value"],decreasing=T),]
  # names of the best and the second best conditions
  best=rownames(pepOrd)[1]
  scnd=rownames(pepOrd)[2]

  ### make a contrast matrix “cMat” with columns matching model coefficients,
  # and 1 row
  ### initialize with 0 for all irrelevant conditions,
  # and use the value 1 for the most expanded condition,
  # and -1 for the second most expanded condition
  cMat <- matrix(rep(0,length(mhcMod$coef)), 1) ## initialize contrast matrix
  # set colnames as names of coefficients
  colnames(cMat) = rownames(coefs)
  ### specify contrast top peptide to second
  cMat[1,c(best,scnd)]=c(1,-1)

  ### fit the contrast using glht function, requires the original model,
  # “mhcMod” and the new contrast matrix
  # get an estimate
  pepComp<- glht(mhcMod, linfct=cMat)  ### fit

  ### extract statistics from the contrast result,
  # in this instance just p-value,
  # but OR between those conditions
  bestCond = gsub("clonec+:pep","",best, fixed = T)
  scndCond = gsub("clonec+:pep","",scnd, fixed = T)
  # convert an estimate to OR
  OR_scnd = setNames(round(exp(summary(pepComp)$test$coeff),3),
                "OR: best_vs_second")
  # p-value
  pval_scnd = setNames(summary(pepComp)$test$pvalues,
                  "pval: best_vs_second")
  # names of the best and the second best conditions
  best_vs_second = setNames(paste0(bestCond,"_vs_",scndCond), "second_comparison")

  # combine output and return it
  res = c(clone = clone,OR,pval, OR_scnd, pval_scnd, best_vs_second)

  return(res)

}


#===========
# wrapper for running the full analysis using negative binomial
# from reading files to output all results
#################
#' @export
#' @title runExperiment
#' @description Reads in files with TCR repertoires from
#' a FEST experiment with replicate samples per condition (stimulating peptide).
#' It fits negative binomial model to find expanded clones
#' comparing to a reference samples.
#' It also compares top conditions to find unique expansions.
#' The results are return and saved in an Excel file.
#' @param files a list of filenames with full paths
#' @param peptides a vector of peptides corresponding to columns in merged data
#' @param ctThresh minimal number of reads required to consider a clone
#' @param control name of control/reference condition
#' @param ORthr threshold for OR to consider a clone expanded
#' @param FDRthr threshold for FDR to consider a clone expanded
#' @param excludeCond a vector of conditions to exclude from the analysis
#' @param xrCond a vector of cross-reactive conditions
#' @param percentThr a threshold for percentage of reads in a sample to consider a clone expanded
#' @param outputFile name of the output file
#' @param saveToFile logical, if TRUE save results to a file
#' @param permute logical, if TRUE permute sample labels to run a permutation test
#' @return a list of all expanded clones, uniquely expanded clones
#' and parameters of the run
#'
#### returns a list of all expanded clones, uniquely expanded clones
#### and parameters of the run

#### requires the following inputs:
#### main arguments are filenames (full paths in current implementation)
#### and a vector of peptide ids for each file.
#### additional arguments include a minimal number of reads required
#### to consider a clone
#### the control peptide ID,
#### a vector of conditions to exclude from the analysis
#### a threshold for OR and FDR to consider a clone expanded
#### a vector of cross-reactive conditions
#### added option for permuting labels to answer to Kellie's request. Need to remove for the app and package

#### 2025-01-29 added comparison to the second best to find unique expansions

runExperiment=function(files, peptides, ctThresh=50, control,
                       ORthr=1, FDRthr = 0.05, excludeCond = NA,
                       xrCond = NA, percentThr = 0,
                       outputFile = "output.xlsx",
                       saveToFile = T, permute = FALSE){


  #### start algorithm read data
  mergeDat=readMergeSave(files, filenames = NULL)$mergedData

  # permute sample labels in mergeDat to run permutation test
  if (permute){
    set.seed(123456)
    # create sampling
    s = sample(1:length(mergeDat),size = length(mergeDat), replace = F)
    # update sample names
    names(mergeDat) = names(mergeDat)[s]
    # update peptides
    peptides = peptides[s]
  }

  # get clones to test
  goodClones = getClonesToTest(mergeDat, ctThresh = ctThresh)
  print(c("good clones #",length(goodClones)))

  # run the analysis for selected clones
  fitResults = fitModelSet(goodClones, mergeDat, peptides,
                       excludeCond = excludeCond,
                       control=control,c.corr=1)
  rownames(fitResults) = fitResults$clone

#browser()
  # get all expanded clones relative to the control
  # find colunms with control
  contCol = grep(control,colnames(fitResults), value = T)
  res_exp = getExpanded(fitResults[,c("clone",contCol)], mergeDat,
                        ORthr = ORthr, FDRthr = FDRthr)

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
  res_exp05 = res_exp[apply(res_exp[,percCol],1,max) >percentThr,]

   #======
  # find uniquely expanded clones by checking the second best clone
  # save the second best comparison results
  screen_scndBest = fitResults[rownames(res_exp05), setdiff(colnames(fitResults),contCol)]
  # check for uniqueness. it should be expanded and significant in comparison to the second best as well
  unique_exp = (screen_scndBest[,grep("OR", colnames(screen_scndBest))]>1 &
    screen_scndBest[,grep("FDR", colnames(screen_scndBest))]<FDRthr)
  # merge the results
  # add columns that indicates how many and what comparisons were significant
  # get T/F matrix for uniqueness comparisons
  # results for the second best comparison
  # the rest of info
  sigComCol = c("clone","n_significant_comparisons","significant_comparisons")
  screen_scndBest = cbind(res_exp05[,sigComCol],
                          unique_exp,
                          res_exp05[,setdiff(colnames(res_exp),sigComCol)],
                          screen_scndBest)


  # get clones expanded in one condiiton comparing to the control
   res_uniq = res_exp05 %>% filter(n_significant_comparisons == 1)
 #=====
  # find cross-reactive clones
  if(!is.na(xrCond))
  {
    res_xr = getXR(res_exp05, peptides, xrCond = xrCond)
    res_uniq = rbind(res_uniq,res_xr)
  }
  #=====
  # get parameters of the run
  # total reads for each sample
  # get total read count for each sample
  totalReadCountPerSample = sapply(mergeDat, sum)
  totalReads = data.frame(parameter = paste(names(totalReadCountPerSample), "total_reads", sep = "_"),
                          value = totalReadCountPerSample)
  # the parameters
  params = data.frame(parameter = c("reference","read_threshold",
                                    "FDR_threshold",
                                    "percent_threshold",
                                    "exclude_conditions",
                                    "cross_reactive_conditions"),
                      value = c(control,ctThresh,FDRthr,percentThr,
                                paste(excludeCond,collapse = ","),
                                paste(xrCond,collapse = ",")))

  output = list(ref_only = res_exp,
                uniquely_exp = res_uniq,
                second_best = screen_scndBest,
                params = rbind(params,totalReads))
  # save into Excel file
  if(saveToFile)
  {
    saveResults(output, outputFile = outputFile)
  } else return(output)
}

#' @export
#' @title getAbundances
#' @description returns the read count for clones of interest in all samples
#' @param clones a vector of clones of interest
#' @param countData a list of counts for all samples
#' @return a matrix of clone abundances across samples in `countData`
# input: a list of merged data, a vector of clones of interest
getAbundances = function(clones,countData)
{
  # create output matrices
  output_counts = matrix(0,nrow = length(clones), ncol = length(countData))
  rownames(output_counts) = clones
  colnames(output_counts) = names(countData)

  # get the read count for clones in each sample
  for (i in names(countData))
  {
    rows = intersect(names(countData[[i]]),rownames(output_counts))
    output_counts[rows,i] = countData[[i]][rows]
  }
  # update colnames to add "abundance"
  colnames(output_counts) = paste(names(countData),'abundance', sep = '_')
  return(output_counts)
}

#' @title getClonesToTest
#' @description function that returns the read counts for clones of interest in all samples
#' @param countDat a list of merged data
#' @param ctThresh a minimal number of reads required to consider a clone
#' @return a vector of clones of interest
# function that returns the read counts for clones of interest in all samples
# input: a list of merged data, a vector of clones of interest
# all clones from all samples
getClonesToTest = function(countDat, ctThresh = 50)
{
  # all clones
  clones=unlist(sapply(countDat,function(x)  return(names(x))))
  # the corresponding counts
  cts=unlist(sapply(countDat,function(x)  return(x)))
  #
  sumCt=tapply(cts,clones,sum)

  # clones that have more than ctThresh reads to run the analysis
  goodClones=names(sumCt)[which(sumCt>ctThresh)]

  return(goodClones)
}

#' fit model for a set of clones
#' @param clones a list of clones to fit the model
#' @param countData a list of counts for all samples
#' @param peptides a vector of peptides corresponding to columns in merged data
#' @param excludeCond a vector of conditions to exclude from the analysis
#' @return a matrix with all ORs, p-values and FDRs
#'
# return a matrix with all ORs, p-values and FDRs
# input: a list of clones, merged data, peptides,

fitModelSet = function(clones, countData, peptides, excludeCond = NA,...)
{
    # run model for "good" clones
  if (length(excludeCond)>0)
  {
    # get indexes of conditions to include in the analysis
    incl = which(!(peptides %in% excludeCond))
    cat("Conditions to include:", peptides[incl],"\n")
    fitResults=lapply(clones,fitModel,countData=countData[incl],
                  peptides=peptides[incl],...)
  }else{
    fitResults=lapply(clones,fitModel,countData=countData,
                  peptides=peptides,...)
  }

 #  browser()

  # remove enties with NULL
  fitResults = fitResults[!sapply(fitResults, is.null)]
  # convert list to a data.frame
  fitResults = as.data.frame(bind_rows(fitResults))

  # add FDR adjustment
  # find columns with p-values
  pvalCol = grep("pval",colnames(fitResults), value = T)
  fdrs = c()
  for(i in pvalCol)
  {
    fdrs = cbind(fdrs,p.adjust(as.numeric(fitResults[,i]), method = "BH"))
  }
  # add colnames for fdrs and add to the fitResults matrix
  colnames(fdrs) = paste0("FDR:",gsub("pval:","",pvalCol))
  fitResults = cbind(fitResults, fdrs)


  return(fitResults)
}

#' function that returns the expanded clones
#' @param fitResults a data frame with ORs, p-values and FDRs
#' @param countData a list of counts for all samples
#' @param ORthr threshold for OR to consider a clone expanded
#' @param FDRthr threshold for FDR to consider a clone expanded
#' @return a data frame with expanded clones
#'
# input: a data frame with ORs, p-values and FDRs
# output: a data frame with expanded clones
getExpanded = function(fitResults, countData, ORthr = 1, FDRthr = 0.05)
{
  # find all significantly expanded clones
  # find comparison names
  comp = grep("OR",colnames(fitResults), value = T)
  comp = gsub("OR: ","",comp)

  expandedClones = c()
  for (i in comp){
    #print(i)
    expandedClones = union(expandedClones,
                           rownames(fitResults)[which(as.numeric(fitResults[,paste0("OR: ",i)]) >= ORthr & # expanded
                                                    as.numeric(fitResults[,paste0("FDR: ",i)]) < FDRthr)])# significant
  }
  # get the results for expanded clones only
  res_exp = fitResults[expandedClones,]
  # add columns that indicates how many and what comparisons were significant
  # get T/F matrix for significant comparisons
  sig = (res_exp[,paste0("FDR: ",comp)] < FDRthr &
           res_exp[,paste0("OR: ",comp)] > ORthr)
 # browser()
  # list significant comparisons
  sigComp = apply(sig, 1, function(x){
    # select significant comparisons and get first condition before "vs"
    s = sapply(strsplit(comp[x], split = "_vs_"), getElement, 1)
    # list significant comparisons using comma
    paste(s, collapse = ",")
    })
  # add the number and the list to the results
  res_exp = cbind(clone = res_exp[,"clone"],
                  n_significant_comparisons = rowSums(sig, na.rm = T),
                  significant_comparisons = sigComp,
                  res_exp[,setdiff(colnames(res_exp),c("clone"))])
   # add abundance and percentage of the top clones in each condition
  # get abundance for the top clones
  abundance = getAbundances(rownames(res_exp), countData)
  # get total read count for each sample
  totalReadCountPerSample = sapply(countData, sum)
  # calculate the percentage of each clone in each sample
  percentage = round(sweep(abundance, 2, totalReadCountPerSample, "/")*100,3)
  colnames(percentage) = paste(names(countData),'percent', sep = '_')

  res_exp = cbind(res_exp,
                  abundance[rownames(res_exp),],
                  percentage[rownames(res_exp),])

  return(res_exp)

}

# function that returns the cross-reactive clones
# input: a data frame with expanded clones
# and a vector of cross-reactive conditions
# output: a data frame with cross-reactive clones

getXR = function(res, conditions, xrCond = xrCond)
{
  # a vector of conditions that shouldn't be cross-reactive
  # find all and take difference
  allCond = res %>% filter(n_significant_comparisons == 1) %>%
    dplyr::select(significant_comparisons) %>% unlist() %>% unique()

  # conditions that shouldn't be cross-reactive
  excludeCond = setdiff(conditions, xrCond)

  # find cross-reactive clones
  res_xr = res %>% filter(n_significant_comparisons > 1)
  # find clone that are cross-reactive fo xrCond
  inclXR = res_xr %>% filter(grepl(paste(xrCond,collapse = "|"),significant_comparisons))

  # exclude conditions that are not in xrCond
  excl = inclXR %>% filter(!grepl(paste(excludeCond,collapse = "|"),significant_comparisons))
  return(excl)
}

# function that saves the results to an excel file
# input: a list of data frames with results

saveResults = function(results, outputFile = "output.xlsx")
{
  library(openxlsx)
  # create a workbook
  wb = createWorkbook()
  # add sheets
  addWorksheet(wb, "expanded")
  addWorksheet(wb, "uniquely_expanded")
  addWorksheet(wb, "second_best_result")
  addWorksheet(wb, "parameters")
  # add data to sheets
  writeData(wb, "expanded", results$ref_only)
  writeData(wb, "uniquely_expanded", results$uniquely_exp)
  writeData(wb, "second_best_result", results$second_best)
  writeData(wb, "parameters", results$params)
  # save the workbook
  saveWorkbook(wb, outputFile, overwrite = TRUE)
}

#' function that extracts condition and replicate information
#' from the file names. The condition and replicate should
#' be separated by "_" and be the last two elements.
#' @param filenames a vector of file names
#' @return a data frame with condition and replicate information

splitFileName = function(filenames)
{
  # split the file names by "_"
  splitNames = strsplit(filenames, "_")
  # get the last two elements
  condRep = sapply(splitNames, function(x) x[(length(x)-2):length(x)])
  # create a data frame with condition and replicate information
  condRep = data.frame(sample = condRep[1,],
                       condition = condRep[2,],
                       replicate = condRep[3,])
  return(condRep)
}
