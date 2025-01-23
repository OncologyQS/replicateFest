# Functions for analysis of Manafest data in shiny app
# using input with replicates
# based on Leslie's development for the Cervical SPORE
#
geti = function(x,i){return(x[i])}


# reads files, removes non-productive sequencies, extracts counts,
# creates all necessary objects, and saves them
# input: a path to folder with input files
# output:
# mergedData is AA level counts (merged by AA sequences)
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


#### function to make the count matrix for 1 clone, from merged data
#### requires 1 input:
###### 1) mergeData=merged data from readMergeSave
#### correct is a parameter to add to all counts to avoid 0s
cTabPR=function(clone,mergeData,correct=.5){
    # replace missing values with 0
     minna=function(x){  ### function to deal with missing values
        if(is.na(x)) x=0
        return(x)}
    # get counts for a clone across all samples
     cts=sapply(mergeData,function(x) return(x[clone]))
     # replace missing values with 0
     cts=sapply(cts,minna)
     # get total read count for each sample minus the count for the clone
     sms=sapply(mergeData,sum)-cts
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
#### master function to perform analysis for one clone.
#### adapted from Leslie's function
#### returns results of regression for the best and the second best clone

#### evaluates model for one clone at a time, use sapply for a set of clones
#### requires the following inputs
###### mergedData=merged data from readMergeSave, a list of clone counts per condition
######.   called mergeData elsewhere, change?
###### peptides=vector of peptides corresponding to columns in merged data
###### control=name of control/reference condition
#### c.corr is a parameter to add to all counts to avoid 0s

fitModel = function(clone,mergedData,peptides,control="NPA",
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
  ctsPR=cTabPR(clone,mergedData,correct=c.corr)
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
  # subset coefficients for interactions
  coefs = coefs[interact,]
  # add condition from which coefficients were extracted
  # and convert coefficients to OR
  condition = gsub("clonec+:pep","",rownames(coefs), fixed = T)
  OR = setNames(round(exp(coefs[,"Estimate"]),3),
                paste0("OR: ", condition, "_vs_",control))
  pval = setNames(coefs[,"Pr(>|z|)"],
                paste0("pval: ", condition, "_vs_",control))

  res = c(clone = clone,OR,pval)

  return(res)

}


#===========
# wrapper for running the full analysis from reading files to output all results
#################
#### this function runs all files in an experiment
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

#### TODO 2025-01-16 add comparison to the second best to find unique clones

runExperiment=function(files, peptides, ctThresh=50,cont="NP",
                       ORthr=1, FDRthr = 0.05, excludeCond = NA,
                       xrCond = NA, percentThr = 0, outputFile = "output.xlsx",
                       saveToFile = T, permute = FALSE){


  #### start algorithm read data
  mergeDat=readMergeSave(files, filenames = NULL)$mergedDat

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
  screen = fitModelSet(goodClones, mergeDat, peptides,
                       excludeCond = excludeCond,
                       control=cont,c.corr=1)
  rownames(screen) = screen$clone
  # get all expanded clones
  res_exp = getExpanded(screen, mergeDat,
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
  res_uniq = res_exp05 %>% filter(n_significant_comparisons == 1)

  #=====
  # find cross-reactive clones
  if(length(xrCond)>0)
  {
    res_xr = getXR(res_exp05, xrCond = xrCond)
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
                      value = c(cont,ctThresh,FDRthr,percentThr,
                                paste(excludeCond,collapse = ","),
                                paste(xrCond,collapse = ",")))

  output = list(ref_only = res_exp,
                uniquely_exp = res_uniq,
                params = rbind(params,totalReads))
  # save into Excel file
  if(saveToFile)
  {
    saveResults(output, outputFile = outputFile)
  } else return(output)
}


# function that returns the read count for clones of interest in all samples
# input: a list of merged data, a vector of clones of interest
getAbundances = function(clones,mergedData)
{
  # create output matrices
  output_counts = matrix(0,nrow = length(clones), ncol = length(mergedData))
  rownames(output_counts) = clones
  colnames(output_counts) = names(mergedData)

  # get the read count for clones in each sample
  for (i in names(mergedData))
  {
    rows = intersect(names(mergedData[[i]]),rownames(output_counts))
    output_counts[rows,i] = mergedData[[i]][rows]
  }
  # update colnames to add "abundance"
  colnames(output_counts) = paste(names(mergedData),'abundance', sep = '_')
  return(output_counts)
}

# function that returns the read count for clones of interest in all samples
# input: a list of merged data, a vector of clones of interest
# all clones from all samples
getClonesToTest = function(mergeDat, ctThresh = 50)
{
  # all clones
  clones=unlist(sapply(mergeDat,function(x)  return(names(x))))
  # the corresponding counts
  cts=unlist(sapply(mergeDat,function(x)  return(x)))
  #
  sumCt=tapply(cts,clones,sum)

  # clones that have more than ctThresh reads to run the analysis
  goodClones=names(sumCt)[which(sumCt>ctThresh)]

  return(goodClones)
}

# fit model for a set of clones and
# return a matrix with all ORs, p-values and FDRs
# input: a list of clones, merged data, peptides,

fitModelSet = function(clones, mergedData, peptides, excludeCond = NA,...)
{
    # run model for "good" clones
  if (length(excludeCond)>0)
  {
    # get indexes of conditions to include in the analysis
    incl = which(!(peptides %in% excludeCond))
    cat("Conditions to include:", peptides[incl],"\n")
    screen=lapply(clones,fitModel,mergedData=mergedData[incl],
                  peptides=peptides[incl],...)
  }else{
    screen=lapply(clones,fitModel,mergedData=mergedData,
                  peptides=peptides,...)
  }

  #  browser()

  # transpose to have clones in rows
  #  screen = data.frame(t(screen), check.names = F)
  # remove enties with NULL
  screen = screen[!sapply(screen, is.null)]
  # convert list to a data.frame
  screen = as.data.frame(bind_rows(screen))

  # add FDR adjustment
  # find columns with p-values
  pvalCol = grep("pval",colnames(screen), value = T)
  fdrs = c()
  for(i in pvalCol)
  {
    fdrs = cbind(fdrs,p.adjust(as.numeric(screen[,i]), method = "BH"))
  }
  # add colnames for fdrs and add to the screen matrix
  colnames(fdrs) = paste0("FDR:",gsub("pval:","",pvalCol))
  screen = cbind(screen, fdrs)


  return(screen)
}

# function that returns the expanded clones
# input: a data frame with ORs, p-values and FDRs
# output: a data frame with expanded clones
getExpanded = function(screen, mergedData, ORthr = 1, FDRthr = 0.05)
{
  # find all significantly expanded clones
  # find comparison names
  comp = grep("OR",colnames(screen), value = T)
  comp = gsub("OR: ","",comp)

  expandedClones = c()
  for (i in comp){
    #print(i)
    expandedClones = union(expandedClones,
                           rownames(screen)[which(as.numeric(screen[,paste0("OR: ",i)]) >= ORthr & # expanded
                                                    as.numeric(screen[,paste0("FDR: ",i)]) < FDRthr)])# significant
  }
  # get the results for expanded clones only
  res_exp = screen[expandedClones,]
  # add columns that indicates how many and what comparisons were significant
  # get T/F matrix for significant comparisons
  sig = (res_exp[,paste0("FDR: ",comp)] < FDRthr &
           res_exp[,paste0("OR: ",comp)] > ORthr)
 # browser()
  # list significant comparisons
  sigComp = apply(sig, 1, function(x){
    # select significant comparisons and get first condition before "vs"
    s = sapply(strsplit(comp[x], split = "_vs_"), geti, 1)
    # list significant comparisons using comma
    paste(s, collapse = ",")
    })
  # add the number and the list to the results
  res_exp = cbind(res_exp,
                  n_significant_comparisons = rowSums(sig, na.rm = T),
                  significant_comparisons = sigComp)
   # add abundance and percentage of the top clones in each condition
  # get abundance for the top clones
  abundance = getAbundances(rownames(res_exp), mergedData)
  # get total read count for each sample
  totalReadCountPerSample = sapply(mergedData, sum)
  # calculate the percentage of each clone in each sample
  percentage = round(sweep(abundance, 2, totalReadCountPerSample, "/")*100,3)
  colnames(percentage) = paste(names(mergedData),'percent', sep = '_')

  res_exp = cbind(res_exp,
                  abundance[rownames(res_exp),],
                  percentage[rownames(res_exp),])

  return(res_exp)

}

# function that returns the cross-reactive clones
# input: a data frame with expanded clones
# and a vector of cross-reactive conditions
# output: a data frame with cross-reactive clones

getXR = function(res, xrCond = xrCond)
{
  # a vector of conditions that shouldn't be cross-reactive
  # find all and take difference
  allCond = res %>% filter(n_significant_comparisons == 1) %>%
    dplyr::select(significant_comparisons) %>% unlist() %>% unique()

  excludeCond = setdiff(unique(timeSamples$Condition), xrCond)

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
  addWorksheet(wb, "parameters")
  # add data to sheets
  writeData(wb, "expanded", results$ref_only)
  writeData(wb, "uniquely_expanded", results$uniquely_exp)
  writeData(wb, "parameters", results$params)
  # save the workbook
  saveWorkbook(wb, outputFile, overwrite = TRUE)
}
