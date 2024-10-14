# Functions for analysis of Manafest data in shiny app
#
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
			print(i)
		  dat = repertoire$data[[i]]
				#count reads of productive sequences only
				#mergedData[[i]] = tapply(dat[[,'Clones']], dat[[,'CDR3.aa']], sum, na.rm = T)
				mergedData[[i]] = tapply(dat$Clones, dat$CDR3.aa, sum, na.rm = T)
				# nucleotide level data
				#ntData[[i]] = tapply(dat[,'Clones'], dat[,'CDR3.nt'], sum, na.rm = T)
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
#### pois in name here and PR elsewhere is legacy of
#### earlier implementation using Poisson regression.
####  should we change to NBR or the like for negative binomial?

### has two modes, a screening mode which produces
#### a simplified result for filtering promising probes.

### a final mode
##### a little bit of a Frankenfunction in the end
##### should probably be broken down further
#### evaluates model for one clone at a time, use sapply for a set of clones
#### requires 3 inputs
###### 1) mergedData=merged data from readMergeSave, a list of clone counts per condition
######.   called mergeData elsewhere, change?
###### 2) peptides=vector of peptides corresponding to columns in merged data
###### 3) control=name of control peptide
#### c.corr?

poisReg=function(clone,mergedData,peptides,control="NPA",
                 c.corr=1,screen=F,printDetail=F){

#### peptides is a vector, equal in length to merged data indicating
####  which peptide is represented in each rep.
#### could extend with a matrix of covariates, rows = length merged data
    require(lme4)
    require(contrast)
    require(multcomp)
    require(MASS)
  require(dplyr)
  
#  print(clone)
### make the count matrix
    ctsPR=cTabPR(clone,mergedData,correct=c.corr)
    ctsPR0=ctsPR-c.corr
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
    coefs=summary(mhcMod)$coef
    # find interactions to include in the output
    interact=grep(":",rownames(coefs),fixed=T)

    # order coefficients by z value to find the best peptide
    pepOrd=interact[order(coefs[interact, "z value"],decreasing=T)]
    # get best clone
    best=pepOrd[1]
    # in what condition is the best peptide?
    bestPep=sub("clonec+:pep","",rownames(coefs)[best],fixed=T)
    # find the second best peptide
    scnd=pepOrd[2]
    
    # calculate the difference in counts between the two best peptides
    bCts=ctsPR[which(peptides==bestPep),,drop=F]
    bRts=sort(bCts[,1]/apply(bCts,1,sum),decreasing=T)
    ctDisc=bRts[2]/bRts[1]
    names(ctDisc)=""
    
#browser()    
    # if running in screening mode, return a simplified result
 if(screen==T){
   # summary of model fitting results
   s = summary(mhcMod)
    ans=c("pep"=rownames(coefs)[best],
          prepStatsPR(s$coef[best,]),
          "coef2"=s$coef[scnd,1],
          "p2"=s$coef[scnd,4],
          "ctDiff"=ctDisc)
    #names(ansP)=c("pep","OR","LCB","UCB","pval","coef2","p2","ctDiff")

 }else{ ### if running in final mode, return a detailed result

     detail=dfResultPR(ctsPR0,coefs[interact,])

     if(printDetail){
     #write.xlsx(detail,file=outFile,sheetName=clone,col.names=F,row.names=F,overwrite=F)

     }

   res=c(rownames(coefs)[best],prepStatsPR(summary(mhcMod)$coef[best,]),
         ctsPR0[,1])
     res[1]=sub("clonec+:pep","",res[1],fixed=T)
     names(res)[1]="peptide"
    ans=list(res,detail)}
    
    return(ans)

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
  ctsPR0=ctsPR-c.corr
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
  coefs=summary(mhcMod)$coef
  # find interactions to include in the output
  interact=grep(":",rownames(coefs),fixed=T)
  
  # order coefficients by z value to find the best peptide
  pepOrd=interact[order(coefs[interact, "z value"],decreasing=T)]
  # get best clone
  best=pepOrd[1]
  # in what condition is the best peptide?
  bestPep=sub("clonec+:pep","",rownames(coefs)[best],fixed=T)
  # find the second best peptide
  scnd=pepOrd[2]
  
  # calculate the difference in counts between the two best peptides
  bCts=ctsPR[which(peptides==bestPep),,drop=F]
  bRts=sort(bCts[,1]/apply(bCts,1,sum),decreasing=T)
  ctDisc=bRts[2]/bRts[1]
  names(ctDisc)=""
  
  #browser()    
    # summary of model fitting results
    s = summary(mhcMod)
    ans=c("pep"=rownames(coefs)[best], # the name of the best peptide
          prepStatsPR(s$coef[best,]), # OR, CI, p-values for the best
          scnd = prepStatsPR(s$coef[scnd,])[c("OR","pval")], # OR and p-value for the second best
          "coef2"=s$coef[scnd,1],
          "p2"=s$coef[scnd,4],
          "ctDiff"=ctDisc) # the difference in counts between the two best peptides

  return(ans)
  
}

#################
#### this function does screening and final reporting using criteria
#### developed for Cervical SPORE project.
#### main arguments are filenames (full paths in current implementation)
#### and a vector of peptide ids for each file.
#### additional arguments include a minimal number of reads required
#### to consider a clone
#### the control peptide ID,
#### and an identifying character string for the output file name.


runRM=function(files, gps,ctThresh=50,cont="NP",outF="manafest"){
 require(openxlsx)
#### initiate output file
    of=paste(outF,"Cands.xlsx",sep="_")
    if (file.exists(of))  file.remove(of)


#### initiate xlsx workbook

    wBook=createWorkbook()

### and main sheet for the workbook
    addWorksheet(wBook,"summary")


#### save wBook for adding detail sheets and acutal summary later



#### start algorthm read data
    mergeDat=readMergeSave(files, filenames = NULL)$mergedDat


    clones=unlist(sapply(mergeDat,function(x)  return(names(x))))
 cts=unlist(sapply(mergeDat,function(x)  return(x)))
 readSums=sapply(mergeDat,sum)
    maxCt=tapply(cts,clones,max)
    sumCt=tapply(cts,clones,sum)

    #cloneCts=table(clones)
    goodClones=names(sumCt)[which(sumCt>ctThresh)]
# run model for "good" clones
    screen=sapply(goodClones,poisReg,mergedData=mergeDat,
                  peptides=gps,control=cont,c.corr=1,screen=T)

    
    screenM=matrix(as.numeric(t(screen[-1,])),ncol=nrow(screen)-1)
    rownames(screenM)=screen[1,]

    colnames(screenM)=rownames(screen)[-1]
    # find results to include in the output
topPbs=which(!rownames(screenM)%in%paste0("clonec+:",cont)&
               as.numeric(screenM[,"OR"])>1 &  # expanded
               p.adjust(as.numeric(screenM[,"pval"]))<0.001 & # significant
               (as.numeric(screenM[,"coef2"])<0 | 
                  as.numeric(screenM[,"p2"])>0.01))# and negative or not significant
               

topPbsF=goodClones[topPbs[order(as.numeric(screenM[topPbs,"ctDiff"])<.1,as.numeric(screenM[topPbs,"pval"]),decreasing=F)]]

 # fullRes=data.frame(t(sapply(topPbsF,poisReg,mergedData=mergeDat,peptides=gps,control="NP",c.corr=1,screen=F,printDetail=T,outFile=of)))

 fullRes=matrix(rep(NA,(length(topPbsF)+1)*(2*length(mergeDat)+5)),
                nrow=length(topPbsF)+1)
 rownames(fullRes)=c("denominators",topPbsF)
 colnames(fullRes)=c("pep","OR","LCB","UCB","pval",
                     paste(names(mergeDat),"abundance", sep="_"),
                     paste(names(mergeDat),"percent", sep="_"))
 fullRes[1,-(1:5)]=as.character(readSums)
#browser() 
 for(cln in topPbsF){
     ans=poisReg(cln,mergedData=mergeDat,peptides=gps,control=cont,c.corr=1,screen=F,printDetail=T)
    # results from the regression and abundances in all conditions
     fullRes[cln,1:length(ans[[1]])]=ans[[1]]
     # add percentage of that clone in each condition
     fullRes[cln,(length(ans[[1]])+1):ncol(fullRes)]=
       round((as.numeric(ans[[1]][6:(length(ans[[1]]))])/readSums)*100,3)
     
#     addWorksheet(wBook,cln)
# writeData(wBook,sheet=cln,x=ans[[2]],colNames=F,rowNames=F)
     }

    writeData(wBook,sheet="summary",x=fullRes,colNames=T,rowNames=T)
    saveWorkbook(wBook, of,overwrite=T)

   # write.xlsx(fullRes,file=of,col.names=T,row.names=T,sheetName="S",append=T)
  return(fullRes)
}

#===========
# new function for running the analysis
#################
#### this function runs all files in an experiment
#### finds expanded clones and saves them in an excel file
#### then finds unqiuely expanded clones and saves them in an excel file

#### main arguments are filenames (full paths in current implementation)
#### and a vector of peptide ids for each file.
#### additional arguments include a minimal number of reads required
#### to consider a clone
#### the control peptide ID,
#### and an identifying character string for the output file name.


runExperiment=function(files, gps, ctThresh=50,cont="NP",outF="manafest"){
  require(openxlsx)
  #### initiate output file
  of=paste(outF,"Cands.xlsx",sep="_")
  if (file.exists(of))  file.remove(of)
  
  
  #### initiate xlsx workbook
#  wBook=createWorkbook()
  
  ### and main sheet for the workbook
#  addWorksheet(wBook,"ref_comparison_only")
  
 #### save wBook for adding detail sheets later
  
  #### start algorithm read data
  mergeDat=readMergeSave(files, filenames = NULL)$mergedDat
  
  
  clones=unlist(sapply(mergeDat,function(x)  return(names(x))))
  cts=unlist(sapply(mergeDat,function(x)  return(x)))
  readSums=sapply(mergeDat,sum)
  maxCt=tapply(cts,clones,max)
  sumCt=tapply(cts,clones,sum)
  
  # clones that have more than ctThresh reads to run the analysis
  goodClones=names(sumCt)[which(sumCt>ctThresh)]
  # run model for "good" clones
  screen=sapply(goodClones,fitModel,mergedData=mergeDat,
                peptides=gps,control=cont,c.corr=1)
  
  
  screenM=matrix(as.numeric(t(screen[-1,])),ncol=nrow(screen)-1)
  # set clones as rownames
  rownames(screenM) = colnames(screen)
  # set colnames
  colnames(screenM)=rownames(screen)[-1]
  # add condition
  screenM = cbind(condition = gsub("clonec+:pep","",screen["pep",], fixed = T), 
                  screenM)
#browser()
  # add FDR adjustment
  screenM = cbind(screenM, FDR =p.adjust(as.numeric(screenM[,"pval"]), method = "BH"))
  # find results to include in the output
  topPbs=which(as.numeric(screenM[,"OR"])>1 &  # expanded
               as.numeric(screenM[,"FDR"])<0.001 & # significant
               (as.numeric(screenM[,"scnd.OR"])<1 | 
                  as.numeric(screenM[,"scnd.pval"])>0.01))# and negative or not significant
  
  res = screenM[topPbs,]
  # order clones
  res = res[order(as.numeric(res[,"ctDiff"])<.1,
                  as.numeric(res[,"pval"]),decreasing=F),]

  return(res)
  
  fullRes=matrix(rep(NA,(length(topPbsF)+1)*(2*length(mergeDat)+5)),
                 nrow=length(topPbsF)+1)
  rownames(fullRes)=c("denominators",topPbsF)
  colnames(fullRes)=c("pep","OR","LCB","UCB","FDR",
                      paste(names(mergeDat),"abundance", sep="_"),
                      paste(names(mergeDat),"percent", sep="_"))
  fullRes[1,-(1:5)]=as.character(readSums)
  #browser() 
  # for(cln in topPbsF){
  #   ans=poisReg(cln,mergedData=mergeDat,peptides=gps,control=cont,c.corr=1,screen=F,printDetail=T)
  #   # results from the regression and abundances in all conditions
  #   fullRes[cln,1:length(ans[[1]])]=ans[[1]]
  #   # add percentage of that clone in each condition
  #   fullRes[cln,(length(ans[[1]])+1):ncol(fullRes)]=
  #     round((as.numeric(ans[[1]][6:(length(ans[[1]]))])/readSums)*100,3)
  #   
  #   #     addWorksheet(wBook,cln)
  #   # writeData(wBook,sheet=cln,x=ans[[2]],colNames=F,rowNames=F)
  # }
  
#  writeData(wBook,sheet="summary",x=fullRes,colNames=T,rowNames=T)
#  saveWorkbook(wBook, of,overwrite=T)
  
  # write.xlsx(fullRes,file=of,col.names=T,row.names=T,sheetName="S",append=T)
}

