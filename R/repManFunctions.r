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


#############################
#### replicate manafest uses negative binomial regression
#### should one line of code even be a function?
#### requires 1 input:
###### 1) mergeData=merged data from readMergeSave
readCount=function(mergeData){  ### calculate total read count for each sample, from merged data
    reads=sapply(mergeData,function(x) return(sum(x)))
    return(reads)
    }


#### function to make the count matrix for 1 clone, from merged data
#### requires 1 input:
###### 1) mergeData=merged data from readMergeSave
cTabPR=function(clone,mergeData,correct=.5){
     minna=function(x){  ### function to deal with missing values
        if(is.na(x)) x=0
        return(x)}

     cts=sapply(mergeData,function(x) return(x[clone]))
     cts=sapply(cts,minna)
     sms=readCount(mergeData)-cts
     ans=cbind(cts,sms)+correct
     return(ans)
}

#### function to make the data frame for regression
#### need to change this to explicitly make the control sample the intercept
#### rather than implicitly using naming ("000")
#### requires 3 inputs:
###### 1) cTab = counts object created by cTabPR
###### 2) peps = vector of peptides corresponding to columns in merged data
###### 3) control=name of control peptide
cDfPR=function(cTab,peps,control="NPA"){
   # if(!(control%in%peps)){ break("control peptide not found")}else{
datPR=data.frame(as.vector(cTab),rep(c("c+","c-"),rep(nrow(cTab),2)),rep(peps,2))
colnames(datPR)=c("cts","clone","pep")
datPR$pep[which(datPR$pep==control)]=paste0("000",control)
return(datPR)
}

#### new function to make data frame for regression
#### sets factor levels to make the control peptide and clone c- the null values
#### rather than implicitly using naming ("000")
#### not tested e.o.d. Aug 6, 2024
#### requires 3 inputs:
###### 1) cTab = counts object created by cTabPR
###### 2) peps = vector of peptides corresponding to columns in merged data
###### 3) control=name of control peptide
cDfPRv2=function(cTab,peps,control="NPA"){
  # if(!(control%in%peps)){ break("control peptide not found")}else{
  datPR=data.frame(as.vector(cTab),
                   factor(rep(c("c+","c-"),rep(nrow(cTab),2)),
                          levels=c("c-","c+")),
                   factor(rep(peps,2),levels=c(control,sort(setdiff(peps,control)))))
  colnames(datPR)=c("cts","clone","pep")
  return(datPR)
}

#### function to put the peptides into order from most expanded to least.
#### requires 1 input:
###### 1) coefs = model coefficeint matrix (mhcMod$coef)
pepOrdPR=function(coefs){
    intCoefLocs=grep("clonec+:",names(coefs),fixed=T)
    names(intCoefLocs)=sub("clonec+:pep","",names(coefs)[intCoefLocs],fixed=T)
    intCoefLocs=intCoefLocs[order(coefs[intCoefLocs],decreasing=T)]
    return(intCoefLocs)}

#### function to convert regression output to OR scale and format
#### requires 1 input:
###### 1) dt = row from mhcMod$coef matrix
######

prepStatsPR=function(dt){  ## dt=row from coefficent matrix
    ##dt=summary(mhcMod)$coef[pepOrd[pep],]
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

#### function to export counts as a vector. Changed form a few times, this version is pretty pointless
#### should one line of code even be a function?
#### still includes code for reporting frequencies instead of or in addition to counts,
## which should be reinstated as agreed in call with Kellie Smith
#### requires 2 inputs
###### 1) cts = counts object created by cTabPR
###### 2) groups = vector of peptides corresponding to columns in merged data,
######.   called 'peps' or 'peptides' in most functions, maybe change?
vResultPR=function(cts,groups){

    #freqs=cts[,1]/apply(cts,1,sum)
    #freqs=paste(cts[,1],apply(cts,1,sum),sep="/")
   # names(freqs)=groups
    return(cts[,1])

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
###### 1) mergedData=merged data from readMergeSave,
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
### make the count matrix
    ctsPR=cTabPR(clone,mergedData,correct=c.corr)
    ctsPR0=ctsPR-c.corr
### make the regression data
    datPR=cDfPR(ctsPR,peps=peptides, control=control)

### perform the regression
    #mhcMod <- glm(cts ~ clone*pep, data = datPR, family = poisson(link = "log"))
mhcMod <- glm.nb(cts ~ clone*pep, data = datPR)
    coefs=summary(mhcMod)$coef

    interact=grep(":",rownames(coefs),fixed=T)

          #best=interact[which(coefs[interact, "z value"]==max(coefs[interact, "z value"]))][1]
    pepOrd=interact[order(coefs[interact, "z value"],decreasing=T)]
    best=pepOrd[1]
    bestPep=sub("clonec+:pep","",rownames(coefs)[best],fixed=T)
    scnd=pepOrd[2]
    cMat <- matrix(rep(0,length(mhcMod$coef)), 1) ## initialize contrast matrix
    cMat[1,c(best,scnd)]=c(1,-1)   ### specify contrast top peptide to second
    pepComp<- glht(mhcMod, linfct=cMat)  ### fit
    compP=summary(pepComp)$test$pvalues

    bCts=ctsPR[which(peptides==bestPep),,drop=F]
    bRts=sort(bCts[,1]/apply(bCts,1,sum),decreasing=T)
    ctDisc=bRts[2]/bRts[1]
    names(ctDisc)=""
    # if running in screening mode, return a simplified result
 if(screen==T){
    ans=c("pep"=rownames(coefs)[best],prepStatsPR(summary(mhcMod)$coef[best,]),"coef2"=summary(mhcMod)$coef[scnd,1],"p2"=summary(mhcMod)$coef[scnd,4],"ctDiff"=ctDisc)
    #names(ansP)=c("pep","OR","LCB","UCB","pval","coef2","p2","ctDiff")

 }else{ ### if running in final mode, return a detailed result

     detail=dfResultPR(ctsPR0,coefs[interact,])

     if(printDetail){
     #write.xlsx(detail,file=outFile,sheetName=clone,col.names=F,row.names=F,overwrite=F)

     }

   res=c(rownames(coefs)[best],prepStatsPR(summary(mhcMod)$coef[best,]),vResultPR(ctsPR0,peptides))
     res[1]=sub("clonec+:pep","",res[1],fixed=T)
     names(res)[1]="peptide"
    ans=list(res,detail)}
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
    of=paste("screen",outF,"Cands.xlsx",sep="_")
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

    screen=sapply(goodClones,poisReg,mergedData=mergeDat,peptides=gps,control=cont,c.corr=1,screen=T)

    screenM=matrix(as.numeric(t(screen[-1,])),ncol=nrow(screen)-1)
    rownames(screenM)=screen[1,]

    colnames(screenM)=rownames(screen)[-1]

topPbs=which(rownames(screenM)%in%c("clonec+:pepE6","clonec+:pepE7","clonec+:pepL2") & as.numeric(screenM[,"OR"])>1 & p.adjust(as.numeric(screenM[,"pval"]))<0.001 & as.numeric(screenM[,"coef2"])>0 & as.numeric(screenM[,"p2"])>0.01)

topPbsF=goodClones[topPbs[order(as.numeric(screenM[topPbs,"ctDiff"])<.1,as.numeric(screenM[topPbs,"pval"]),decreasing=F)]]

 # fullRes=data.frame(t(sapply(topPbsF,poisReg,mergedData=mergeDat,peptides=gps,control="NP",c.corr=1,screen=F,printDetail=T,outFile=of)))

 fullRes=matrix(rep(NA,(length(topPbsF)+1)*(length(gps)+5)),nrow=length(topPbsF)+1)
 rownames(fullRes)=c("denominators",topPbsF)
 colnames(fullRes)=c("pep","OR","LCB","UCB","pval",gps)
 fullRes[1,-(1:5)]=as.character(readSums)
 for(cln in topPbsF){
     ans=poisReg(cln,mergedData=mergeDat,peptides=gps,control="NP",c.corr=1,screen=F,printDetail=T)
     fullRes[cln,]=ans[[1]]
     addWorksheet(wBook,cln)
 writeData(wBook,sheet=cln,x=ans[[2]],colNames=F,rowNames=F)
     }
    # wb=loadWorkbook(file=of)
   # sheets=getSheets(wb)
                                        # sumSheet=sheets[[1]]

    writeData(wBook,sheet="summary",x=fullRes,colNames=T,rowNames=T)
    saveWorkbook(wBook, of,overwrite=T)

   # write.xlsx(fullRes,file=of,col.names=T,row.names=T,sheetName="S",append=T)
return(fullRes)


}


###################### usage

# vdjFiles=list.files("RR_fest_Processed_data",full.name=T)
# files003Wk11=grep("1_003_wk11",vdjFiles,value=T)
# gp003Wk11=c("CEF","CEF","E6","E6","E6","E7","E7","E7","HIV","HIV","L2","L2","L2","NP","NP")
#
# out003=runRM(files003Wk11, gps=gp003Wk11,outF="Sample-1-003-Wk11")
