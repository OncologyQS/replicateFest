

# Functions for analysis of Manafest data in shiny app
# input is experiment without replicate
# the analysis is done using pairwise Fisher's test
#

#' @title runFisher
#' @description Runs Fisher's exact test for a pair of samples
#' @param pair a pair of samples to compare.
#' The reference sample is the second sample in pair if dir = T
#' (suggesting expansion compare to the reference and
#' corresponding OR > 1)
#' @param mergedData a list of data frames with read counts for each sample
#' @param nReadFilter a vector with the number of reads for each sample in the pair
#' @param dir a logical value indicating the direction of comparison
#' @param clones a vector with clones to analyze. If NULL, all clones
#' with the number of reads >= nReadFilter[1] are analyzed
#' @return a data frame with p-value, FDR, odds ratio, and confidence interval
#' @export
# Use filter for the number of reads (nReadFilter) if a set of clones to analyze is not provided
runFisher = function(pair, mergedData,
                     nReadFilter = c(10,0), dir = T,
                     clones = NULL)
{
#  browser()
  # calculate the number of reads for each sample
  totalReadCounts = sapply(mergedData, sum)
  # get sample IDs to compare
	sampID1 = pair[1]
	sampID2 = pair[2]
	# get clones to test
	if (dir) {
	  clonesToTest = setdiff(names(mergedData[[sampID1]])[which(mergedData[[sampID1]] >= nReadFilter[1])],"")
	}else{
		clonesToTest = union(names(mergedData[[sampID1]])[which(mergedData[[sampID1]] >= nReadFilter[1])],
			names(mergedData[[sampID2]])[which(mergedData[[sampID2]] >= nReadFilter[2])])
		clonesToTest = setdiff(clonesToTest,"")
	}

#print(length(clonesToTest))

# if set of clones is specified take intersection
	if(!is.null(clones)) clonesToTest = intersect(clonesToTest, clones)

	if (length(clonesToTest) == 0)
	{
		print('There are no clones to test')
		return(NULL)
	}

	dat = matrix(0,ncol = 2, nrow = length(clonesToTest))
	colnames(dat) = c(sampID1,sampID2)
	rownames(dat) = clonesToTest
	s = intersect(clonesToTest,names(mergedData[[sampID1]]))
	dat[s,sampID1] = mergedData[[sampID1]][s]
	s = intersect(clonesToTest,names(mergedData[[sampID2]]))
	dat[s,sampID2] = mergedData[[sampID2]][s]
#	print(head(dat))
#	print(dim(dat))

	res = matrix(nrow = length(clonesToTest), ncol = 4)
	rownames(res) = clonesToTest
#print(sampID1)
	for (clone in clonesToTest)
	{
#print(clone)
		counts1 = dat[clone,sampID1]
		counts2 = dat[clone,sampID2]
#		browser()
		tab = rbind(c(counts1,counts2),c(totalReadCounts[sampID1]-counts1, totalReadCounts[sampID2]-counts2))
		fres = fisher.test(tab)
#print(fres)
		res[clone,] = c(fres$p.value,fres$estimate,fres$conf.int)
	}
#print(dim(res))
	if(length(clonesToTest)>1)
	{
		res = cbind(res[,1],p.adjust(res[,1], method = 'BH'),res[,c(2:4)])
	}else{
		res = data.frame(matrix(c(res[,1],p.adjust(res[,1], method = 'BH'),res[,c(2:4)]), nrow = 1))
		rownames(res) = clonesToTest
	}
#print(dim(res))
	colnames(res) = c('p.val','FDR','odds.ratio','CI_low','CI_up')
	res = res[order(res[,'FDR']),]
	return(res)
}

#' createResTable
#' @description Creates a table with significant clones and conditions
#' @param res a table with results of Fisher's test
#' @param mergedData a list of data frames with read counts for each sample
#' @param orThr a threshold for odds ratio
#' @param fdrThr a threshold for FDR
#' @param saveCI a logical value indicating if confidence intervals should be saved
#' @param significanceTable a logical value indicating if a table with significant clones should be returned
#' @return a data frame with significant clones and the corresponding
#'  OR, FDR, counts, and percentages for all conditions
#' @export
# write output
# p-value, OR from fisher test + abundance + frequency
createResTable = function(res,mergedData,
                          orThr = 1,
                          fdrThr = 0.05,
                          percentThr = 0,
                          saveCI = T,
                          significanceTable = F)
{
  # get all clones for output
  clones = c()
	notNullComp = names(res)[!sapply(res,is.null)]
	# get all clones for output
	for (i in notNullComp){
#print(i)
		clones = union(clones,rownames(res[[i]])[which(as.numeric(res[[i]][,'odds.ratio']) >= orThr &
			as.numeric(res[[i]][,'CI_low']) > 1 & as.numeric(res[[i]][,'FDR']) < fdrThr)])
	}
	if(length(clones)==0){print('There is no significant clones'); return(NULL)}
	# create output matrices
	output_fdr = output_OR = output_CI = matrix(nrow = length(clones), ncol = length(res))
	colnames(output_fdr) = colnames(output_OR) =
	  colnames(output_CI) = names(res)
	rownames(output_fdr) = rownames(output_OR) =
	  rownames(output_CI) = clones
	for (i in notNullComp)
	{
		rows = intersect(rownames(res[[i]]),rownames(output_fdr))
		output_fdr[rows,i] = res[[i]][rows,'FDR']
		output_OR[rows,i] = round(res[[i]][rows,'odds.ratio'],3)
		if (saveCI) {output_CI[rows,i] = paste(round(res[[i]][rows,'CI_low'],3),'-',round(res[[i]][rows,'CI_up'],3), sep = '')}
	}
	colnames(output_fdr) = paste('FDR:',names(res))
	colnames(output_OR) = paste('OR:',names(res))
	colnames(output_CI) = paste('CI95%:',names(res))
	if(saveCI) {
		output_OR_CI = c()
		for(i in 1:ncol(output_OR))
		{
		  output_OR_CI = cbind(output_OR_CI, output_OR[,i],output_CI[,i])
		}
		colnames(output_OR_CI) = paste(rep(c('OR:','CI95%:'),length(res)),rep(names(res),each = 2))
	}else {output_OR_CI = output_OR}
	# calculate the number of significant comaprisons
	significant_comparisons = apply((as.numeric(output_fdr) < fdrThr &
	                                   output_OR > orThr),1,sum,na.rm = T)


	# get abundances and percentages
	output_counts_percent = getCountsPercent(clones, mergedData,
	                         samp = names(mergedData))
	# grep columns with percentage
	percCol = grep("percent",colnames(output_counts_percent), value = T)
	# get clones with maximum percentage higher than specified threshold
	output_counts_percent = output_counts_percent[apply(output_counts_percent[,percCol],1,max) >percentThr,]

	# update clones to output
	clones = rownames(output_counts_percent)

	outTab = data.frame(clone = clones,
	                    n_significant_comparisons =significant_comparisons[clones],
	                    output_fdr[clones,],
	                    output_OR_CI[clones,],
	                    output_counts_percent, check.names = F)
	# remove not significant clones
	outTab = outTab[which(outTab[,'n_significant_comparisons'] > 0),]

	# if significanceTable return a binary table clones vs conditions specifying which clone is significant in what condition
	if (significanceTable)
	{
		signTab = matrix(as.numeric(output_fdr[rownames(outTab),]) < fdrThr & output_OR[rownames(outTab),] > orThr,
			nrow = nrow(outTab), ncol = length(res),
			dimnames = list(rownames(outTab),names(res)))
		return(signTab)
	} else return(outTab)
}


#' getCountsPercent
#' @description Returns read counts and frequencies for selected clones and samples
#' @param clones a vector with clones to analyze
#' @param mergedData a list of data frames with read counts for each sample
#' @param samp a vector with sample IDs
#' @return a data frame with read counts and frequencies for selected clones and samples
#' @export
# returns read counts and frequencies for selected clones and samples
getCountsPercent = function(clones, mergedData, samp = names(mergedData))
{
  # calculate the total number of reads for each sample
  totalReadCountPerSample = sapply(mergedData, sum)

  output_counts = output_freq = matrix(0,nrow = length(clones), ncol = length(samp))
	rownames(output_counts) = rownames(output_freq) = clones
	colnames(output_counts) = colnames(output_freq) = samp
	for (i in samp)
	{
		rows = intersect(names(mergedData[[i]]),clones)
		output_counts[rows,i] = mergedData[[i]][rows]
		output_freq[rows,i] = (mergedData[[i]][rows]/totalReadCountPerSample[i])*100
	}
	colnames(output_counts) = paste(samp,'abundance', sep = '_')
	colnames(output_freq) = paste(samp,'percent', sep = '_')
	tab = cbind(output_counts,output_freq)
	return(tab)
}


# returns difference in percent for selected clones and samples
getDiff = function(clones, mergedData, samp = names(mergedData), refSamp)
{
  totalReadCountPerSample = sapply(mergedData, sum)
  output_freq = matrix(0,nrow = length(clones), ncol = length(samp)+1, dimnames = list(clones, c(samp,refSamp)))
	# get frequincies for all samples
	for (i in c(samp,refSamp))
	{
		rows = intersect(names(mergedData[[i]]),clones)
		output_freq[rows,i] = (mergedData[[i]][rows]/totalReadCountPerSample[i])*100
	}
	if (ncol(output_freq) == 2)
	{
		output_diff = data.frame(output_freq[,samp] - output_freq[,refSamp])
		colnames(output_diff) = paste0('DIFF:',samp)
	}else{
		output_diff = apply(output_freq[,samp],2,function(x) x-output_freq[,refSamp])
		colnames(output_diff) = paste0('DIFF:',colnames(output_diff))
	}
	return(output_diff)
}

#' @title getFreq
#' @description Returns frequencies in percent for selected clones and samples
#' @param clones a vector with clones, for which frequencies should be calculated
#' @param mergedData a list of data frames with read counts for each sample
#' @param samp a vector with sample IDs
#' @param colSuf a suffix for column names
#' @param minRead a threshold for the number of reads
#' @return a data frame with frequencies in percent for selected clones and samples
#' @export
#'
getFreq = function(clones, mergedData, samp = names(mergedData), colSuf = 'percent', minRead = 0)
{
  totalReadCountPerSample = sapply(mergedData, sum)
  output_freq = matrix(minRead,nrow = length(clones), ncol = length(samp))
	rownames(output_freq) = clones
	colnames(output_freq) = samp
	for (i in samp)
	{
		rows = intersect(names(mergedData[[i]]),clones)
		output_freq[rows,i] = mergedData[[i]][rows]
		output_freq[,i] = (output_freq[,i]/totalReadCountPerSample[i])*100
	}
	if (colSuf != '') colnames(output_freq) = paste(samp,colSuf, sep = '_')
	else colnames(output_freq) = samp
	return(data.frame(output_freq,check.names = F))
}

#' @title getFreqOrCount
#' @description
#' Returns frequencies or abundances in percent for selected clones and samples
#' @param clones a vector with clones, for which frequencies should be calculated
#' @param mergedData a list of data frames with read counts for each sample
#' @param samp a vector with sample IDs
#' @param colSuf a suffix for column names
#' @param minRead a threshold for the number of reads
#' @param returnFreq a logical value indicating if frequencies should be returned
#' @return a data frame with frequencies or abundances in percent for selected clones and samples
#' @export
getFreqOrCount = function(clones, mergedData,
                          samp = names(mergedData),
                          colSuf = 'percent',
                          minRead = 0,
                          returnFreq = T)
{
  totalReadCountPerSample = sapply(mergedData, sum)
  output_freq = matrix(minRead,nrow = length(clones), ncol = length(samp))
	rownames(output_freq) = clones
	colnames(output_freq) = samp
	for (i in samp)
	{
		rows = intersect(names(mergedData[[i]]),clones)
		output_freq[rows,i] = mergedData[[i]][rows]
		if (returnFreq) output_freq[,i] =
		  (output_freq[,i]/totalReadCountPerSample[i])*100
	}
	if (colSuf != '') colnames(output_freq) = paste(samp,colSuf, sep = '_')
	else colnames(output_freq) = samp
	return(data.frame(output_freq,check.names = F))
}

#' return fold change of frequencies for selected clones and samples
#' @param clones a vector of clones to calculate FC
#' @param mergedData a list of clones and
#' the corresponding read for all samples in the experiment
#' @param refSamp one sample from mergedData relative to which
#' the fold change should be calculated
#' @param samp a vector of sample names for which FC should be calculated
#' @param colSuf a string to be used at the beginning of columns with FC

getFC = function(clones, mergedData, refSamp,
                 samp = names(mergedData), colSuf = 'FC:')
{
  totalReadCountPerSample = sapply(mergedData, sum)
  samp = setdiff(samp,refSamp)
	freq = getFreq(clones, mergedData, union(refSamp,samp), colSuf = '', minRead = 1)
	# reference samples frequencies
	ref = freq[,refSamp]
	# divide to reference samples frequencies
	fc = round(sweep(data.frame(freq[,samp]),1, ref, "/"))
	colnames(fc) = paste0(colSuf,samp)
	rownames(fc) = clones
	return(fc)
}

# returns clones that unique for a selected condition
getUniqueClones = function(samp, mergedData, readCountThr = 0)
{
	res = names(mergedData[[samp]])[which(mergedData[[samp]]>= readCountThr)]
	for(i in setdiff(names(mergedData),samp)) { res = setdiff(res,names(mergedData[[i]])[which(mergedData[[i]]> readCountThr)])}
	return(res)
}

#' @title makeHeatmaps
#' @description Makes heatmaps of frequencies of each element
#' of input list of clones and samples
#' and plots fold change if refSamp is specified
#' @param listOfclones a list of clones to plot
#' @param mergedData a list of data frames with read counts for each sample
#' @param samp a vector with sample IDs
#' @param refSamp a reference sample ID
#' @param fileName a name of the output file
#' @param size a height and width of heatmap in the output PDF file
# make heatmaps of frequencies (if refSamp is not specified) of each element of input list of clones and samples and plot FC if refSamp is specified
makeHeatmaps = function(listOfclones, mergedData,
                        samp = names(mergedData), refSamp = NULL,
                        fileName = 'heatmap.pdf',size = 7)
{
  totalReadCountPerSample = sapply(mergedData, sum)
	if (length(listOfclones)==0)
	{
		print('There are no clones to plots')
		return;
	}
	require(gplots)
	pdf(fileName, width = size, height = size)
	samp = setdiff(samp, refSamp)
	for(i in 1:length(listOfclones))
	{
		clones = listOfclones[[i]]
		if(length(clones)< 2) next;
		plotData = getFreq(clones, mergedData, samp)
		if (!is.null(refSamp))
		{
			freqRef = getFreq(clones, mergedData, refSamp)
			freqRef[freqRef == 0] = min(plotData[plotData > 0]) - 1e-10
			plotData = sweep(plotData,1,unlist(freqRef),'/')
		}
#		par(oma = c(.1,.1,.1,3))
		heatmap(as.matrix(plotData), col = bluered(100),
		        labCol = lapply(strsplit(colnames(plotData),'_percent'),getElement,1),
			scale = 'row', cexCol = .5, cexRow = .3)
		title(main = names(listOfclones)[[i]], cex=0.2)
	}
	dev.off()
}

runSingleFisher = function(clone, pair, mergedData)
{
  totalReadCounts = sapply(mergedData, sum)
  samp1 = pair[1]
	samp2 = pair[2]
	counts1 = mergedData[[samp1]][clone]
	counts2 = mergedData[[samp2]][clone]
	if(is.na(counts2)) counts2 = 0
	tab = rbind(c(counts1,counts2),c(totalReadCounts[samp1]-counts1, totalReadCounts[samp2]-counts2))
	fres = fisher.test(tab)
	return(c(fres$p.value,fres$estimate,fres$conf.int))
}

#' getPositiveClones
#' @description
#' returns a vector with positive clones as names and conditions,
#' in which a clone is significant, as values.
#' For every clone, checks top two conditions with the highest number of reads
#' to see if a clone is uniquely expanded
#' @param analysisRes a list with results of Fisher's test
#' for pairwise comparisons
#' @param mergedData a list of data frames with read counts
#' for each sample
#' @param samp a vector with sample IDs to analyze
#' @param orThr a threshold for odds ratio
#' @param fdrThr a threshold for FDR
#' @param nReads a threshold for the number of reads
#' @param percentThr a threshold for the percentage of reads
#' @return a vector with positive clones as names and conditions,
#' in which a clone is significant, as values
#' @export
#'
getPositiveClones = function(analysisRes, mergedData,
                             samp = names(mergedData),
                             orThr = 1,
                             fdrThr = 0.05,
                             percentThr = 0)
{
  totalReadCounts = sapply(mergedData, sum)
  resTable = createResTable(analysisRes, mergedData, orThr = orThr,
				fdrThr=fdrThr, saveCI =F)
	if(is.null(resTable)) return(NULL)

	# find significant expansions in one condition
	clones = rownames(resTable)[which(resTable[,'n_significant_comparisons'] == 1)]
	#================
	# find condition in which it's expanded
	signTable = createResTable(analysisRes, mergedData, orThr = orThr,
				fdrThr=fdrThr, saveCI =F, significanceTable = T)
#	signTable = data.frame(signTable[clones,])
	signMatrix = matrix(signTable[clones,],
			nrow = length(clones),
			ncol = ncol(signTable),
			dimnames = list(clones,sapply(strsplit(colnames(signTable),'_vs_'), function(x)x[1])))

	# returns condition in which a clone is significant
	signCond = apply(signMatrix,1,function(x) colnames(signMatrix)[which(x)])
	#=====================
	# if we have only one comparison
	if (length(analysisRes) == 1)
	{
		# get name of the condition
		fdrCol = grep('FDR: ', colnames(resTable), fixed = T, value = T)
		n = unlist(strsplit(fdrCol,'FDR: '))[2]
		n = unlist(strsplit(n,'_vs_'))[1]
		return(data.frame(clone = clones,
		                  significant_condition = rep(n, length(clones))))
	}
	#=============================
	# check if they are unique by checking top two conditions with the highest number of reads
	freqMatrix = getFreqOrCount(clones,mergedData,samp,
	                            colSuf = '',
	                            returnFreq = T)

	#==============================
	# select clones with max percentage more than percentThr
	 freqMatrix = freqMatrix[apply(freqMatrix,1,max) > percentThr,]

	#===============================
	# compare with the second highest
	fishRes1 = getFisherForNclone(freqMatrix, rownames(freqMatrix),2,mergedData)
	# compare with the third highest
	fishRes2 = getFisherForNclone(freqMatrix, rownames(freqMatrix),3,mergedData)
	# combine results
	fishResComb = cbind(fishRes1[,'FDR'],fishRes2[,'FDR'], fishRes1[,'odds.ratio'],fishRes2[,'odds.ratio'])

	# select clones that have significant FDRs and OR higher than threshold
	# meaning that a clone is significant and unique expansion
	# also select clones that have NAs in FDR and OR, which means
	# that this clone appears in only one condition and
	# there is nothing to compare
	# check if there is a condition that is also significantly expanded
	fdrClones2 = apply(fishResComb,1,function(x) any(as.numeric(x[1:2])>fdrThr|as.numeric(x[3:4])<orThr))

	# select clones that have maximum frequency across conditions
	# more than the threshold
	posClones = names(fdrClones2)[which(!fdrClones2|is.na(fdrClones2))]
	if(length(posClones)>0)
	{
	# get corresponding condition
#		return(signCond[posClones])
	  return(data.frame(clone = posClones,
	                    significant_condition = signCond[posClones]))
	}else{return(NULL)}
}

#' @title getPositiveClonesFromTopConditions
#' @description
#' Selects clones that have significant FDRs and OR higher than threshold
#' meaning that a clone is significant and uniquely expanded
#' also select clones that have NAs in FDR and OR, which means that
#' this clone appears in only one condition and there is nothing to compare
#' checks if there is a condition that is also significantly expanded.
#' This function takes a table with FDRs, ORs, and condition in which
#' the clone is the most abundant
#' Column 1 and 2 are FDRs, 3 and 4 - ORs, and 5 is condition
#' It returns a vector with positive clones as names and conditions,
#' in which a clone is significant, as values
#' @param fisherResTable a table with FDRs, ORs, and condition
#' in which the clone is the most abundant.
#' This table is an output from compareWithOtherTopConditions function
#' @param orThr a threshold for odds ratio
#' @param fdrThr a threshold for FDR
#' @param percentThr a threshold for percentage
#' @return a vector with positive clones as names and conditions,
#' in which a clone is significant, as values
#' @export
getPositiveClonesFromTopConditions = function(fisherResTable,
                                              orThr = 1,
                                              fdrThr = 0.05,
                                              percentThr = 0, ...)
{
  # apply FDR and OR thresholds
  fdrClones2 = apply(fisherResTable,1,function(x) any(as.numeric(x[1:2])>fdrThr|as.numeric(x[3:4])<orThr))
  # find positive clones
  posClones = names(fdrClones2)[which(!fdrClones2|is.na(fdrClones2))]

  #=============================
  # check if they are unique by checking top two conditions with the highest number of reads
  freqMatrix = getFreqOrCount(posClones,
                              colSuf = '',
                              returnFreq = T, ...)

  #==============================
  # select clones with max percentage more than percentThr
  freqMatrix = freqMatrix[apply(freqMatrix,1,max) > percentThr,]
  posClones = rownames(freqMatrix)

  if(length(posClones)>0)
  {
    # return conditions of positive clones
    #		return(fisherResTable[posClones,'condition'])
    return(data.frame(clone = posClones,
                      significant_condition = fisherResTable[posClones,'condition']))

  }else{return(NULL)}
}

# runs Fisher's test for the nth clone with the highest frequency/the number of reads
getFisherForNclone = function(freq, clones, n = 2,mergedData)
{
  productiveReadCounts = sapply(mergedData, sum)
	fishRes = c()
	freq[freq==0] = NA
	for( i in clones)
	{
		r = unlist(freq[i,])
		r = r[order(r,decreasing = T)]
		if (!is.na(r[n]))
		  {
		    res = runSingleFisher(i,names(r)[c(1,n)],mergedData)
		    }else {
		    res = rep(NA,4)
		    }
		fishRes = rbind(fishRes,c(names(r)[1],res, n, names(r)[n]))
	}
	rownames(fishRes) = clones
	fishResMatrix = matrix(ncol = ncol(fishRes)+1, nrow = nrow(fishRes),
		dimnames = list(rownames(fishRes),c('condition','p.val','FDR','odds.ratio','CI_low','CI_up','nComp','condToComp')))
	fishResMatrix[,'FDR'] = p.adjust(fishRes[,2], method = 'BH')
	fishResMatrix[,setdiff(colnames(fishResMatrix),'FDR')] = fishRes
	return(fishResMatrix)
}

#
# creates output tables to be saved in Excel
createPosClonesOutput = function(posClones,
                                 mergedData,
                                 refSamp = NULL,
                                 replicates = FALSE)
{
  output = vector(mode = 'list')

  # get per sample summary of positive clones
  # a table with the number of positive clones per sample
  # and the sum of their frequencies
	output$condition_summary = getPerSampleSummary(posClones,
	                                               mergedData,
	                                               replicates)

	# get significant conditions
	signCond = unique(posClones$significant_condition)

	clones = posClones$clone

	# write clone-level summary
		if(!is.null(refSamp))
		{
		  # find replicates of refSamp in merged data
		  refSamp = grep(refSamp, names(mergedData), value = T)
      # a table with a per clone information
		  tab = data.frame(posClones)
	    # in case of several samples in the reference condition
	    # iterate through the samples and add FC for every sample
	    for( i in refSamp)
	    {
  		  # calculate FC relative to reference
  		  fc_ref = getFC(unique(clones),mergedData,i, signCond, "")
  			tab = data.frame(tab, getFreq(clones,mergedData,i),
  				sapply(1:nrow(posClones),
  				       function(j)
  				         fc_ref[posClones[j,"clone"],posClones[j,"significant_condition"]]),check.names = F)
  			colnames(tab)[ncol(tab)] = paste0('FC:', i)
	    }
		}else{
			tab = posClones
		}

	output$positive_clones_summary = tab

	# extended information per clone for all samples
	# if data is with replicates, remove duplications in clones and
	# keep only significant comparison, not sample
	if (replicates)
	{
	  posClones = posClones[which(!duplicated(posClones$clone)),c("clone","significant_comparison")]
	}
	# add per sample percentage
	tab = getCountsPercent(posClones$clone, mergedData)
#browser()

  output$positive_clones_all_data = data.frame(posClones,
                                               tab,
                                               check.names = F)
	return(output)
}

#' getPerSampleSummary
#' returns per sample summary of positive clones a table with the
#' number of positive clones per sample and the sum of their frequencies
#' @params posClones a data frame with positive clones in the first column
#' and their conditions to be summarized in the second.
#' In the replication version, these should be samples
#' @params mergedData a list of data frames with read counts for each sample
#' @params replicates a logical value indicating if there are replicates
getPerSampleSummary = function(posClones, mergedData, replicates = FALSE)
{
  library(dplyr)
  clones = unique(posClones[,1])
  # write peptide summary of positive clones
  # get frequencies of positive clones across all samples
  freqMatrix = getFreq(clones,mergedData,names(mergedData), colSuf = '')

  # make a list of conditions with positive clones in each condition
  peptLevelList = tapply(posClones[,1],posClones[,2],
                         FUN = function(x)return(x))
  peptideTab = data.frame(condition = names(peptLevelList),
      n_positive_clones = sapply(peptLevelList,length),
      sum_freq = sapply(names(peptLevelList),
                          function(x) sum(freqMatrix[peptLevelList[[x]],x])))
  # if there are replicates, further summarize frequency by conditions
  if (replicates)
  {
    # extract conditions from the file names
    sampAnnot = splitFileName(names(mergedData))
    # add peptide name to peptideTab to be able to summarize by it
    tab = merge(sampAnnot[,c("file","condition")],peptideTab,
                by.y = "condition", by.x = "file")
    # summarize
    tab = tab %>% group_by(condition, n_positive_clones) %>% summarise(sum_freq = sum(sum_freq))
    # convert to data.frame
    peptideTab = as.data.frame(tab)
  }

  return(peptideTab)
}
#' getFreqThreshold
# calculate frequency threshold using the number of cells and probability
# (1-((1-p)^(1/n)))
# where n=number of cells per well and p=selected probability
getFreqThreshold = function(n, p)
{
	return (1-((1-p)^(1/n)))
}

#' @title compareWithOtherTopConditions
#' @description
#' selects clones that pass the nReads threshold and compares top conditions
#' with the second and third top conditions
#' and returns a table with FDRs and ORs for those comparisons
#' that will be used to select positive clones using FDR and OR threshold later
#' there is no comparison to the reference condition
#' @param mergedData a list of data frames with read counts for each sample
#' @param sampForAnalysis a vector with sample IDs
#' @param nReads a threshold for the number of reads. Default is 10
#' @param clones an optional vector with clones to analyze.
#' Default is NULL, which means that all clones with the number of reads more than nReads will be analyzed
#' @return a table with FDRs and ORs for comparisons with the second and third top conditions
compareWithOtherTopConditions = function(mergedData,
                                         sampForAnalysis,
                                         nReads = 10,
                                         clones = NULL)
{
  productiveReadCounts = sapply(mergedData, sum)
  # combine clones from all conditions that have the number of clones more than nReads
	allClones = sapply(mergedData[sampForAnalysis],
	                   function(x) setdiff(names(x)[which(x >= nReads)],""))
	clonesToTest = c()
	for(i in names(allClones))
	{
		clonesToTest = union(clonesToTest, allClones[[i]])
	}

# if set of clones is specified take intersection
	if(!is.null(clones)) clonesToTest = intersect(clonesToTest, clones)

	if (length(clonesToTest) == 0)
	{
		print('There are no clones to test')
		return(NULL)
	}

	# get frequency matrix
	freqMatrix = getFreqOrCount(clonesToTest,mergedData,sampForAnalysis, colSuf = '', returnFreq = T)
	#print(countMatrix)
	# compare with the second highest
	fishRes1 = getFisherForNclone(freqMatrix, rownames(freqMatrix),2,mergedData)
	# compare with the third highest
	fishRes2 = getFisherForNclone(freqMatrix, rownames(freqMatrix),3,mergedData)
	# combine results
	fishResComb = cbind(FDR1 = fishRes1[,'FDR'], FDR2 = fishRes2[,'FDR'], OR1 = fishRes1[,'odds.ratio'], OR2 = fishRes2[,'odds.ratio'], condition = fishRes1[,'condition'])

	return(fishResComb)
}



#===========
# wrapper for running the full analysis using Fisher's test
# from reading files to output all results
#################
#' @title runExperimentFisher
#' @description Reads in files with TCR repertoires from
#' a FEST experiment without replicate samples. It
#' performs pairwise Fisher's test find expanded clones
#' comparing to a reference samples.
#' It also compares top conditions to find unique expansions.
#' The results are return and saved in an Excel file.
#' @param files a list of file names with read counts
#' @param refSamp a reference sample ID
#' @param nReads a threshold for the number of reads
#' @param fdrThr a threshold for FDR
#' @param orThr a threshold for OR
#' @param percentThr a threshold for percentage
#' @param excludeSamp a sample ID to exclude from analysis
#' @param compareToRef a logical value indicating if the comparison to the reference sample should be performed
#' @param outputFile a name of the output file
#' @param saveToFile a logical value indicating if the output should be saved to a file
#' @return a list with output tables
#' @export

runExperimentFisher=function(files,
                             refSamp,
                             nReads = 50,
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
                            c(excludeSamp,refSamp))
  productiveReadCounts = sapply(mergedData, sum)

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
    # create comparing pairs (to refSamp)
    compPairs = cbind(sampForAnalysis,rep(refSamp,
                                          length(sampForAnalysis)))
    # run pair-wise Fisher's test
#browser()
      fisherRes = apply(compPairs,1,runFisher,mergedData,
                          clones = clonesToTest,
                          nReadFilter = c(nReads,0))
    if (!is.null(fisherRes))
    {
     # add names of compared conditions
      names(fisherRes) = apply(compPairs,1,paste,collapse = '_vs_')
      # select positive clones with specified thresholds
      posClones = getPositiveClones(fisherRes, mergedData,
                                    samp = sampForAnalysis,
                                    orThr = orThr,
                                    fdrThr=fdrThr,
                                    percentThr = percentThr)
    }	else{
      print('There are no clones to analyze. Try to reduce confidence or the number of templates 1')
      return(NULL)
    }


  }else{ # if there is no comparison to ref sample

    fisherRes = compareWithOtherTopConditions(mergedData,
                                              sampForAnalysis,
                                              nReads = nReads,
                                              clones = NULL)

    # select positive clones with specified thresholds
    posClones = getPositiveClonesFromTopConditions(fisherRes,
                                                   orThr = orThr,
                                                   fdrThr = fdrThr,
                                                   percentThr = percentThr,
                                                   mergedData,
                                                   sampForAnalysis)
  }


  # create  output
  tablesToXls = createPosClonesOutput(posClones, mergedData,
                                      refSamp)
  #===================
  # add the ref_comparison_only sheet
  #===================
  if(compareToRef) #if there is comparison to the ref sample
  {
    # create a table with results
    resTable = createResTable(fisherRes,mergedData,
                              orThr = orThr,
                              fdrThr = fdrThr, saveCI = F)
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
  s = names(productiveReadCounts)
  param = c('Reference sample',
            'Excluded samples','Compare to reference',
            'n template threshold','FDR threshold',
            'OR threshold','percent threshold',
            'Nucleotide level analysis',
            'n analyzed samples',
            paste(s, 'n templates',sep = '_'))
  value = c(toString(refSamp),
            paste(excludeSamp, collapse = ', '),
            compareToRef, nReads,fdrThr,
            orThr, percentThr,
            FALSE,
            length(s), productiveReadCounts[s])

  tablesToXls$parameters = data.frame(param, value)

  # save into Excel file
  if(saveToFile)
  {
    saveResults(tablesToXls,outputFile)
  } else return(tablesToXls)

}
