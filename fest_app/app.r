# 2017-11-14 (v6)
# add fields for the number of cells and probability to calculate baseline frequency
# to filter out clones with low frequency
# this version uses only clones that appear in baseline above calculated frequency
# previously, we used 10 reads as a threshold and do not restrict to baseline clones

# 2018-01-23 (v7)
# Added a threshold for the number of reads into interface
# Allowed user to ignore baseline threshold
# Make only one heatmap with all positive clones
# Output all significant clones relative to the reference even there are no positive clones

# 2018-03-07 (v7.1)
# Trying to solve out of memory issue when saving large table to Excel
# Solved by using WriteXLS and saving all results at once
# fixed problem with using AA level data when selecting clones to analyze with Use nucleotide level option selected

# 2018-04-18 (v8)
# Added a tab with user's manual to UI
# Added baseline percent and the corresponding number of templates in baseline to output parameter sheet

# 2018-06-07 (v9)
# switch to the parser from the tcR package to eliminate sequencing error on 5' end of the nt sequence

# 2018-07-20 (v10)
# fixed reading function. It doesn't save mergedData and ntData to the disc;
# it returns a list with two element - mergedData and ntData
# fixed creating Excel file with no significant clones

# 2019-06-04
# fixed output parameters issue
# add the OR threshold when compare to all other conditions (previously, it was only the FDR threshold)

# 2019-06-14 (v12)
# improved an input process. Now a user can select what file format to upload: Adaptive, VDJtools, or the previously saved R object
# add excluded samples to output
# remove nRead requirement when compare with other wells

# 2020-05-03 (v13)
# switch to immunarch: can't install on shiny server
# make comparison with reference as an option (checkbox)
# if off, then each condition will be compared with top two other conditions

# 2023-12-01
# switched to immunarch

# 2025-02-05 (v14)
# Split loading the data and analysis to the different tabs
# added percentage threshold to the Fisher's test version
# fixed output to heatmap when there are no or one clone
# added all functions to the replicateFest package

# 2025-05-06 (v15)
# made load data, run analysis, and save the results to the different tabs
# added analysis with replicates
# increased the default for the number of reads
# removed the baseline option
# switched to openxlsx for saving results
# made a list with analysis results and parameters that were used for the analysis
# made custom separator for file names to separate replicates

# TODO:
# - add custom separator into interface
#   Couldn't easily add, because that the getPerSampleSummary function
#   also requires separator, but it's not straightforward to pass it
#   Probably, need to save that as another element in analysisRes object

# - In cases where there are not replicates,
# can you add a filter to only show clones that are detected
# in X% of other wells? For example, sometimes when there are
# a lot of clones that are only found in one well, I implement
# a threshold of 20%, meaning the clone should be detected
# in at least 20% of wells. detected would mean that there is
# a non-zero abundance

# - what do you think about removing the “compare to reference* button
# in favor of  testing reference==’none’?
# (Maybe want to reinitialize when new data is loaded
# to avoid errors)

# - A message saying “Analysis in progress” would help,
# to be replaced with “Analysis is done” , maybe also a
# note about expected wait

#==============================
# User interface
#==============================
ui <- fluidPage(
  #  title = "FEST data analysis",
  headerPanel("FEST data analysis"),
  # add description of the app and steps
  tags$div(
    tags$p("The FEST (Functional Expansion of Specific T-cells)
           application for analysis of TCR sequencing data of short-term,
           peptide-stimulated cultures to identify antigen-specific
           clonotypic amplifications. The app loads TCR sequencing data
           in the tab-delimited format, performs the analysis, and
           visualizes and saves results. For analysis, we use only
           productive clones and summarize template counts for
           nucleotide sequences that translated into the same amino acid
           sequence."),
    tags$h3("Steps"),
    tags$ol(
      tags$li(HTML("<b>Load data</b>: upload FEST files or a previously saved R
              object with data")),
      tags$li(HTML("After the data is successfully loaded,
              select the <b>Run analysis</b>
                   tab to run the analysis.
                   <br><b>Important!</b> Please make sure to specify
                   whether the input data has replicates or not.
                   If the input data has replicates, the conditions
                   will be extracted from the input file names.
                   To be able to correctly extract the conditions, file names should
                   follow the format 'sampleID_condition_replicate.ext'.
                   E.g. 'sample1_HIV_1.tsv'.
                   <br>If the input data has replicates,
                   the negative binomial model will be used in the analysis;
                   if the input data does not have replicates,
                   the Fisher's exact test will be used.")),
      tags$li(HTML("After the analysis is done,
              select the <b>Save results</b> tab
                   to save the results.")),
    ),
    tags$p(HTML("For more details, please refer to the
                <b>User's manual</b> tab"))
  ),

  # layout with tabs
  tabsetPanel(
    #============================
    # tab for loading the data
    tabPanel("Load data",
             sidebarLayout(
               # the left side panel
               sidebarPanel(
                 # select the format of input files
                 selectInput('inputFiles', 'Select an input format',
                             choices = c('FEST files','an R object with data'),
                             selected = NULL, multiple = FALSE,selectize = TRUE,
                             width = NULL, size = NULL),

                 # if a previously saved R object with the input data will be uploaded
                 conditionalPanel(
                   condition = "input.inputFiles == 'an R object with data'",
                   fileInput('inputObj', 'Upload an R object with data (.rda)',
                             multiple = FALSE,
                             accept=c('rda', '.rda'))
                 ),
                 # show the fileInput control depending on the selected values of inputFiles
                 # if files with raw data will be uploaded
                 conditionalPanel(
                   condition = "input.inputFiles == 'FEST files'",
                   fileInput('sourceFiles', 'Upload FEST files (.tsv)',
                             multiple = TRUE,
                             accept=c('text/tsv', 'text/comma-separated-values,text/plain', '.tsv')),
                   downloadButton('saveInputObj', 'Save Input Object')
                 ), # end conditionalPanel
              ), # end sidebarPanel

              # the main panel
              mainPanel(
                htmlOutput("message_load")
              )

          ) # end sidebarLayout
    ), # end tabPanel with loading data

    #==============================
    # tab for analysis
    #==============================
  	tabPanel("Run analysis",
  	 sidebarLayout(
  	   # the left side panel
  		sidebarPanel(
  		  # check box that specifies the input data format
    		checkboxInput('replicates','Analyze with replicates',
    		              value = FALSE),

    		# if analysis with replicates, then specify separator
    		# conditionalPanel(
    		#   condition = "replicates == TRUE",
    		#   tags$span("Specify a separator:"),
    		#   textInput('separator', NULL,
    		#             value = "_", width = "3em",
    		#             placeholder = NULL),
    		# ),


    		# check box that controls the type of analysis - if the comparison with a reference should be performed or not
    		shinyjs::useShinyjs(),
    		checkboxInput('compareToRef','Compare to reference', value = TRUE),

    		  		# specify if analysis should be performed on nucleotide level data
    		checkboxInput('nuctleotideFlag','Use nucleotide level data'),

      	# drop down menu to select a reference
    		selectInput('refSamp', 'Select a reference', choices = list('None'), selected = NULL, multiple = FALSE,
    			selectize = TRUE, width = NULL, size = NULL),


    		selectInput('excludeSamp', 'Select conditions to exclude', choices = list('None'), selected = NULL, multiple = TRUE,
    			selectize = TRUE, width = NULL, size = NULL),
    		textInput('nReads', 'Specify the minimal number of templates (increase in case of large files)',
    		          value = "50", width = NULL, placeholder = NULL),

   		  actionButton('runAnalysis', 'Run Analysis', width = '60%'),

  		),

  		# the main panel
  	  mainPanel(
  		 # textOutput('contents'),
  		  htmlOutput("message_analysis")
  		)
  	 )
	 ),# end tabPanel Analysis

    #=============================
    # tab for saving results of the analysis
    #=============================
    tabPanel("Save results",
         sidebarLayout(
           # the left side panel
           sidebarPanel(
             tags$p(HTML("<b>Specify thresholds:</b>")),
             numericInput('fdrThr', 'for FDR ', value = 0.05,
                          max = 1, step= 0.01),
             textInput('orThr', 'for odds ratio', value = "5",
                       width = NULL, placeholder = NULL),
             # keep clones that have maximal abundance above that values
             textInput('percentThr', 'for clone abundance (in percent)',
                       value = "0"),
             # a threshold for percent of wells where a clone has non-zero abundance
             textInput('wellsThr', 'for wells with non-zero abundance (in percent)',
                       value = "0"),
             downloadButton('saveResults', 'Download Results'),
             downloadButton('saveHeatmaps', 'Create heatmaps')
           ), # end sidebarPanel

           # the main panel
           mainPanel(
             htmlOutput("save_results")
           )

         ) # end sidebarLayout
    ), # end tabPanel with saving data


    #=============================
    # the tab with manual
    #=============================
	 tabPanel("User\'s manual",
	          includeHTML('userManual.html'))

 ) # end tabsetPanel
)# end of the UI


#====================
# Server
#====================
server <- function(input, output,session) {
  #rm(list = ls())
  # increase file upload limit to 100M
  options(shiny.maxRequestSize=100*1024^2, java.parameters = "-Xmx8000m")
  library(shiny)
  library(tools)
  library(gplots)
  library("Matrix")
  library(DT)
  library(kableExtra)
  library(dplyr)
  library(shinyjs)
  if (!require(openxlsx)) install.packages("openxlsx")
  if(!require(immunarch)) install.packages("immunarch")
  if (!require(devtools)) install.packages("devtools")
  if (!require(contrast)) install.packages("contrast")
  if (!require(multcomp)) install.packages("multcomp")
  if (!require(replicateFest)) devtools::install_github("OncologyQS/replicateFest")
  library(replicateFest)

#  source("R/manafest_shiny_functions.r")
#  source("R/repManFunctions.r")

  #===================
  # read input files
  observeEvent(input$sourceFiles,{
    output$message_load = renderUI({
    # check if there is file to analyze
    if (length(input$sourceFiles) == 0)
    {
      HTML('There are no files to analyze. Please upload files')
    }else if(length(input$sourceFiles$name) == 1){
      HTML('There should be at least two files to analyze')
    }else if(length(input$sourceFiles$name) > 1){
      # read loaded files in and save into inputData object

      # if there are previously loaded data, remove it
      if (exists('mergedData', envir = .GlobalEnv))
      {
        rm(list = c('mergedData','productiveReadCounts','ntData'), envir = .GlobalEnv)
        updateSelectInput(session, "refSamp", choices=list('None'))
#        updateSelectInput(session, "baselineSamp", choices=list('None'))
        updateSelectInput(session, "excludeSamp", choices=list('None'))
      }
      # specify input file format
      if (input$inputFiles == 'FEST files') fileFormat = 'fest'
      # extract path to a folder with input file to pass into reading function
      # supplied source names of loaded files to be used as names for data objects
      #  and read in input files
      res = readMergeSave(files = dirname(input$sourceFiles$datapath[1]),
                          filenames = unlist(input$sourceFiles$name))

      if(!is.null(res))
      {
        assign('mergedData', res$mergedData, envir = .GlobalEnv)
        assign('ntData', res$ntData, envir = .GlobalEnv)
        assign('productiveReadCounts', sapply(res$mergedData,sum), envir = .GlobalEnv)

        if(file.exists('inputData.rda')) file.remove('inputData.rda')
        updateSelectInput(session, "refSamp", choices=c('None',names(mergedData)))
        updateSelectInput(session, "excludeSamp", choices=names(mergedData))

        HTML(c('Uploaded files:<br/>',paste(names(mergedData), collapse = '<br/>')))
      }else{
        HTML('Error in reading files')
      }
    }

  })
  })

  #===================
  # load RDA file with the input data
  observeEvent(input$inputObj,{
    output$message_load = renderUI({
      # check if there is a file to analyze
      if (length(input$inputObj) == 0)
      {
        HTML('There are no files to analyze. Please upload files')
      }else if(length(input$inputObj$name) > 0 )
      {
        # if there are previously loaded data, remove it
        if (exists('mergedData', envir = .GlobalEnv))
        {
          rm(list = c('mergedData','productiveReadCounts','ntData'), envir = .GlobalEnv)
          updateSelectInput(session, "refSamp", choices=list('None'))
#          updateSelectInput(session, "baselineSamp", choices=list('None'))
          updateSelectInput(session, "excludeSamp", choices=list('None'))
        }
        if (tolower(file_ext(input$inputObj$name)) == 'rda'){
          # upload an existing input object
          load(input$inputObj$datapath, envir = .GlobalEnv)
          if (exists('mergedData', envir = .GlobalEnv))
          {
            updateSelectInput(session, "refSamp", choices=c('None',names(mergedData)))
#            updateSelectInput(session, "baselineSamp", choices=c('None',names(mergedData)))
            updateSelectInput(session, "excludeSamp", choices=names(mergedData))
            HTML(c('Uploaded samples:<br/>',paste(names(mergedData), collapse = '<br/>')))
          }else{ HTML('There are no input data to analyze. Please check the input file')}
        }else{
          HTML('This is not an .rda file. Please check the input file')
        }
      }
    })
  })

  #===================
  # save the inputData object with uploaded files for further analysis when the Save Input Object button is clicked
  output$saveInputObj <- downloadHandler(
    filename='inputData.rda',
    content=function(file){
      #		load('inputData.rda')
      if (exists('mergedData', envir = .GlobalEnv)) {
        save(mergedData,productiveReadCounts,ntData,file=file)
      }else{
        print('No data to save')
        output$message_load = renderText('No data to save')
        emptyObj = NULL
        save(emptyObj,file=file)
      }
    })

   #==================
  # an observerEvent for checking the The input with replicates checkbox
  observeEvent(input$replicates, {
    if(exists('mergedData', envir = .GlobalEnv))
    {

      if (input$replicates == TRUE){
        # if the input data has replicates, then the reference sample should be selected
        # update conditions
        # check if there are objects to run analysis
         # extract conditions from the file names
          sampAnnot = splitFileName(names(mergedData))
        # check if there are conditions
        if (all(is.na(sampAnnot$condition)))
        {
          output$message_analysis = renderText('There are no conditions to analyze.
                                                 Please check the input file names')
        }else{
          # output a table with sample annotations
          output$message_analysis = renderTable({
            # Create a table using kable
            kable(sampAnnot, format = "html") %>%
              kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
          }, sanitize.text.function = function(x) x)

          # update inputs in the dropdowns
          updateSelectInput(session, "refSamp",
                            choices=c('None',sampAnnot$condition))
          updateSelectInput(session, "excludeSamp", choices=sampAnnot$condition)

        }
        shinyjs::disable('compareToRef')

      # if The input with replicates is unchecked
      }else{
        # if the input data does not have replicates, then the reference sample should be selected
        # update conditions
        updateSelectInput(session, "refSamp", choices=c('None',names(mergedData)))
        updateSelectInput(session, "excludeSamp", choices=names(mergedData))
        shinyjs::enable('compareToRef')
        # output sample names in the main panel
        output$message_analysis = renderText(
          HTML(c('Uploaded samples:<br/>',paste(names(mergedData),
                                                collapse = '<br/>'))))

      }
    }
  } )

  #================
  # run analysis with the Run Analysis button is clicked
  observeEvent(input$runAnalysis,{
    output$message_analysis = renderText('Analysis is running...')
    # remove results of the previous analysis
    if(exists('analysisRes', envir = .GlobalEnv)) rm('analysisRes', envir = .GlobalEnv)
    # check if there are objects to run analysis
    if(exists('mergedData', envir = .GlobalEnv) && exists('ntData', envir = .GlobalEnv))
    {
      # run analysis on aa level data
      # load object with input data
      if (!input$nuctleotideFlag){
        sampNames = names(mergedData)
        obj = mergedData
      }else{ # run analysis on nucleotide level data
        #if 'Use nucleotide level data' checkbox is selected
        sampNames = names(ntData)
        obj = ntData
      }

      # create an object with results
      # it's a list with two elements:
      # one is the analysis output
      # and the other one is the parameters of the analysis
      analysisRes = list(res = NULL,
                         params = reactiveValuesToList(input))
#browser()
      #================================
      ## analysis without replicates
      #================================
      if (input$replicates == FALSE)
      {
        # list samples to analyze that excludes reference and other
        # samples that should be excluded from the analysis
        sampForAnalysis = setdiff(sampNames,
                                  c(input$excludeSamp,input$refSamp))

        # select clones to test
        clonesToTest = NULL
        # if the comparison to reference should be included in the analysis
        if(input$compareToRef)
        {
          if (input$refSamp == 'None')
          {
            output$message_analysis = renderText('There is no reference. Please select a reference sample')
          }else{

            # create comparing pairs (to refSamp)
            compPairs = cbind(sampForAnalysis,rep(input$refSamp,length(sampForAnalysis)))

            # run Fisher's test
            if(nrow(compPairs) == 1)
            {
              analysisRes$res = list(runFisher(compPairs,obj, clones = clonesToTest, nReadFilter = c(as.numeric(input$nReads),0)))
            }else{
              analysisRes$res = apply(compPairs,1,runFisher,obj, clones = clonesToTest, nReadFilter = c(as.numeric(input$nReads),0))
            }
            #			browser()
            if (!is.null(analysisRes$res))
            {
              names(analysisRes$res) = apply(compPairs,1,paste,collapse = '_vs_')
           }
          }
        }else{
          # if there is no comparison with the reference, then compare only within conditions
#          if(input$ignoreBaseline) clonesToTest = NULL
          # take only clones that have the number of reads more then nReads and compare with top second and top third conditions
          analysisRes$res = compareWithOtherTopConditions(obj,
                                                    productiveReadCounts,
                                                    sampForAnalysis,
                                                    nReads = as.numeric(input$nReads),
                                                    clones = clonesToTest)
       }
      } else {# end of analysis without replicates


      #================================
      ## analysis with replicates
      #================================
        # extract conditions from the file names
        sampAnnot = splitFileName(names(obj))
        # add that to the analysisRes object to be used in generating output
        analysisRes$sampAnnot = sampAnnot

        # check if there are conditions
        if (all(is.na(sampAnnot$condition)))
        {
          output$message_analysis = renderText('There are no conditions to analyze.
                                               Please check the input file names')
        }else{
          # output a table with sample annotations
          output$message_analysis = renderTable({
            # Create a table using kable
            kable(sampAnnot, format = "html") %>%
              kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))
          }, sanitize.text.function = function(x) x)

            # update inputs in the dropdowns
            updateSelectInput(session, "refSamp", choices=c('None',sampAnnot$condition))
#            updateSelectInput(session, "baselineSamp", choices=c('None',sampAnnot$condition))
            updateSelectInput(session, "excludeSamp", choices=sampAnnot$condition)

            # get clones to test
            clonesToTest = getClonesToTest(obj, nReads = as.numeric(input$nReads))

            # check if there are enought clones to analyze
            if (length(clonesToTest)<1)
            {
              output$message_analysis = renderText('There are no clones to analyze. Try to reduce confidence or the number of templates 1')
              return()
            }
            # run the analysis for selected clones
            analysisRes$res = fitModelSet(clonesToTest, obj,
                                     sampAnnot$condition,
                                     excludeCond = input$excludeSamp,
                                     refSamp=input$refSamp,
                                     c.corr=1)
            rownames(analysisRes$res) = analysisRes$res$clone

      }
    }# end of analysis with replicates
     # check if there are any results to save
    if (!is.null(analysisRes$res))
    {
        output$message_analysis = renderText('Analysis is done. Go to the Save results tab')
        assign('analysisRes',analysisRes, envir = .GlobalEnv)
    }	else{
        output$message_analysis = renderText('There are no clones to analyze. Try to reduce the number of templates')
    }


    } else{
      output$message_analysis = renderText('There are no data to analyze. Please load files')
    }
  })


  # save results with selected thresholds when the Download Results button is clicked
  output$saveResults <- downloadHandler(
    filename=function() {
      paste0('analysisRes_',Sys.Date(),'.xlsx')
    },
    contentType = "text/xlsx",
    content=function(file){
      if (exists('analysisRes', envir = .GlobalEnv))
      {
        # check if analysis was done on aa or nt level
        if (input$nuctleotideFlag) obj = ntData else obj = mergedData
        sampForAnalysis = setdiff(names(obj), c(input$excludeSamp,input$refSamp))
        #======================
        # get positive clones
        posClones = NULL
        if(input$compareToRef) #if there is comparison to the ref sample
        {

         # if analyzed data with replicate
          if(input$replicates)
          {
            posClones = getPositiveClonesReplicates(analysisRes$res,
                                                    obj,
                                                    refSamp = analysisRes$params$refSamp,
                                                    samp = sampForAnalysis,
                                                    excludeCond = analysisRes$params$excludeSamp,
                                                    orThr = as.numeric(input$orThr),
                                                    fdrThr = as.numeric(input$fdrThr),
                                                    percentThr = as.numeric(input$percentThr))
            }else{ # if analyzed data without replicates
              posClones = getPositiveClones(analysisRes$res, obj,
                                        samp = sampForAnalysis,
                                        orThr = as.numeric(input$orThr),
                                        fdrThr = as.numeric(input$fdrThr),
                                        percentThr = as.numeric(input$percentThr))
            }
        }else{ # if there is no comparison to ref sample
          posClones = getPositiveClonesFromTopConditions(analysisRes$res,
                                                         orThr = as.numeric(input$orThr),
                                                         fdrThr=as.numeric(input$fdrThr),
                                                         percentThr = as.numeric(input$percentThr),
                                                         mergedData = obj,
                                                         samp = sampForAnalysis)
        }

        #===============================
        # change "None" to NULL for reference sample
        refSamp = analysisRes$params$refSamp
        if(input$refSamp == 'None') refSamp = NULL

        #========================
        # create object with results to write to Excel
        #=======================
        tablesToXls = vector(mode = 'list')
        # if there is no positive clones
        if (nrow(posClones)==0)
        {
          output$save_results = renderText('There are no positive clones. Try to adjust thresholds')
          # add a sheet to the output with significant clones comparing to reference
          #clones = rownames(resTable)[which(resTable[,'significant_comparisons'] == 1)]
          tablesToXls$summary = data.frame('There are no positive clones', row.names = NULL, check.names = F)
        }else{
          # if there are positive clones, save them in to Excel file
          # create table with results
          tablesToXls = createPosClonesOutput(posClones,
                                              obj,
                                              analysisRes$params$refSamp,
                                              replicates = analysisRes$params$replicates)
        }
        #===================
        # add the ref_comparison_only sheet
        #===================
        if(analysisRes$params$compareToRef) #if there is comparison to the ref sample
        {
          # if analyzing data with replicate
          if(analysisRes$params$replicates)
          {
            resTable = createResTableReplicates(analysisRes$res,
                                                obj,
                                      orThr = as.numeric(input$orThr),
                                      fdrThr = as.numeric(input$fdrThr),
                                      percentThr = as.numeric(input$percentThr))

          }else{ # if without replicates
          # create a table with results
            resTable = createResTable(analysisRes$res,obj,
                                    orThr = as.numeric(input$orThr),
                                    fdrThr = as.numeric(input$fdrThr),
                                    percentThr = as.numeric(input$percentThr),
                                    saveCI = F)
          }
           #===========================
           # check if there are no clones to save
           if (is.null(resTable))
             tablesToXls$ref_comparison_only =
               data.frame(res = 'There are no significant clones')
           # if there are clones
           if(nrow(resTable) > 0)
             tablesToXls$ref_comparison_only = resTable
          }else{
            tablesToXls$ref_comparison_only =
              data.frame(res = 'There are no significant clones')
          }
        #============
        # add a sheet with parameters
        #============
        s = names(productiveReadCounts)
        param = c("Data with replicates",
                  'Reference sample',
                  'Excluded samples',
                  'Compare to reference',
                  'n template threshold',
                  'FDR threshold',
                  'OR threshold',
                  'percent threshold',
                  'Nucleotide level analysis',
                  'n samples',
                  paste(s, 'n templates',sep = '_'))
        value = c(analysisRes$params$replicates,
                  analysisRes$params$refSamp,
                  paste(analysisRes$params$excludeSamp, collapse = ', '),
                  analysisRes$params$compareToRef,
                  analysisRes$params$nReads,
                  input$fdrThr,
                  input$orThr,
                  input$percentThr,
                  analysisRes$params$nuctleotideFlag,
                  length(s), productiveReadCounts[s])

        tablesToXls$parameters = data.frame(param, value)

        #========
        # save results to xlsx
        saveResults(tablesToXls,file)
        output$save_results = renderText('The results are saved')
      }else{

        output$save_results = renderText('Please click the Run Analysis button')

      }
    })

  # save heatmaps with selected thresholds when the Download Heatmaps button is clicked
  output$saveHeatmaps <- downloadHandler(
    filename=function(){
      paste0('heatmap_',Sys.Date(),'.pdf')},
    contentType = "image/pdf",
    content=function(file){
      if (exists('analysisRes', envir = .GlobalEnv))
      {
        # check if analysis was done on aa or nt level
        if (input$nuctleotideFlag) obj = ntData else obj = mergedData
        sampForAnalysis = setdiff(names(obj),
                                  c(analysisRes$params$excludeSamp,analysisRes$params$refSamp))
        posClones = getPositiveClones(analysisRes$res, obj,
                                      samp = sampForAnalysis,
                                      orThr = as.numeric(analysisRes$params$orThr),
                                      fdrThr=as.numeric(analysisRes$params$fdrThr),
                                      percentThr = as.numeric(analysisRes$params$percentThr))
        # if there is no positive clones,
        # with a message about that
        if (length(posClones)==0)
        {
          # message
          m = 'There are no positive clones. Try to adjust thresholds'
          # show in the app
          output$save_results = renderText(m)
          # save in pdf
          pdf(file)
          plot.new()
          text(0.5,0.5,m)
          dev.off()
        }
        if(length(posClones)==1)
        {
          m = 'There is only one positive clone. Try to adjust thresholds to get more clones to plot'
          # show in the app
          output$save_results = renderText(m)
          # save in pdf
          pdf(file)
          plot.new()
          text(0.5,0.5,m)
          dev.off()
        }
        if (length(posClones)>1)
        {
          # create a table with results
          resTable = createResTable(analysisRes$res,obj,
                                    orThr = as.numeric(analysisRes$params$orThr),
                                    fdrThr = as.numeric(analysisRes$params$fdrThr),
                                    percentThr = as.numeric(analysisRes$params$percentThr),
                                    saveCI = F)
          # make heatmap with all significant clones
          posClones = list(rownames(resTable))
          names(posClones) = 'Positive_clones'

          makeHeatmaps(posClones, obj,
                       samp = sampForAnalysis,
                       refSamp = analysisRes$params$refSamp,
                       fileName = file, size = 7)
          output$save_results = renderText('The heat map is saved')
        }

      }
    })
}

shinyApp(ui = ui, server = server)
