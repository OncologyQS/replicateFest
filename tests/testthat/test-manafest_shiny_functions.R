test_that("reading input file", {
  expect_equal(readMergeSave(NULL), NULL)
  #expect_equal(readMergeSave(''), NULL)
})

test_that("running analysis without replicates", {
  # find folder with input data
  dir = testthat::test_path("testdata", "no_replicates")
  # list paths to files with data
  files = list.files(dir, full.names = T,
                     pattern = "tsv", recursive = TRUE)
  res = runExperimentFisher(files,
                            refSamp = "sample1_control",
                            nReads = 10,
                            fdrThr = .05,
                            orThr = 1,
                            percentThr = 0,
                            condsThr = 0,
                            excludeSamp = "",
                            compareToRef = TRUE,
                            saveToFile = F)
  # load previously saved results
  testRes = readRDS(file.path(dir,"no_replicate_results.rds"))
  # compare two results
  expect_equal(res, testRes)
  #expect_equal(readMergeSave(''), NULL)
})

test_that("running analysis without replicates and no comparison to reference", {
  # find folder with input data
  dir = testthat::test_path("testdata", "no_replicates")
  # list paths to files with data
  files = list.files(dir, full.names = T,
                     pattern = "tsv", recursive = TRUE)
  res = runExperimentFisher(files,
                            refSamp = "",
                            nReads = 10,
                            fdrThr = .05,
                            orThr = 1,
                            percentThr = 0,
                            condsThr = 0,
                            excludeSamp = "sample1_control",
                            compareToRef = FALSE,
                            saveToFile = F)
  # load previously saved results
  testRes = readRDS(file.path(dir,"no_replicate_no_ref_results.rds"))
  # compare two results
  expect_equal(res, testRes)
})

test_that("running analysis with replicates", {
  # find folder with input data
  dir = testthat::test_path("testdata", "with_replicates")
  # list paths to files with data
  # list paths to files with data
  files = list.files(dir, full.names = T,
                     pattern = "tsv", recursive = TRUE)
  filenames = file_path_sans_ext(basename(files))
  sampAnnot = splitFileName(filenames)

  # run all clones in a patient and time point and return the results
  res = runExperiment(files,
                      peptides = sampAnnot$condition,
                      refSamp= "control",
                      fdrThr = 0.05,
                      nReads = 20,
                      percentThr = 0,
                      xrCond = NULL,
                      saveToFile = F)

   # load previously saved results
  testRes = readRDS(file.path(dir,"replicate_results.rds"))
  # compare two results
  expect_equal(res, testRes)
  #expect_equal(readMergeSave(''), NULL)
})


# make heatmap test
test_that("making a heatmap", {
  # find folder with input data
  dir = testthat::test_path("testdata", "no_replicates")
  # load input data
  load(file = file.path(dir,"no_replicate_inputData.rda"))
  # load previously saved results
  testRes = readRDS(file.path(dir,"no_replicate_no_ref_results.rds"))
  posClones = rownames(testRes$positive_clones_all_data)
  makeHeatmaps(posClones, mergedData = mergedData,
               refSamp = "sample1_control")
  expect_equal(makeHeatmaps(posClones, mergedData = mergedData,
                            refSamp = "sample1_control"), NULL)
  #expect_equal(readMergeSave(''), NULL)
})



