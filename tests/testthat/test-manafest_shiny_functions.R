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
                            excludeSamp = "",
                            compareToRef = TRUE,
                            saveToFile = F)
  # load previously saved results
  testRes = readRDS(file.path(dir,"no_replicate_results.rds"))
  # compare two results
  expect_equal(res, testRes)
  #expect_equal(readMergeSave(''), NULL)
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
                      nReads = 10,
                      percentThr = 0,
                      xrCond = NULL,
                      saveToFile = F)

   # load previously saved results
  testRes = readRDS(file.path(dir,"replicate_results.rds"))
  # compare two results
  expect_equal(res, testRes)
  #expect_equal(readMergeSave(''), NULL)
})

