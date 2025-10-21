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
