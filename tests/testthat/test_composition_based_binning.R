test_that("test parameters",{
    expect_error(composition_based_binning())
    expect_error(composition_based_binning(1))
    expect_error(composition_based_binning("thisdoesnotexist"))
    expect_error(composition_based_binning(fastaFile, kMerSizeCB = 0))
    expect_error(composition_based_binning(fastaFile, kMerSizeCB = "some"))
    expect_error(composition_based_binning(fastaFile, numOfClustersCB = 1))
    expect_error(composition_based_binning(fastaFile, numOfClustersCB = "some"))
    expect_error(composition_based_binning(fastaFile, outputCB = 1))
    expect_error(composition_based_binning(fastaFile, outputCB = ""))
    expect_error(composition_based_binning(fastaFile, keepQuality = 10))
    expect_error(composition_based_binning(fastaFile, keepQuality = "some"))
    expect_error(composition_based_binning(fastaFile, numOfThreads = 0))
    expect_error(composition_based_binning(fastaFile, numOfThreads = "some"))
})

test_that("test return is a list",{
    expect_type(
        composition_based_binning(
            system.file("extdata", "reads.metagenome.fasta.gz",
                        package = "metabinR"),
            dryRun = TRUE, kMerSizeCB = 4, numOfClustersCB = 2),
        "list")
})
