test_that("test parameters",{
    expect_error(abundance_based_binning())
    expect_error(abundance_based_binning(1))
    expect_error(abundance_based_binning("thisdoesnotexist"))
    expect_error(abundance_based_binning(fastaFile, eMin = 0))
    expect_error(abundance_based_binning(fastaFile, eMin = ""))
    expect_error(abundance_based_binning(fastaFile, eMax = -1))
    expect_error(abundance_based_binning(fastaFile, eMax = ""))
    expect_error(abundance_based_binning(fastaFile, kMerSizeAB = 0))
    expect_error(abundance_based_binning(fastaFile, kMerSizeAB = "some"))
    expect_error(abundance_based_binning(fastaFile, numOfClustersAB = 1))
    expect_error(abundance_based_binning(fastaFile, numOfClustersAB = "some"))
    expect_error(abundance_based_binning(fastaFile, outputAB = 1))
    expect_error(abundance_based_binning(fastaFile, outputAB = ""))
    expect_error(abundance_based_binning(fastaFile, keepQuality = 10))
    expect_error(abundance_based_binning(fastaFile, keepQuality = "some"))
    expect_error(abundance_based_binning(fastaFile, numOfThreads = 0))
    expect_error(abundance_based_binning(fastaFile, numOfThreads = "some"))
})

test_that("test return is a list",{
    expect_type(
        abundance_based_binning(
            system.file("extdata", "reads.metagenome.fasta.gz",
                        package = "metabinR"),
            dryRun = TRUE, kMerSizeAB = 8, numOfClustersAB = 2),
    "list")
})
