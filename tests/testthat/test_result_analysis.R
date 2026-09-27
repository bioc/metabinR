make_analysis_result <- function(algorithm, assignments) {
    methods::new(
        "MetabinResult",
        assignments = S4Vectors::DataFrame(assignments, check.names = FALSE),
        parameters = list(), inputs = character(), algorithm = algorithm
    )
}

test_that("bin_summary counts observed bins and assigned distances", {
    x <- make_analysis_result("CB", data.frame(
        read_id = c("a", "b", "c"), CB = c(2L, 1L, 2L),
        CB.1 = c(0.8, 0.1, 0.9), CB.2 = c(0.2, 0.4, 0.3)
    ))
    result <- as.data.frame(bin_summary(x))
    expect_equal(result$bin, c("2", "1"))
    expect_equal(result$n_reads, c(2L, 1L))
    expect_equal(result$proportion, c(2 / 3, 1 / 3))
    expect_equal(result$mean_distance, c(0.25, 0.1))
    expect_equal(result$median_distance, c(0.25, 0.1))
})

test_that("ambiguous_reads handles ties and the margin boundary", {
    x <- make_analysis_result("AB", data.frame(
        read_id = c("tie", "edge", "clear"), AB = c(1L, 2L, 1L),
        AB.1 = c(0.2, 0.5, 0.1), AB.2 = c(0.2, 0.25, 0.8)
    ))
    result <- as.data.frame(ambiguous_reads(x, margin = 0.25))
    expect_equal(result$read_id, c("tie", "edge"))
    expect_equal(result$bin, c("1", "2"))
    expect_equal(result$distance_margin, c(0, 0.25))
    expect_error(ambiguous_reads(x, margin = -1))
})

test_that("hierarchical ambiguity uses only finite candidate distances", {
    x <- make_analysis_result("ABxCB", data.frame(
        read_id = c("a", "b", "c"), ABxCB = c(1L, 3L, 1L),
        ABxCB.1 = c(0.2, NA_real_, 0.1),
        ABxCB.2 = c(0.22, NA_real_, NA_real_),
        ABxCB.3 = c(NA_real_, 0.3, NA_real_),
        ABxCB.4 = c(NA_real_, 0.8, NA_real_)
    ))
    expect_equal(as.data.frame(ambiguous_reads(x, margin = 0.05))$read_id, "a")
    expect_equal(as.data.frame(bin_summary(x))$n_reads, c(2L, 1L))
})
