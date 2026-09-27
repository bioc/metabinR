make_evaluation_result <- function(bins, ids = paste0("r", seq_along(bins)),
                                   algorithm = "CB") {
    methods::new(
        "MetabinResult",
        assignments = S4Vectors::DataFrame(
            read_id = ids, CB = bins
        ),
        parameters = list(), inputs = character(), algorithm = algorithm
    )
}

test_that("evaluation matches read IDs and ignores cluster label names", {
    x <- make_evaluation_result(c(2L, 2L, 1L, 1L))
    truth <- data.frame(
        anonymous_read_id = c("r4", "extra", "r2", "r1", "r3"),
        genome_id = c("A", "unused", "B", "B", "A")
    )
    result <- evaluate_bins(x, truth, id = "anonymous_read_id")
    expect_equal(unname(result$confusion), matrix(c(2L, 0L, 0L, 2L), 2L))
    expect_equal(as.data.frame(result$per_bin)$purity, c(1, 1))
    expect_equal(as.data.frame(result$per_origin)$best_bin_recovery, c(1, 1))
    expect_equal(result$overall$adjusted_rand_index, 1)
})

test_that("mixed bins have the expected read-level scores", {
    x <- make_evaluation_result(c(1L, 2L, 1L, 2L))
    truth <- data.frame(read_id = paste0("r", 1:4),
                        genome_id = c("A", "A", "B", "B"))
    result <- evaluate_bins(x, truth)
    expect_equal(result$overall$weighted_purity, 0.5)
    expect_equal(result$overall$weighted_recovery, 0.5)
    expect_equal(result$overall$adjusted_rand_index, -0.5)
})

test_that("evaluation rejects ambiguous or incomplete truth mappings", {
    x <- make_evaluation_result(c(1L, 2L))
    truth <- data.frame(read_id = "r1", genome_id = "A")
    expect_error(evaluate_bins(x, truth), "missing 1 result read")
    truth <- data.frame(read_id = c("r1", "r1", "r2"),
                        genome_id = c("A", "A", "B"))
    expect_error(evaluate_bins(x, truth), "unique")
    expect_error(evaluate_bins(make_evaluation_result(c(1L, 2L),
                                                      ids = c("r1", "r1")),
                               truth), "unique")
    expect_error(evaluate_bins(x, truth, label = "missing"), "lacks")
})

test_that("ARI handles one read and trivial identical partitions", {
    one <- make_evaluation_result(1L)
    expect_true(is.na(evaluate_bins(
        one, data.frame(read_id = "r1", genome_id = "A")
    )$overall$adjusted_rand_index))
    three <- make_evaluation_result(c(1L, 1L, 1L))
    truth <- data.frame(read_id = paste0("r", 1:3),
                        genome_id = c("A", "A", "A"))
    expect_equal(evaluate_bins(three, truth)$overall$adjusted_rand_index, 1)
    expect_equal(metabinR:::.adjusted_rand_index(
        matrix(c(50000L, 0L, 0L, 50000L), nrow = 2L)
    ), 1)
})
