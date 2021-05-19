##########################
# LabelRows Test
##########################

context("labelRows + full-workflow")

test_that("row annotation tests", {
    data("plasma20")
    data("plasma30")

    p20 = metabData(plasma20, zero = TRUE, samples = "CHEAR")
    p30 = metabData(plasma30, zero = TRUE, samples = "CHEAR")

    p = metabCombiner(p30,p20, binGap = 0.0075)
    p = selectAnchors(p, windx = 0.03, windy = 0.02)

    set.seed(100)
    p = fit_gam(p, k = seq(14,20,2), family = "gaussian", iterFilter = 2,
                method = "GCV.Cp")

    scores = evaluateParams(p, A = seq(60,100,10), B = seq(10,15),
                            C = seq(0,0.5,0.1), minScore = 0.7)
    p = calcScores(p, A = scores$A[1], B = scores$B[1], C = scores$C[1])

    p = labelRows(p, minScore = 0.5, maxRankX = 2, maxRankY = 2,
                    method = "score", delta = 0.1, remove = TRUE)

    p.output = combinedTable(p)

    labels = c("", "CONFLICT", "IDENTITY")
    testthat::expect_equal(sort(unique(p.output[["labels"]])), labels)
    testthat::expect_equal(sum(p.output[["labels"]] == "IDENTITY"), 527)

    p.output.2 = dplyr::filter(p.output, .data$labels != "IDENTITY")
    testthat::expect_equal(max(p.output.2[["rankX"]]), 2)
    testthat::expect_equal(max(p.output.2[["rankY"]]), 2)
    testthat::expect_false(any(p.output.2[["score"]] < 0.5))

    testthat::expect_error(labelRows(p, method = "mzrt", delta = 0.1))
    testthat::expect_error(labelRows(p, minScore = 1, delta = 0.1))
    testthat::expect_error(labelRows(p, maxRankX = 0, delta = 0.1))
})



