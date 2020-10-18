##########################
# Tests for selectAnchors
##########################

context("selectAnchors tests")

test_that("anchor selection", {
    data("plasma20")
    data("plasma30")

    p20 = metabData(plasma20, zero = TRUE, samples = "CHEAR")
    p30 = metabData(plasma30, zero = TRUE, samples = "CHEAR")

    p = metabCombiner(p30,p20, binGap = 0.0075)
    p = selectAnchors(p, windx = 0.03, windy = 0.02)
    anchors = getAnchors(p)
    testthat::expect_false(any(abs(anchors[["mzx"]]-anchors[["mzy"]]) > 0.003))
    testthat::expect_false(any(abs(anchors[["Qx"]]-anchors[["Qy"]]) > 0.3))

    rts = sort(anchors[["rtx"]])
    rtdiffs = rts[seq(2,length(rts))] - rts[seq(1,length(rts)-1)]
    testthat::expect_true(min(rtdiffs) > 0.03)

    rts = sort(anchors[["rty"]])
    rtdiffs = rts[seq(2,length(rts))] - rts[seq(1,length(rts)-1)]
    testthat::expect_true(min(rtdiffs) > 0.02)

    p = selectAnchors(p, useID = TRUE)
    anchors = getAnchors(p)
    testthat::expect_equal(sum(anchors[["labels"]] == "I"), 136)
})



