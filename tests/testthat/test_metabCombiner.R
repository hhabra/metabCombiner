###############################################
# Tests for the metabCombiner object and methods
###############################################

context("metabCombiner tests")

test_that("metabCombiner creation", {
    mzx = c(seq(50, 50 + 0.01*500, 0.01), seq(50.002, 50.002 + 0.01*499, 0.01))
    rtx = 0.5 + seq(0, 10,0.01)
    sx = seq(1001.1,1.1, -1)
    datax = data.frame(mz = mzx, rt = rtx, sx = sx)
    dx = suppressWarnings(metabData(datax, duplicate = c(0,0)))

    mzy = sort(c(seq(60, 60 + 0.01*400, 0.01),
                 seq(50.002, 50.002 + 0.01*499, 0.01)))
    rty = 0.7 +  seq(0,9,0.01)
    sy = seq(901.1,1.1,-1)
    datay = data.frame(mz = mzy, rt = rty, sy = sy)
    dy = suppressWarnings(metabData(datay, duplicate = c(0,0)))
    dxdy = metabCombiner(dx, dy, binGap = 0.005)

    testthat::expect_equal(getSamples(dxdy), "sx")
    testthat::expect_equal(getSamples(dxdy, data = "y"), "sy")
    testthat::expect_equal(getStats(dxdy)[["nGroups"]], 500)
})


