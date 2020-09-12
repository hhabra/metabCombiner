###############################################
# Tests for the metabData object and methods
###############################################

context("metabData tests")

test_that("field detection tests",{
    data("plasma20")

    p20 = metabData(plasma20, samples = "Red", extra = "CHEAR")

    expect_length(getSamples(p20), 5)
    expect_length(getExtra(p20), 5)
    expect_false(any(getSamples(p20) %in% getExtra(p20)))

    data = getData(p20)
    expect_equal(nrow(plasma20), nrow(data))

    p20_2 = metabData(plasma20, samples = "RedCross", extra = "CH")
    expect_equivalent(p20, p20_2)

    plasma20$rt[25] = plasma20$rt[25] * -1
    expect_error(metabData(plasma20, samples = "RedCross", extra = "CHEAR"))

    plasma20$rt[25] = plasma20$rt[25] * -10000
    expect_warning(metabData(plasma20, samples = "RedCross", extra = "CHEAR"))

    plasma20$rt[25] = plasma20$rt[25] * 1/10000
    plasma20$mz[3] = plasma20$mz[3] * 100
    expect_warning(metabData(plasma20, samples = "RedCross", extra = "CHEAR"))

    plasma20$mz[3] = plasma20$mz[3] / 100
    p20 = suppressWarnings(metabData(plasma20))
    expect_equivalent(grep("POOL|CHEAR|Red|Bl", names(plasma20), value = TRUE),
                      getSamples(p20))

    expect_error(metabData(plasma20, mz = "mz", rt = "mz"))
    expect_error(metabData(plasma20, rt = "id"))
})

test_that("RT filtering test", {
    data("plasma20")
    p20 = metabData(plasma20, samples = "CHEAR", misspc = 15, rtmax = 17.25)

    expect_lt(max(getData(p20)[["rt"]]), 17.25)
    expect_gt(getStats(p20)[["filtered_by_rt"]], 0)
})


test_that("filter_missingness test", {
    data("plasma20")

    plasma20[["CHEAR.20min.2"]] = NA
    plasma20[["CHEAR.20min.2"]][100] = 1000
    plasma20[["CHEAR.20min.3"]][100] = NA

    p20 = metabData(plasma20, samples = "CHEAR", misspc = 15)

    expect_equal(getStats(p20)[["final_count"]], 0)

})

test_that("filter_duplicates test", {
    mz = seq(100, 110, 0.01)
    rt = rep(10,1001)
    s = seq(1001.1,1.1, -1)
    data = data.frame(mz = mz, rt = rt, s = s)
    d = suppressWarnings(metabData(data, duplicate = c(0.01,0.01)))
    expect_equal(getStats(d)[["filtered_as_duplicates"]], 500)
})





