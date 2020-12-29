##########################
# fit_gam Test
##########################

test_that("fit_gam outlier detection", {
      #test set-up: m/z, rt, id and values
      mz1 <- seq(100, 150, 0.005)
      n <- length(mz1)
      set.seed(100)
      mz2 <- seq(100, 150, 0.005) + rnorm(n, sd = 0.015)
      rt1 <- seq(0.5, 10.5, 0.001)
      rt2 <- c(seq(0.4,3.05, 0.001), seq(3.05, 12, 0.0033333),
               seq(12,16,0.002), seq(16, 20, 0.0015025))
      id1 <- rep("", n)
      id1[seq(1, n, 100)] <- paste("C", seq(1, n, 100), sep = "")
      id2 <- id1
      set.seed(100); val1 <- runif(n, min = 1000, max = 1000000)
      set.seed(101); val2 <- runif(n, min = 1000, max = 100000)
      data1 <- data.frame(mz = mz1, rt = rt1, id = id1, val = val1)
      data2 <- data.frame(mz = mz2, rt = rt2, id = id2, val = val2)

      #metabCombiner objects
      md1 <- metabData(data1, samples = "val")
      md2 <- metabData(data2, samples = "val")
      mc <- metabCombiner(md1, md2, binGap = 0.005)
      mc <- selectAnchors(mc, useID = TRUE)

      #adding outliers to set of ordered pairs
      anchors <- getAnchors(mc)
      anchors$rty[seq(0,100,25)] <- anchors$rty[seq(0,100,25)] + 2.5
      anchors$rty[seq(0,100,10)] <- anchors$rty[seq(0,100,10)] - 1.5
      anchors$rty[seq(0,100,5)] <- anchors$rty[seq(0,100,5)] + 0.5
      mc@anchors <- anchors

      ##fit_gam tests
      mc <- fit_gam(mc, outlier = "MAD", coef = 2, family = "gaussian",
                    message = FALSE, iterFilter = 2, k = seq(6,12,2))
      weights <- getModel(mc)[["weights"]]
      testthat::expect_equal(sum(seq(6,101,5) %in% which(weights == 0)), 20)

      mc <- fit_gam(mc, outlier = "boxplot", coef = 1.5, family = "gaussian",
                    message = FALSE, iterFilter = 3, bs = "ps", k = seq(7,12))
      weights <- getModel(mc)[["weights"]]
      testthat::expect_equal(sum(seq(6,101,5) %in% which(weights == 0)), 20)
})







