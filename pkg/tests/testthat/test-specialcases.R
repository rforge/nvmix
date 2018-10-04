context("Special Cases")

test_that("pNorm() deals correctly with d = 1, upper = Inf, lower = -Inf", {
  expect_equal(pnvmix(upper = Inf, lower = -Inf, qmix = "constant"), expect = 1, check.attributes = FALSE)
}
)

