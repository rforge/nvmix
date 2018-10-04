context("Special Cases")

## Some inputs that can be used:
d <- 3
n <- 10
set.seed(157)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A %*% t(A))
b <-  3 * runif(d) * sqrt(d) # random upper limit



test_that("pnvmix() deals correctly with 'upper = rep(Inf,d)', 'lower = rep(-Inf,d)'.", {
  expect_equal(pnvmix(upper = rep(Inf, d), qmix = "constant"), expect = 1, check.attributes = FALSE)
})


test_that("pnvmix() throws error when 'scale' not positive definite.", {
  ## Define singular 'scale' matrix 
  scale <- matrix(c(1,0,0,
                    0,1,0,
                    0,0,0), ncol = 3, byrow = TRUE)
  expect_error(pnvmix(upper = c(1,2,3), qmix = "constant", scale = scale))
})


test_that("pStudent deals with 'df = Inf'.", {
  ## Reset seed to get exactly the same results
  set.seed(1)
  p.St.dfInf <- pStudent(b, df = Inf, scale = P)
  set.seed(1)
  p.Norm <- pNorm(b, scale = P)
  expect_equal(p.St.dfInf, p.Norm, check.attributes = FALSE)
})


test_that("rnvmix can have non-square 'factor'.", {
  A <- matrix( c(1,0,
                 0,1,
                 0,1), ncol = 3)
  
  expect_equal( dim(rnvmix(n, rmix = "constant", factor = A)), c(n, 3))
  expect_equal( dim(rnvmix(n, rmix = "constant", factor = t(A))), c(n, 2))
})