library(Matrix)

# simulates (peaks x samples) counts with a controlled sample-specific bias
simBias <- function(N=20000, S=8, gcEffect=0, accEffect=0, seed=1){
  set.seed(seed)
  gc <- pmin(0.8, pmax(0.25, rnorm(N, 0.45, 0.08)))
  lam <- 10^rnorm(N, 1.5, 0.5)
  z <- scale(log10(lam))[,1]
  gz <- scale(gc)[,1]
  beta <- matrix(rnorm(2*S), 2, S)   # sample-specific slopes
  lf <- outer(gz, beta[1,]*gcEffect) + outer(z, beta[2,]*accEffect)
  mu <- lam * exp(lf)
  m <- matrix(rnbinom(N*S, mu=mu, size=20), N, S)
  m[rowSums(m)==0, 1] <- 1
  list(counts=m, gc=gc)
}

test_that("testBiasDimensions detects the absence of bias", {
  d <- simBias()
  res <- testBiasDimensions(d$counts, bias=d$gc, folds=4)
  expect_s3_class(res, "biasDimTest")
  expect_true(all(res$models$gain[res$models$stratum=="all"] < 0.002))
  expect_output(print(res), "Cross-validated")
})

test_that("testBiasDimensions identifies GC-only and accessibility bias", {
  d <- simBias(gcEffect=0.15)
  res <- testBiasDimensions(d$counts, bias=d$gc, folds=4)
  imp <- res$importance[res$importance$stratum=="all",]
  imp <- setNames(imp$share, imp$dimension)
  expect_gt(imp[["GC"]], 0.8)
  expect_lt(imp[["accessibility"]], 0.2)
  d <- simBias(accEffect=0.15)
  res <- testBiasDimensions(d$counts, bias=d$gc, folds=4)
  imp <- res$importance[res$importance$stratum=="all",]
  imp <- setNames(imp$share, imp$dimension)
  expect_gt(imp[["accessibility"]], 0.8)
  expect_lt(imp[["GC"]], 0.2)
})

test_that("testBiasDimensions supports flbias, grouping and SE input", {
  d <- simBias(N=6000, gcEffect=0.1)
  fl <- rnorm(nrow(d$counts))
  res <- testBiasDimensions(d$counts, bias=d$gc, flbias=fl, folds=3,
                            grouping=rep(1:2, each=4))
  expect_true(all(c("accessibility","GC","flbias") %in% res$importance$dimension))
  expect_equal(nrow(res$models), 4*7)
  se <- SummarizedExperiment::SummarizedExperiment(list(counts=d$counts),
                                                  rowData=S4Vectors::DataFrame(bias=d$gc))
  res2 <- testBiasDimensions(se, folds=3, useWidthAdj=FALSE)
  expect_s3_class(res2, "biasDimTest")
})
