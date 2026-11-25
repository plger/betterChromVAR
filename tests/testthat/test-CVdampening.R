simDamp <- function(effect=c("none","wide"), fc=2, N=15000, S=8, seed=1){
  effect <- match.arg(effect); set.seed(seed)
  gc <- pmin(0.8, pmax(0.25, rnorm(N, 0.45, 0.08)))
  lam <- 10^rnorm(N, 1.5, 0.5)
  depth <- runif(S, 0.7, 1.4)
  mu <- outer(lam, depth)
  grp <- rep(c("A","B"), each=S/2)
  if(effect=="wide"){   # group effect spanning the top of the accessibility range
    idx <- which(lam > quantile(lam, 0.7))
    mu[idx, grp=="B"] <- mu[idx, grp=="B"]*fc
  }
  m <- matrix(rnbinom(N*S, mu=mu, size=20), N, S)
  m[rowSums(m)==0, 1] <- 1
  list(m=m, gc=gc, grp=grp)
}

test_that("no warning when the groups do not drive the bin-level bias", {
  d <- simDamp("none")
  expect_no_warning(a <- CVnorm(d$m, bias=d$gc, useWidthAdj=FALSE, 
                                grouping=d$grp))
  # the warning does not change the output
  b <- CVnorm(d$m, bias=d$gc, useWidthAdj=FALSE, grouping=d$grp, 
              warnWeight=NULL)
  expect_equal(a, b)
})

test_that("a group effect spanning the accessibility range triggers a warning", {
  d <- simDamp("wide")
  expect_warning(
    a <- CVnorm(d$m, bias=d$gc, useWidthAdj=FALSE, grouping=d$grp),
    "mean dampening weight")
  expect_no_warning(
    b <- CVnorm(d$m, bias=d$gc, useWidthAdj=FALSE, grouping=d$grp, 
                warnWeight=NULL))
  expect_equal(a, b)
  # no grouping, no dampening, no warning
  expect_no_warning(CVnorm(d$m, bias=d$gc, useWidthAdj=FALSE))
})
