set.seed(42)
N <- 3000; S <- 8
gc <- runif(N, 0.3, 0.7)
lam <- 10^rnorm(N, 1.2, 0.5)
depth <- runif(S, 0.6, 1.6)
m <- matrix(rpois(N*S, outer(lam, depth) * exp(outer((gc-0.5)*6, rnorm(S, 0, 0.3)))), N, S)
grp <- rep(c("A","B"), each=4)
# continuous version (no ties), for properties that hold exactly only without ties
mc <- matrix(rgamma(N*S, shape=2, rate=0.1/outer(lam, depth)), N, S)

test_that("GCQN equalizes the distributions within each GC bin", {
  o <- GCQN(mc, bias=gc, g=5, round=FALSE)
  expect_equal(dim(o), dim(mc))
  cuts <- unique(quantile(gc, (1:5)/5, type=1, names=FALSE))
  bins <- pmin(length(cuts), findInterval(gc, cuts, left.open=TRUE)+1L)
  for(b in unique(bins)){
    s <- apply(o[bins==b,], 2, sort)
    expect_equal(s, matrix(s[,1], nrow(s), ncol(s)), ignore_attr=TRUE)
  }
  expect_true(all(o>=0))
  expect_true(all(GCQN(m, bias=gc)==round(GCQN(m, bias=gc))))
})

test_that("GCQN with median and mean references", {
  ref <- function(m, f){
    q <- apply(m, 2, sort); r <- apply(q, 1, f)
    apply(m, 2, function(x) r[rank(x, ties.method="min")])
  }
  o <- GCQN(m, bias=gc, g=1, summary="median", round=FALSE)
  expect_equal(o, ref(m, median), ignore_attr=TRUE)
  o <- GCQN(m, bias=gc, g=1, summary="mean", round=FALSE)
  expect_equal(o, ref(m, mean), ignore_attr=TRUE)
  expect_error(GCQN(m, bias=gc, summary="mode"))
})

test_that("smooth GCQN (mean) reproduces qsmooth", {
  skip_if_not_installed("qsmooth")
  o <- GCQN(m, bias=gc, g=1, grouping=grp, summary="mean", round=FALSE)
  q <- qsmooth::qsmooth(m, group_factor=factor(grp))@qsmoothData
  expect_equal(unname(o), unname(q), tolerance=1e-8)
})

test_that("grouping preserves part of the differences between groups", {
  m2 <- mc; m2[,5:8] <- m2[,5:8]*1.5   # global group effect
  plain <- GCQN(m2, bias=gc, g=1, round=FALSE)
  smooth <- GCQN(m2, bias=gc, g=1, grouping=grp, round=FALSE)
  d <- function(x) mean(x[,5:8]) - mean(x[,1:4])
  expect_lt(abs(d(plain)), 1e-6)
  expect_gt(d(smooth), 0.1*d(m2))
  # a single group is the same as no grouping
  expect_equal(GCQN(m, bias=gc, grouping=rep("A", S), summary="median"), 
               GCQN(m, bias=gc, summary="median"))
})

test_that("GCQN supports SummarizedExperiment objects and fails properly", {
  se <- SummarizedExperiment::SummarizedExperiment(list(counts=m), 
    rowData=S4Vectors::DataFrame(bias=gc), 
    colData=S4Vectors::DataFrame(grp=grp))
  o <- GCQN(se, grouping="grp")
  expect_true("GCQN" %in% SummarizedExperiment::assayNames(o))
  expect_equal(SummarizedExperiment::assay(o, "GCQN"), 
               GCQN(m, bias=gc, grouping=grp), ignore_attr=TRUE)
  expect_error(GCQN(m), "bias")
  mm <- m; mm[1,1] <- NA
  expect_error(GCQN(mm, bias=gc), "missing")
})
