library(Matrix)
set.seed(123)

# peaks with gamma-distributed accessibility and sample-specific variation
simulateCounts <- function(meanCount, nPeaks, nSamples){
  acc <- rgamma(nPeaks, shape=2, rate=2) * rgamma(nPeaks*nSamples, 5, 5)
  counts <- matrix(rpois(nPeaks*nSamples, meanCount*acc), nrow=nPeaks)
  counts[rowSums(counts)>0,]
}

# ratio of the analytic background SD of random peak sets to the SD of their
# counts over sampled background sets, for each set and sample
backgroundSDratio <- function(counts, nSets=10, setSize=50, niterations=500){
  bg <- getBackgroundBins(counts, bias=runif(nrow(counts)), bs=20,
                          verbose=FALSE)
  bg <- computeBackgrounds(counts, bg)
  bgPeaks <- sampleBackgroundPeaks(bg, niterations=niterations)
  V <- as.matrix(bg@V)[bg@peak2bin,]
  counts <- as.matrix(counts)
  unlist(lapply(seq_len(nSets), function(i){
    peaks <- sample(nrow(counts), setSize)
    empiricalSD <- apply(counts, 2, function(x){
      sd(colSums(matrix(x[bgPeaks[peaks,]], nrow=setSize)))
    })
    sqrt(colSums(V[peaks,]))/empiricalSD
  }))
}

test_that("analytic background SD matches sampled background sets", {
  # median ratios should be within 3% of 1 (within 1% across seeds; 1.1.10
  # gave about 2.6 for pseudobulk-like and 1.07 for single-cell-like counts)
  # pseudobulk-like depth
  r <- backgroundSDratio(simulateCounts(200, 2000, 10))
  expect_lt(abs(median(r)-1), 0.03)
  # single-cell-like depth
  r <- backgroundSDratio(as(simulateCounts(0.05, 2000, 200), "CsparseMatrix"))
  expect_lt(abs(median(r)-1), 0.03)
})

test_that("the bin-to-bin probabilities are balanced", {
  counts <- simulateCounts(200, 4000, 6)
  bias <- runif(nrow(counts))
  getBins <- function(...) getBackgroundBins(counts, bias=bias, bs=20, 
                                             verbose=FALSE, ...)
  # row sums (the weights of the background of each bin) and column sums (the 
  # weights with which each bin is used in the background of others)
  sums <- function(bg){
    n <- bg@binDensity
    list(rows=as.numeric(bg@binBinProbs %*% n)[n>0],
         cols=as.numeric(Matrix::crossprod(bg@binBinProbs, n))[n>0])
  }
  s <- sums(getBins())
  expect_lt(max(abs(s$rows-1)), 1e-6)
  expect_lt(max(abs(s$cols-1)), 1e-3)
  expect_lt(max(abs(sums(getBins(tol=1e-4))$cols-1)), 1e-4)
  # without balancing, the adaptive kernel leaves them unbalanced
  s0 <- sums(getBins(tol=NULL))
  expect_lt(min(s0$cols), 0.5)
  # also without adaptive kernel
  expect_lt(min(sums(getBins(minCount=NULL, tol=NULL))$cols), 0.8)
  s <- sums(getBins(minCount=NULL))
  expect_lt(max(abs(s$rows-1)), 1e-4)
  expect_lt(max(abs(s$cols-1)), 1e-3)
  expect_error(getBins(tol=-1))
  expect_error(getBins(tol=c(1e-3, 1e-4)))
  expect_warning(getBins(tol=1e-12), "could not be balanced")
  # the sampling of background peaks still works
  expect_equal(dim(sampleBackgroundPeaks(getBins(), niterations=3)), 
               c(nrow(counts), 3))
})

test_that("the balancing and the ratio correction remove the bias of random sets", {
  counts <- simulateCounts(200, 4000, 6)
  bias <- runif(nrow(counts))
  nSets <- 60
  sets <- sparseMatrix(i=unlist(lapply(seq_len(nSets), function(i) 
                                       sample.int(nrow(counts), 1500))),
                       j=rep(seq_len(nSets), each=1500), x=1,
                       dims=c(nrow(counts), nSets))
  bgb <- function(...) getBackgroundBins(counts, bias=bias, bs=20, 
                                         verbose=FALSE, ...)
  # without balancing, peaks are not used as background as often as they are 
  # targets; with it, they are
  colUse <- function(bg) 
    as.numeric(Matrix::crossprod(bg@binBinProbs, bg@binDensity))[bg@binDensity>0]
  expect_gt(max(abs(colUse(bgb(tol=NULL))-1)), 0.1)
  expect_lt(max(abs(colUse(bgb())-1)), 0.01)
  dev <- function(...){
    bg <- computeBackgrounds(counts, bgb(...))
    assay(computeDeviationsAnalytic(counts, bg, sets), "deviations")
  }
  # the ratio correction makes the deviations of random sets unbiased even 
  # when the bins are not balanced
  expect_lt(abs(mean(dev(tol=NULL))), 0.002)
  expect_lt(abs(mean(dev())), 0.002)
})


test_that("minCount widens the kernel of sparse bins only, and maxShare is deprecated", {
  counts <- simulateCounts(200, 4000, 6)
  bias <- runif(nrow(counts))
  getBins <- function(...) getBackgroundBins(counts, bias=bias, bs=20, tol=NULL,
                                             verbose=FALSE, ...)
  b0 <- getBins(minCount=NULL)
  b1 <- getBins(minCount=10)
  n <- b0@binDensity
  changed <- Matrix::rowSums(abs(b0@binBinProbs-b1@binBinProbs)) > 1e-12
  expect_true(any(changed))
  # only the occupied bins with fewer than minCount peaks are modified
  expect_true(all(n[changed] > 0 & n[changed] < 10))
  expect_identical(b0@peak2bin, b1@peak2bin)
  # widening increases the effective background size of the changed bins
  ess <- function(b) 1/as.numeric((b@binBinProbs*b@binBinProbs) %*% b@binDensity)
  expect_true(all(ess(b1)[changed] >= ess(b0)[changed]))
  # minCount=0 is the same as NULL
  expect_equal(getBins(minCount=0)@binBinProbs, b0@binBinProbs)
  expect_warning(bm <- getBins(maxShare=0.1), "deprecated")
  expect_equal(bm@binBinProbs, b1@binBinProbs)
})
