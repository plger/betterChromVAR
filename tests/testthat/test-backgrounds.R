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
  expect_lt(min(sums(getBins(maxShare=NULL, tol=NULL))$cols), 0.8)
  s <- sums(getBins(maxShare=NULL))
  expect_lt(max(abs(s$rows-1)), 1e-4)
  expect_lt(max(abs(s$cols-1)), 1e-3)
  expect_error(getBins(tol=-1))
  expect_error(getBins(tol=c(1e-3, 1e-4)))
  expect_warning(getBins(tol=1e-12), "could not be balanced")
  # the sampling of background peaks still works
  expect_equal(dim(sampleBackgroundPeaks(getBins(), niterations=3)), 
               c(nrow(counts), 3))
})

test_that("the balancing removes the bias of random sets of peaks", {
  counts <- simulateCounts(200, 4000, 6)
  bias <- runif(nrow(counts))
  nSets <- 60
  sets <- sparseMatrix(i=unlist(lapply(seq_len(nSets), function(i) 
                                       sample.int(nrow(counts), 1500))),
                       j=rep(seq_len(nSets), each=1500), x=1,
                       dims=c(nrow(counts), nSets))
  dev <- function(...){
    bg <- getBackgroundBins(counts, bias=bias, bs=20, verbose=FALSE, ...)
    bg <- computeBackgrounds(counts, bg)
    assay(computeDeviationsAnalytic(counts, bg, sets), "deviations")
  }
  # positive bias without balancing (about 0.8%), absent with it
  expect_gt(mean(dev(tol=NULL)), 0.004)
  expect_lt(abs(mean(dev())), 0.002)
})
