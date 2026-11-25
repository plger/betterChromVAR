library(BiocParallel)
library(Matrix)
library(SummarizedExperiment)
set.seed(123)
attach(getDummyData())

test_that("the finite-population correction does not affect the deviations", {
  bg <- computeBackgrounds(counts, getBackgroundBins(counts, verbose=FALSE))
  d1 <- computeDeviationsAnalytic(counts, bg, motifMatches)
  d0 <- computeDeviationsAnalytic(counts, bg, motifMatches, fpc=FALSE)
  expect_identical(assay(d1, "deviations"), assay(d0, "deviations"))
  # it reduces the variance, i.e. increases the z-scores
  expect_true(all(abs(assay(d1, "z")) >= abs(assay(d0, "z"))))
  expect_true(all(abs(assay(d1, "z")) > abs(assay(d0, "z"))))
})

test_that("the finite-population variance matches a peak-level computation", {
  bins <- getBackgroundBins(counts, verbose=FALSE)
  bg <- computeBackgrounds(counts, bins)
  st <- betterChromVAR:::.fpcStats(motifMatches, bg)
  # variance of sum_q (a_q - g_b(q)) x_q with independent peaks of variance V_b
  b <- bg@peak2bin
  direct <- sapply(seq_len(ncol(motifMatches)), function(m){
    cq <- motifMatches[,m] - st$g[b,m]
    colSums(cq^2 * as.matrix(bg@V)[b,,drop=FALSE])
  })
  expect_equal(t(direct), unname(betterChromVAR:::.fpcVariance(motifMatches, bg, bg@V)),
               ignore_attr=TRUE)
  # effective number of peaks, from the peak level
  binTotals <- numeric(length(bg@binDensity))
  binTotals[] <- rowsum(bg@expectation, b)[match(seq_along(binTotals), sort(unique(b)))]
  binTotals[is.na(binTotals)] <- 0
  Eb <- as.numeric(bg@binBinProbs %*% binTotals)
  neff <- sapply(seq_len(ncol(motifMatches)), function(m){
    a <- motifMatches[,m]; cq <- a - st$g[b,m]
    sum(a) * sum(a*cq*Eb[b])^2 / (sum(cq^2*Eb[b]) * sum(a*Eb[b]))
  })
  expect_equal(unname(betterChromVAR:::.effectiveN(motifMatches, bg, bg@expectation)),
               neff)
  expect_true(all(neff < colSums(motifMatches)))
  # without correction of the variance the z-scores are smaller
  expect_true(all(betterChromVAR:::.effectiveN(motifMatches, bg, bg@expectation, 
                                               fpc=FALSE) <= neff))
  # it is stored in the object, and used for the normalization
  dev <- betterChromVAR(counts, motifMatches)
  expect_equal(rowData(dev)$Neff, neff, ignore_attr=TRUE)
  dev <- normalizeDevsForSize(dev)
  expect_equal(assay(dev, "norm"), 
               assay(dev, "z")*sqrt(median(neff)/neff), ignore_attr=TRUE)
  # objects without Neff (e.g. from earlier versions) use N
  rowData(dev)$Neff <- NULL
  expect_equal(assay(normalizeDevsForSize(dev), "norm"), 
               assay(dev, "z")*sqrt(median(rowData(dev)$N)/rowData(dev)$N),
               ignore_attr=TRUE)
})

test_that("the finite-population correction is applied consistently", {
  dev0 <- betterChromVAR(counts, motifMatches, fpc=FALSE)
  dev1 <- betterChromVAR(counts, motifMatches)
  expect_identical(assay(dev0, "deviations"), assay(dev1, "deviations"))
  expect_false(isTRUE(all.equal(assay(dev0, "z"), assay(dev1, "z"))))
  # betterChromVAR is equivalent to the individual steps
  bg <- computeBackgrounds(counts, getBackgroundBins(counts, verbose=FALSE))
  expect_equal(assays(dev1), assays(computeDeviationsAnalytic(counts, bg,
                                                       motifMatches)))
  # as is the weighted version when applied to unweighted counts
  mc <- Matrix::t(motifMatches) %*% assay(counts)
  for(fpc in c(TRUE, FALSE)){
    dw <- computeDeviationsWeighted(mc, counts, motifMatches, fpc=fpc,
                                    verbose=FALSE)
    expect_equal(assay(dw, "z"), assay(betterChromVAR(counts, motifMatches,
                                                      fpc=fpc), "z"))
  }
})

test_that("estimateNoiseKappa works", {
  res <- estimateNoiseKappa(counts, motifMatches, pairs=FALSE)
  expect_s3_class(res, "data.frame")
  expect_equal(colnames(res), c("kappa", "reliability", "reliability.se"))
  expect_equal(rownames(res), colnames(motifMatches))
  expect_true(all(is.finite(res$kappa)) && all(res$kappa>0))
  expect_true(all(res$reliability>=0 & res$reliability<=1))
  kappa <- res
  # reproducible when setting the seed
  set.seed(1)
  k1 <- estimateNoiseKappa(counts, motifMatches, pairs=FALSE, chunkSize=3)
  set.seed(1)
  k2 <- estimateNoiseKappa(counts, motifMatches, pairs=FALSE, chunkSize=3)
  expect_identical(k1, k2)
  set.seed(2)
  expect_false(identical(k1, estimateNoiseKappa(counts, motifMatches, 
                                                pairs=FALSE, chunkSize=3)))
  # also reproducible when running in parallel
  skip_on_os("windows")
  set.seed(1)
  k3 <- estimateNoiseKappa(counts, motifMatches, pairs=FALSE, chunkSize=3,
                           nthreads=2L)
  set.seed(1)
  k4 <- estimateNoiseKappa(counts, motifMatches, pairs=FALSE, chunkSize=3,
                           nthreads=2L)
  expect_identical(k3, k4)
  expect_equal(dim(k3), dim(k1))
  # wrong input
  expect_error(estimateNoiseKappa(counts, motifMatches, nthreads="a"))
  expect_error(estimateNoiseKappa(assay(counts), motifMatches))
})

test_that("paired insertions are detected and thinned in pairs", {
  se <- counts
  assay(se) <- assay(counts)*2L
  expect_message(estimateNoiseKappa(se, motifMatches), "paired")
  expect_no_message(estimateNoiseKappa(counts, motifMatches, pairs=FALSE))
  y <- as(as(assay(se), "dMatrix"), "CsparseMatrix")
  sp <- betterChromVAR:::.thinCounts(y, 0.5, pairs=TRUE)
  expect_equal(as.matrix(sp$A + sp$B), as.matrix(y), ignore_attr=TRUE)
  expect_true(all(as.matrix(sp$A)%%2==0))
  expect_true(all(as.matrix(sp$B)%%2==0))
  sp <- betterChromVAR:::.thinCounts(y, 0.5, pairs=FALSE)
  expect_equal(as.matrix(sp$A + sp$B), as.matrix(y), ignore_attr=TRUE)
  expect_false(all(as.matrix(sp$A)%%2==0))
})

test_that("the noise calibration is reasonable for sparse Poisson counts", {
  set.seed(42)
  np <- 3000; nc <- 300
  mu <- rgamma(np, shape=2, scale=0.05)
  cnt <- matrix(rpois(np*nc, rep(mu, nc)), nrow=np)
  cnt <- cnt[rowSums(cnt)>0,]
  np <- nrow(cnt)
  se <- SummarizedExperiment(list(counts=cnt),
                             rowData=DataFrame(bias=runif(np, 0.3, 0.7)))
  ann <- Matrix(sample(c(0,1), np*6, replace=TRUE, prob=c(0.9, 0.1)),
                nrow=np, sparse=TRUE)
  colnames(ann) <- paste0("m", 1:6)
  kappa <- estimateNoiseKappa(se, ann, pairs=FALSE, bs=10)$kappa
  expect_true(all(kappa > 0.6 & kappa < 1.4))
})

test_that("reliability identifies motifs with real cell-level variation", {
  set.seed(7)
  np <- 4000; nc <- 600
  mu <- rgamma(np, shape=2, scale=0.05)
  ann <- Matrix(sample(c(0,1), np*6, replace=TRUE, prob=c(0.9, 0.1)),
                nrow=np, sparse=TRUE)
  colnames(ann) <- paste0("m", 1:6)
  # the first motif is more accessible in half of the cells
  up <- sample(c(FALSE, TRUE), nc, replace=TRUE)
  rate <- matrix(rep(mu, nc), nrow=np)
  rate[ann[,1]>0, up] <- rate[ann[,1]>0, up]*1.5
  cnt <- matrix(rpois(np*nc, rate), nrow=np)
  keep <- rowSums(cnt)>0
  se <- SummarizedExperiment(list(counts=cnt[keep,]),
                             rowData=DataFrame(bias=runif(sum(keep), 0.3, 0.7)))
  res <- estimateNoiseKappa(se, ann[keep,], pairs=FALSE, bs=10, ncells=600)
  expect_gt(res$reliability[1], 0.4)
  expect_lt(max(res$reliability[-1]), 0.25)
  expect_gt(res$reliability[1], max(res$reliability[-1]) + 2*res$reliability.se[1])
  # unweighted version also works
  res2 <- estimateNoiseKappa(se, ann[keep,], pairs=FALSE, bs=10, ncells=600,
                             cellWeights=FALSE)
  expect_gt(res2$reliability[1], max(res2$reliability[-1]))
  expect_error(estimateNoiseKappa(se, ann[keep,], pairs=FALSE, cellWeights="x"))
})

test_that("smoothDeviationsKNN works", {
  set.seed(1)
  nC <- 120
  nM <- 8
  grp <- rep(1:4, each=30)
  truth <- matrix(rnorm(nM*4, sd=0.4), nM, 4)[,grp] +
    matrix(rnorm(nM*nC, sd=0.3), nM, nC)
  sdn <- matrix(runif(nM*nC, 0.2, 0.6), nM, nC)
  dv <- truth + matrix(rnorm(nM*nC), nM, nC)*sdn
  se <- SummarizedExperiment(list(deviations=dv, z=dv/sdn))
  coords <- cbind(grp, 0) + matrix(rnorm(nC*2, sd=0.1), nC, 2)
  knn <- BiocNeighbors::findKNN(coords, k=15)

  sm <- smoothDeviationsKNN(se, knn)
  expect_true(all(c("deviations", "z", "smoothed", "smoothed.z") %in%
                    assayNames(sm)))
  expect_equal(dim(assay(sm, "smoothed")), dim(dv))
  expect_false(anyNA(assay(sm, "smoothed")))
  expect_false(anyNA(assay(sm, "smoothed.z")))
  expect_true("tau2" %in% colnames(rowData(sm)))
  expect_true(all(rowData(sm)$tau2 >= 0))
  # smaller error to the truth than the raw deviations
  expect_lt(mean((assay(sm, "smoothed")-truth)^2), 0.75*mean((dv-truth)^2))
  # the matrix of indices is also accepted
  expect_equal(assay(smoothDeviationsKNN(se, knn$index), "smoothed"),
               assay(sm, "smoothed"))
  # a cell being its own neighbor is ignored
  expect_equal(assay(smoothDeviationsKNN(se, cbind(seq_len(nC), knn$index)),
                     "smoothed"), assay(sm, "smoothed"))
  # k is taken into account
  sm5 <- smoothDeviationsKNN(se, knn, k=5)
  expect_false(isTRUE(all.equal(assay(sm5, "smoothed"), assay(sm, "smoothed"))))
  # smoothing is reduced when the estimated biological variance is increased
  sm2 <- smoothDeviationsKNN(se, knn, tauScale=1e6)
  expect_lt(mean((assay(sm2, "smoothed")-dv)^2), mean((assay(sm, "smoothed")-dv)^2))
  expect_equal(assay(sm2, "smoothed"), dv, tolerance=1e-3)
  # cells with a lower noise keep more of their own value
  smw <- smoothDeviationsKNN(se, knn, retWeights=TRUE)
  w <- assay(smw, "own.weight")
  expect_true(all(w>=0 & w<=1))
  expect_gt(cor(as.vector(w), -as.vector(sdn)), 0)
  # kappa
  expect_error(smoothDeviationsKNN(se, knn, kappa=rep(1,3)))
  expect_error(smoothDeviationsKNN(se, knn, kappa=-1))
  kp <- runif(nM, 0.5, 1)
  sm <- smoothDeviationsKNN(se, knn, kappa=kp)
  expect_false(anyNA(assay(sm, "smoothed")))
  expect_equal(assay(smoothDeviationsKNN(se, knn, kappa=data.frame(kappa=kp)),
                     "smoothed"), assay(sm, "smoothed"))
  # blocks do not affect the results
  expect_equal(assay(smoothDeviationsKNN(se, knn, blockSize=3), "smoothed"),
               assay(smoothDeviationsKNN(se, knn), "smoothed"))
  # invalid input
  expect_error(smoothDeviationsKNN(se, knn$index[1:10,]))
  expect_error(smoothDeviationsKNN(se, knn, k=100))
  expect_error(smoothDeviationsKNN(se[,1:10], knn))
  expect_error(smoothDeviationsKNN(SummarizedExperiment(list(deviations=dv)), knn))
})

test_that("smoothDeviationsKNN handles undefined noise values", {
  set.seed(2)
  nC <- 60
  dv <- matrix(rnorm(3*nC, sd=0.3), 3, nC)
  z <- dv/0.3
  z[1,5] <- 0         # undefined noise variance
  z[2,7] <- Inf       # zero noise variance
  se <- SummarizedExperiment(list(deviations=dv, z=z))
  knn <- BiocNeighbors::findKNN(matrix(rnorm(nC*2), nC, 2), k=10)
  sm <- smoothDeviationsKNN(se, knn)
  expect_false(anyNA(assay(sm, "smoothed")))
  expect_true(all(is.finite(assay(sm, "smoothed.z"))))
})

test_that("computeMotifVariability gives the reliability", {
  set.seed(3)
  nC <- 400
  nM <- 4
  # the first motif varies between cells (sd 0.5), the others do not
  sdn <- 0.4
  truth <- matrix(0, nM, nC)
  truth[1,] <- rnorm(nC, sd=0.5)
  dv <- truth + matrix(rnorm(nM*nC, sd=sdn), nM, nC)
  rownames(dv) <- paste0("m", 1:nM)
  r <- computeMotifVariability(dv/sdn)
  expect_true(all(r$reliability>=0 & r$reliability<=1))
  expect_equal(r$reliability, pmax(0, 1-1/r$variability^2))
  # true reliability: 0.25/(0.25+0.16)
  expect_equal(r$reliability[1], 0.25/(0.25+sdn^2), tolerance=0.1)
  expect_lt(max(r$reliability[-1]), 0.2)
  # is also in the output of the main functions
  dev <- betterChromVAR(counts, motifMatches)
  expect_true("reliability" %in% colnames(rowData(dev)))
})

test_that("computeDeviationsWeighted is consistent with the unweighted version", {
  bg <- computeBackgrounds(counts, getBackgroundBins(counts, verbose=FALSE))
  mc <- as.matrix(Matrix::t(motifMatches) %*% assay(counts))
  ref <- computeDeviationsAnalytic(counts, bg, motifMatches)
  w <- computeDeviationsWeighted(mc, counts, motifMatches, bg=bg)
  expect_equal(assay(w, "z"), assay(ref, "z"))
  expect_equal(assay(w, "deviations"), assay(ref, "deviations"))
  expect_equal(rowData(w)$Neff, rowData(ref)$Neff)
  # neither the z-scores nor the deviations depend on the magnitude of the 
  # weights
  w2 <- computeDeviationsWeighted(2.5*mc, counts, motifMatches, bg=bg)
  expect_equal(assay(w2, "z"), assay(ref, "z"))
  expect_equal(assay(w2, "deviations"), assay(ref, "deviations"))

})

.vr <- function(se) as.data.frame(rowData(computeMotifVariability(se)))

test_that("the significance of the variability is calibrated and powerful", {
  set.seed(11)
  nC <- 300
  nM <- 400
  # noise only, with noise variances varying between cells and motifs, and 
  # heavy-tailed noise
  sdn <- matrix(runif(nM*nC, 0.2, 1.5), nM, nC) * rep(exp(rnorm(nC, sd=0.5)), each=nM)
  noise <- matrix(rt(nM*nC, df=8)/sqrt(8/6), nM, nC)*sdn
  se0 <- SummarizedExperiment(list(deviations=noise, z=noise/sdn))
  rel0 <- .vr(se0)
  expect_true(all(rel0$var.pval >= 0 & rel0$var.pval <= 1))
  # roughly uniform p-values under the null hypothesis
  expect_lt(mean(rel0$var.pval < 0.05), 0.09)
  expect_gt(mean(rel0$var.pval < 0.05), 0.02)
  expect_lt(mean(rel0$var.adjPval < 0.05), 0.02)
  # motifs whose activity truly varies are detected, in proportion to their 
  # reliability
  truth <- matrix(0, nM, nC)
  truth[1:20,] <- rnorm(20*nC, sd=1)
  truth[21:40,] <- rnorm(20*nC, sd=0.5)
  truth[41:60,] <- rnorm(20*nC, sd=0.25)
  se1 <- SummarizedExperiment(list(deviations=noise+truth, 
                                   z=(noise+truth)/sdn))
  rel1 <- .vr(se1)
  expect_true(all(rel1$var.adjPval[1:40] < 1e-3))
  expect_gt(mean(rel1$var.adjPval[41:60] < 0.05), 0.2)
  expect_lt(mean(rel1$var.adjPval[-(1:60)] < 0.05), 0.02)
  expect_gt(mean(rel1$reliability[1:20]), mean(rel1$reliability[21:40]))
  expect_gt(mean(rel1$reliability[21:40]), mean(rel1$reliability[41:60]))
})
