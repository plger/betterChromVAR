library(Matrix)
library(SummarizedExperiment)
set.seed(123)

# dummy data with two overlapping motifs (motif2 shares half of the peaks of 
# motif1), an unrelated one, and an effect on the peaks of motif1 in half of 
# the samples
dat <- getDummyData(nRegions=3000, nSamples=20, nMotifs=3)
counts <- dat$counts
mm <- as.matrix(dat$motifMatches)
i1 <- which(mm[,1]==1)
mm[,2] <- 0
mm[c(i1[seq_len(floor(length(i1)/2))], 
     sample(which(mm[,1]==0), floor(length(i1)/2))), 2] <- 1
motifMatches <- as(Matrix(mm, sparse=TRUE), "CsparseMatrix")
null <- counts
grp <- rep(1:2, each=10)
assay(counts)[i1, grp==2] <- 2L*assay(counts)[i1, grp==2]

test_that("no ridge, or an infinite one, gives the standard results", {
  bg <- computeBackgrounds(counts, getBackgroundBins(counts, verbose=FALSE))
  d0 <- computeDeviationsAnalytic(counts, bg, motifMatches)
  d1 <- computeDeviationsAnalytic(counts, bg, motifMatches, ridge=Inf)
  expect_identical(assay(d0, "deviations"), assay(d1, "deviations"))
  expect_identical(assay(d0, "z"), assay(d1, "z"))
  # a very large ridge converges to it
  d2 <- computeDeviationsAnalytic(counts, bg, motifMatches, ridge=1e8)
  expect_equal(assay(d0, "deviations"), assay(d2, "deviations"), 
               tolerance=1e-6)
  expect_equal(assay(d0, "z"), assay(d2, "z"), tolerance=1e-6)
  expect_error(computeDeviationsAnalytic(counts, bg, motifMatches, ridge=-1))
  expect_error(computeDeviationsAnalytic(counts, bg, motifMatches, ridge=1:2))
})

test_that("the motif overlap matches a peak-level computation", {
  bg <- computeBackgrounds(counts, getBackgroundBins(counts, verbose=FALSE))
  O <- betterChromVAR:::.motifOverlap(motifMatches, bg)
  st <- betterChromVAR:::.fpcStats(motifMatches, bg)
  # sum_q e_q c_qi c_qj, with c_q = a_q - g_b(q)
  cq <- as.matrix(motifMatches) - st$g[bg@peak2bin,]
  expect_equal(O, crossprod(cq, cq*bg@expectation), ignore_attr=TRUE)
  # dense annotations are processed by blocks of rows
  x <- motifMatches
  x@x <- runif(length(x@x))
  e <- bg@expectation
  expect_equal(betterChromVAR:::.weightedCrossprod(x, e, maxDensity=0, 
                                                   blockCells=1000),
               as.matrix(crossprod(x, x*e)), ignore_attr=TRUE)
  expect_equal(betterChromVAR:::.weightedCrossprod(as.matrix(x), e),
               as.matrix(crossprod(x, x*e)), ignore_attr=TRUE)
  op <- betterChromVAR:::.ridgeOperator(motifMatches, bg, 1)
  expect_equal(dim(op$L), rep(ncol(motifMatches), 2))
  expect_true(all(op$varRatio>0))
  # without penalty, the excess counts produced by a change in the peaks of 
  # one motif (a column of O) are entirely attributed to that motif
  op <- betterChromVAR:::.ridgeOperator(motifMatches, bg, 1e-10)
  expect_equal(op$L %*% O, diag(diag(O)), ignore_attr=TRUE, tolerance=1e-6)
})

test_that("the ridge reduces the spill-over to overlapping motifs", {
  d0 <- betterChromVAR(counts, motifMatches)
  d1 <- betterChromVAR(counts, motifMatches, ridge=0.1)
  eff <- function(d){
    x <- assay(d, "deviations")
    rowMeans(x[,grp==2]) - rowMeans(x[,grp==1])
  }
  e0 <- eff(d0)
  e1 <- eff(d1)
  # the overlapping motif has a clear spurious effect, which is mostly removed
  expect_gt(e0[2], 0.1)
  expect_lt(abs(e1[2]), abs(e0[2])/3)
  # while the effect on the true motif is retained
  expect_gt(e1[1], 0.8*e0[1])
  # and the same holds for the z-scores
  zeff <- function(d){
    x <- assay(d, "z")
    rowMeans(x[,grp==2]) - rowMeans(x[,grp==1])
  }
  expect_lt(abs(zeff(d1)[2]), abs(zeff(d0)[2])/3)
  expect_gt(zeff(d1)[1], 0.7*zeff(d0)[1])
})

test_that("the z-scores remain on the same scale in the absence of signal", {
  d0 <- betterChromVAR(null, motifMatches)
  d1 <- betterChromVAR(null, motifMatches, ridge=0.3)
  r <- matrixStats::rowSds(assay(d1, "z"))/matrixStats::rowSds(assay(d0, "z"))
  expect_true(all(r > 0.7 & r < 1.4))
})

test_that("the ridge gives the same results when cells are split in chunks", {
  big <- getDummyData(nRegions=1000, nSamples=120, nMotifs=4)
  d1 <- betterChromVAR(big$counts, big$motifMatches, ridge=1)
  d2 <- betterChromVAR(big$counts, big$motifMatches, ridge=1, nthreads=2)
  expect_equal(assay(d1, "deviations"), assay(d2, "deviations"))
  expect_equal(assay(d1, "z"), assay(d2, "z"))
})

test_that("scores are graded as expected", {
  lenient <- sparseMatrix(i=c(1,2,3,4, 1,2), j=c(1,1,1,1, 2,2), 
                          x=c(4,6,8,10, 5,7), dims=c(5,3))
  stringent <- sparseMatrix(i=c(3,4), j=c(1,1), x=c(8,10), dims=c(5,3))
  g <- betterChromVAR:::.gradeScores(lenient, stringent, minWeight=0)
  expect_equal(as.matrix(g)[,1], c(0, 0.5, 1, 1, 0))
  # no stringent match: no weight
  expect_true(all(as.matrix(g)[,2:3]==0))
  g <- betterChromVAR:::.gradeScores(lenient, stringent, minWeight=0.6)
  expect_equal(as.matrix(g)[,1], c(0, 0, 1, 1, 0))
  # identical thresholds give binary matches
  g <- betterChromVAR:::.gradeScores(stringent, stringent)
  expect_equal(as.matrix(g)[,1], c(0, 0, 1, 1, 0))
})

test_that("graded matches can be obtained from motifmatchr", {
  skip_if_not_installed("motifmatchr")
  skip_if_not_installed("TFBSTools")
  pfm <- TFBSTools::PFMatrix(ID="m1", name="m1", profileMatrix=matrix(
    c(18,1,1,0, 0,0,20,0, 20,0,0,0, 0,0,0,20, 19,0,1,0, 17,1,1,1), nrow=4,
    dimnames=list(c("A","C","G","T"), NULL)))
  pwms <- TFBSTools::PWMatrixList(m1=TFBSTools::toPWM(pfm))
  seqs <- vapply(seq_len(2000), FUN.VALUE=character(1), function(i)
    paste(sample(c("A","C","G","T"), 200, replace=TRUE), collapse=""))
  seqs <- Biostrings::DNAStringSet(seqs)
  bin <- motifmatchr::motifMatches(motifmatchr::matchMotifs(pwms, seqs,
                                                            p.cutoff=5e-5))
  g <- getGradedMatches(pwms, seqs, pCutoffs=c(5e-3, 5e-5))
  expect_equal(dim(g), dim(bin))
  expect_equal(colnames(g), "m1")
  expect_true(all(g@x > 0 & g@x <= 1))
  # stringent matches have a weight of 1, and there are additional weak ones
  expect_true(all(g[which(bin[,1]),1]==1))
  expect_gt(sum(g[,1]>0), sum(bin[,1]))
  expect_true(all(g[which(!bin[,1]),1] < 1))
})

test_that("annotations that are not sparse are processed by blocks", {
  x <- motifMatches
  x@x <- runif(length(x@x))
  sparseCounts <- as(Matrix(assay(counts), sparse=TRUE), "CsparseMatrix")
  for(cnts in list(assay(counts), sparseCounts)){
    ref <- as.matrix(crossprod(x, cnts))
    expect_equal(betterChromVAR:::.annotCrossprod(x, cnts, maxDensity=0,
                                                  blockCells=1000), ref)
    expect_equal(betterChromVAR:::.annotCrossprod(x, cnts), ref)
    expect_equal(betterChromVAR:::.annotCrossprod(as.matrix(x), cnts), ref)
  }
  # logical matches
  expect_equal(betterChromVAR:::.annotCrossprod(as(motifMatches, "lMatrix"), assay(counts),
                                                maxDensity=0),
               as.matrix(crossprod(motifMatches, assay(counts))))
  # a dense annotation gives the same results whether given as a sparse (block
  # processing) or as a base (direct) matrix
  dense <- Matrix(runif(nrow(counts)*4), ncol=4, sparse=TRUE,
                  dimnames=list(NULL, paste0("m", 1:4)))
  dense <- as(dense, "CsparseMatrix")
  expect_false(betterChromVAR:::.isSparseAnnotation(dense))
  d1 <- betterChromVAR(counts, dense, ridge=3)
  d2 <- betterChromVAR(counts, as.matrix(dense), ridge=3)
  expect_equal(assay(d1, "deviations"), assay(d2, "deviations"))
  expect_equal(assay(d1, "z"), assay(d2, "z"))
  expect_equal(rowData(d1)$Neff, rowData(d2)$Neff)
})
