library(BiocParallel)
library(Matrix)
set.seed(123)
attach(getDummyData())
counts$groups <- factor(rep(LETTERS[1:2], each=5))

checkDevOutput <- function(dev, checkAssayNames=TRUE){
  expect_true(is(dev, "SummarizedExperiment"))
  if(checkAssayNames){
    expect_all_true(c("deviations", "z") %in% assayNames(dev))
  }
  expect_false(any(is.na(assay(dev))))
}

test_that("betterChromVAR runs", {
  dev <- betterChromVAR(counts, motifMatches)
  checkDevOutput(dev)
})

test_that("dgcMatrix supported", {
  assay(counts) <- as(assay(counts), "sparseMatrix")
  dev <- betterChromVAR(counts, motifMatches)
  checkDevOutput(dev)
})


bg <- getBackgroundBins(counts)

test_that("multithreading works", {
  bp <- SnowParam(2)
  dev <- betterChromVAR(counts, motifMatches, grouping=counts$groups,
                        nthreads=bp)
  checkDevOutput(dev)
})

test_that("CVnorm works", {
  counts <- CVnorm(counts, grouping=counts$groups)
  checkDevOutput(counts, FALSE)
  expect_true("corrected" %in% assayNames(counts))
})

test_that("Fragment length bias works", {
  rowData(counts)$flbias <- pmax(rnorm(nrow(counts), 2.5, 0.25),0.5)
  background <- getBackgroundBins(counts, bs=c(10,10,4))
  expect_equal(nrow(background@binBinProbs), (10*10*4))
})

test_that("Sampling background peaks works", {
  background <- getBackgroundBins(counts)
  bg_peaks <- sampleBackgroundPeaks(background, niterations=10)
  expect_all_true(dim(bg_peaks)==c(nrow(counts), 10L))
})


set.seed(123)
z <- matrix(rnorm(mean=rnorm(10), sd=runif(10, max=2), 1000), nrow=10)

compareVarDat <- function(x,y,isP=FALSE){
  if(isP) return(1-cor(x,y))
  mean(abs(x-y)/x)
}

test_that("computeMotifVariability works", {
  v0 <- computeMotifVariability(z, confInt = 0.6)
  v1 <- computeMotifVariability(z, confInt = 0.6, method="normal")
  v2 <- computeMotifVariability(z, confInt = 0.6, method="bootstrap")
  diff <- vapply(colnames(v0)[1:4], FUN.VALUE=numeric(2), \(x){
    c(compareVarDat(v0[[x]], v1[[x]], grepl("pval", x)),
      compareVarDat(v0[[x]], v2[[x]], grepl("pval", x)))
  })
  expect_all_true(as.numeric(diff)<0.05)
})

