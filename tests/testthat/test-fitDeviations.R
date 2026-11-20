test_that("fitDeviations fits a limma model on the deviations", {
  dd <- getDummyData()
  dev <- betterChromVAR(dd$counts, dd$motifMatches)
  dev$group <- rep(LETTERS[1:2], each=5)
  design <- model.matrix(~group, data=colData(dev))
  fit <- fitDeviations(dev, design)
  expect_s4_class(fit, "MArrayLM")
  expect_equal(nrow(fit), nrow(dev))
  expect_equal(colnames(fit$coefficients), colnames(design))
  # the default uses the size-normalized z-scores
  expect_equal(fit$coefficients, 
               limma::lmFit(assay(normalizeDevsForSize(dev), "norm"), 
                            design)$coefficients)
  # a formula is evaluated on the colData
  expect_equal(fitDeviations(dev, ~group), fit)
  fit2 <- fitDeviations(dev, design, type="deviations")
  expect_equal(fit2$coefficients, 
               limma::lmFit(assay(dev, "deviations"), design)$coefficients)
  # contrasts
  cm <- matrix(c(0,1), ncol=1, dimnames=list(colnames(design), "B-A"))
  fit3 <- fitDeviations(dev, design, contrasts=cm)
  expect_equal(colnames(fit3$coefficients), "B-A")
  expect_error(fitDeviations(dev[,1:3], design))
})

test_that("fitDeviations uses a trend with log(Neff) when there are enough motifs", {
  data("NR3C1example", package="betterChromVAR", envir=environment())
  dev <- betterChromVAR(peakSE, motifMatches)
  expect_gte(nrow(dev), 20)
  design <- model.matrix(~group, data=colData(dev))
  fit <- fitDeviations(dev, design)
  ref <- limma::eBayes(limma::lmFit(assay(normalizeDevsForSize(dev), "norm"), 
                                    design), 
                       trend=log(pmax(rowData(dev)$Neff, 1)))
  expect_equal(fit$s2.prior, ref$s2.prior)
  expect_equal(fit$p.value, ref$p.value)
  # not with a small number of motifs (as there is not enough to fit a trend)
  fit2 <- fitDeviations(dev[1:10,], design)
  expect_length(fit2$s2.prior, 1L)
  expect_length(fit$s2.prior, nrow(dev))
})
