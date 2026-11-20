local({
  data("NR3C1example", package="betterChromVAR", envir=environment())
  idev <- getInteractionsDeviations(peakSE, motifMatches, bait="GCR",
                                    minCount=10)

  test_that("getInteractionsDeviations returns a SummarizedExperiment", {
    expect_s4_class(idev, "SummarizedExperiment")
    expect_true(all(c("z","deviations","sd") %in% assayNames(idev)))
    expect_true(all(is.finite(assay(idev, "sd"))) && all(assay(idev, "sd")>0))
    expect_equal(sum(rowData(idev)$isBait), 1L)
    expect_equal(row.names(idev)[rowData(idev)$isBait], "GCR")
  })

  test_that("getInteractionsDeviations defines disjoint sets", {
    rd <- rowData(idev)
    expect_equal(sort(unique(as.character(rd$set))),
                 c("bait","baitOnly","both","partnerOnly"))
    p <- "JUND"
    a <- assay(motifMatches)>0
    n <- rd$N[match(paste0(p, c(".both",".baitOnly",".partnerOnly")),
                    row.names(idev))]
    expect_equal(n, c(sum(a[,"GCR"] & a[,p]), sum(a[,"GCR"] & !a[,p]),
                      sum(!a[,"GCR"] & a[,p])))
    # null sd is dev/z
    i <- which(row.names(idev)=="JUND.both")
    expect_equal(assay(idev,"sd")[i,],
                 abs(assay(idev,"deviations")[i,]/assay(idev,"z")[i,]))
  })

  test_that("getInteractionsDeviations adds overlap metadata columns", {
    rd <- rowData(idev)
    expect_true(all(c("overlap.percentOfBait","overlap.jaccard",
                      "overlap.log2Enr","partner") %in% colnames(rd)))
    expect_equal(unique(rd$overlap.percentOfBait[rd$partner=="JUND"]),
                 round(100*rd["JUND.both","N"]/rd["GCR","N"], 2))
  })

  test_that("getInteractionsDeviations errors on bad input", {
    expect_error(getInteractionsDeviations(peakSE, motifMatches,
                                           bait="notAMotif"))
    expect_error(getInteractionsDeviations(peakSE, motifMatches, bait="GCR",
                                           minCount=1e5))
  })

  test_that("discoverMotifInteractions returns the expected data.frame", {
    res <- discoverMotifInteractions(idev, "group")
    expect_s3_class(res, "data.frame")
    expect_true(all(c("contrast","motif","difference","P.Value","adj.P.Val",
                      "overlap.N","diff.both","diff.baitOnly",
                      "diff.partnerOnly") %in% colnames(res)))
    expect_true(all(res$adj.P.Val >= 0 & res$adj.P.Val <= 1))
    expect_false(is.unsorted(res$P.Value))
    expect_setequal(res$motif, setdiff(unique(rowData(idev)$partner), "GCR"))
    # the interaction is (approximately, since the sets are weighted 
    # differently) the difference between the response of the sets
    expect_gt(cor(res$difference,
                  res$diff.both-res$diff.baitOnly-res$diff.partnerOnly), 0.99)
    # the response of the bait's own sites is positive (GCR activation)
    expect_true(median(res$diff.baitOnly)>0)
  })

  test_that("discoverMotifInteractions works with global tests, non-syntactic
            groups, covariates and without weights", {
    expect_s3_class(discoverMotifInteractions(idev, "group", global=TRUE),
                    "data.frame")
    i2 <- idev
    i2$g <- c("0h","0h","0h","4h-DEX","4h-DEX","4h-DEX")
    expect_equal(unique(discoverMotifInteractions(i2, "g")$contrast),
                 "4h-DEX-0h")
    i2$g3 <- c("a","a","b","b","c","c")
    i2$cv <- c(1,2,3,2,1,2)
    res <- discoverMotifInteractions(i2, "g3", covar="cv", global=TRUE)
    expect_true("F" %in% colnames(res))
    expect_equal(length(grep("^difference\\.", colnames(res))), 2L)
    res <- discoverMotifInteractions(i2, "g3", weights="none")
    expect_equal(length(unique(res$contrast)), 2L)
  })

  test_that("discoverMotifInteractions rejects objects from older versions", {
    old <- idev
    assay(old, "sd") <- NULL
    expect_error(discoverMotifInteractions(old, "group"), "getInteractionsDev")
  })

  test_that("peakLevelMotifInteractions returns the expected data.frame", {
    res <- peakLevelMotifInteractions(peakSE, motifMatches, bait="GCR",
                                      group="group", minCount=10)
    expect_s3_class(res, "data.frame")
    expect_true(all(c("contrast","motif","difference","se","P.Value",
                      "adj.P.Val","bait.diff","partner.diff","tau.both",
                      "overlap.N") %in% colnames(res)))
    expect_true(all(res$se>0) && all(res$P.Value>=0 & res$P.Value<=1))
    expect_false(is.unsorted(res$P.Value))
    expect_true(median(res$bait.diff)>0)
    expect_error(peakLevelMotifInteractions(peakSE, motifMatches, bait="GCR",
                                      group="group", minCount=1e5))
    # without the variance inflation the errors can only be smaller
    r2 <- peakLevelMotifInteractions(peakSE, motifMatches, bait="GCR",
                                     group="group", minCount=10, vif=FALSE)
    r2 <- r2[match(res$motif, r2$motif),]
    expect_true(all(r2$se<=res$se+1e-12))
    expect_equal(r2$difference, res$difference)
  })

  test_that("peakLevelMotifInteractions recovers a planted interaction", {
    set.seed(1)
    a <- as.matrix(assay(motifMatches)>0)
    pm <- as(cbind(a[,c("GCR","JUND")],
                   rnd=runif(nrow(a))<0.1), "CsparseMatrix")
    se <- peakSE
    cn <- assay(se)
    treated <- se$group==levels(factor(se$group))[2]
    w <- which(a[,"GCR"] & pm[,"rnd"])
    cn[w,treated] <- rpois(length(w)*sum(treated),
                           as.vector(cn[w,treated])*3)
    assay(se) <- cn
    res <- peakLevelMotifInteractions(se, pm, bait="GCR", group="group",
                                      minCount=10)
    expect_equal(res$motif[1], "rnd")
    expect_true(res$difference[1]>1)
    expect_true(res$adj.P.Val[1]<0.01)
  })

  test_that("motifCoOccurence distances and regions are right", {
    skip_if_not_installed("TFBSTools")
    skip_if_not_installed("Rsamtools")
    set.seed(1)
    # (C/G-only background, so that the motifs only match where planted)
    sq <- paste(sample(c("C","G"), 20000, TRUE), collapse="")
    mA <- "ACGTGCAT"
    mB <- "TTGGCCAA"
    pos <- list(c(1050,mA),c(1100,mB), c(3050,mA),c(3250,mB), c(5100,mA),
                c(7010,mA),c(7018,mB), c(9100,mB), c(12100,mA),c(12400,mB),
                c(13900,mA),c(14500,mB))
    for(p in pos) substr(sq, as.integer(p[1]), as.integer(p[1])+7) <- p[2]
    fa <- tempfile(fileext=".fa")
    writeLines(c(">chr1", sq), fa)
    Rsamtools::indexFa(fa)
    pfm <- function(x, id){
      m <- sapply(strsplit(x,"")[[1]], \(b) 1000*(c("A","C","G","T")==b))
      rownames(m) <- c("A","C","G","T")
      TFBSTools::PFMatrix(ID=id, name=id, profileMatrix=m)
    }
    mot <- TFBSTools::PFMatrixList(A=pfm(mA,"A"), B=pfm(mB,"B"))
    names(mot) <- c("A","B")
    # regions: 1: d=50, 2: d=200, 3: A only, 4: d=8, 5: B only, 
    # 6: B is outside (400bp, the region being 300bp-wide), 
    # 7: has no match, but lies between an A and a B match
    regions <- GRanges("chr1", IRanges(c(1000,3000,5000,7000,9000,12000,14000),
                                       width=300))
    w <- \(x) lapply(x, \(m) which(as.matrix(m)[,1]))
    pr <- list(c("A","B"))
    
    # fixed thresholds
    expect_equal(w(motifCoOccurence(mot, pr, regions, fa, minDist=5, 
                                    maxDist=300, nDistQuantiles=NULL))[[1]],
                 c(1,2,4,6))
    expect_equal(w(motifCoOccurence(mot, pr, regions, fa, minDist=10, 
                                    maxDist=300, nDistQuantiles=NULL))[[1]],
                 c(1,2,6))
    # (region 7 is between the matches, but contains none)
    expect_equal(w(motifCoOccurence(mot, pr, regions, fa, minDist=0, 
                                    maxDist=800, nDistQuantiles=NULL))[[1]],
                 c(1,2,4,6))
    expect_equal(w(motifCoOccurence(mot, pr, regions, fa, minDist=5, 
                                    maxDist=100, nDistQuantiles=NULL))[[1]],
                 c(1,4))
    
    # distances between the closest edges instead of the centers (matches are
    # 8bp-wide: gaps of 42, 192, 0 and 292 in regions 1, 2, 4 and 6)
    ed <- \(minD, maxD) w(motifCoOccurence(mot, pr, regions, fa, minDist=minD, 
                      maxDist=maxD, centerDist=FALSE, nDistQuantiles=NULL))[[1]]
    expect_equal(ed(5, 300), c(1,2,6))
    expect_equal(ed(0, 45), c(1,4))
    expect_equal(ed(0, 40), 4)
    expect_equal(ed(5, 100), 1)
    
    # quantile bins, using the shortest distance in each region
    q <- motifCoOccurence(mot, pr, regions, fa, minDist=5, maxDist=500, 
                          nDistQuantiles=2)
    expect_equal(unname(attr(q, "breaks")[[1]]), c(8,125,300))
    expect_equal(w(q), list(bin1=c(1,4), bin2=c(2,6)))
    # minDist excludes the close pair
    q <- motifCoOccurence(mot, pr, regions, fa, minDist=10, maxDist=500, 
                          nDistQuantiles=2)
    expect_equal(w(q), list(bin1=c(1,2), bin2=6))
    expect_equal(attr(q, "breaks")[[1]], c(50,200,300))
    # restricting to the regions removes the match outside of region 6
    q <- motifCoOccurence(mot, pr, regions, fa, minDist=5, maxDist=500, 
                          nDistQuantiles=2, restrictToRegions=TRUE)
    expect_equal(w(q), list(bin1=c(1,4), bin2=2))
    expect_equal(attr(q, "breaks")[[1]], c(8,50,200))
    
    expect_error(motifCoOccurence(mot, pr, regions, fa, minDist=-1))
  })
})

local({
  data("NR3C1example", package="betterChromVAR", envir=environment())
  
  # differential accessibility analysis of the peaks
  .dea <- function(se){
    design <- model.matrix(~group, data=colData(se))
    v <- limma::voom(as.matrix(assay(se)), design, plot=FALSE)
    fit <- limma::eBayes(limma::lmFit(v, design))
    list(design=design, 
         dea=limma::topTable(fit, coef=2, number=Inf, sort.by="none"))
  }
  d <- .dea(peakSE)
  
  test_that("peakLevelMotifTests returns the expected data.frame", {
    res <- peakLevelMotifTests(peakSE, d$dea, motifMatches, design=d$design,
                               minCount=10)
    expect_s3_class(res, "data.frame")
    expect_true(all(c("contrast","motif","N","difference","se","t","P.Value",
                      "adj.P.Val","tau","rho") %in% colnames(res)))
    expect_true(all(res$se>0) && all(res$P.Value>=0 & res$P.Value<=1))
    expect_false(is.unsorted(res$P.Value))
    a <- assay(motifMatches)
    expect_equal(res$N, unname(Matrix::colSums(a)[res$motif]))
    expect_error(peakLevelMotifTests(peakSE, d$dea, motifMatches, 
                                     design=d$design, minCount=1e5))
    # GCR and MCR sites increase the most with GCR activation
    expect_true(all(res$motif[1:2] %in% c("GCR","MCR")))
    expect_true(all(res$difference[1:2]>0))
    expect_true(any(res$rho>0))
    # without the variance inflation the errors can only be smaller
    r2 <- peakLevelMotifTests(peakSE, d$dea, motifMatches, minCount=10, 
                              vif=FALSE)
    r2 <- r2[match(res$motif, r2$motif),]
    expect_true(all(r2$se<=res$se+1e-12))
    expect_equal(r2$difference, res$difference)
  })
  
  test_that("peakLevelMotifTests recovers a planted effect and is calibrated", {
    set.seed(1)
    n <- nrow(peakSE)
    pm <- as(cbind(rnd=runif(n)<0.05, rnd2=runif(n)<0.05, rnd3=runif(n)<0.05),
             "CsparseMatrix")
    se <- peakSE
    cn <- assay(se)
    treated <- se$group==levels(factor(se$group))[2]
    w <- which(pm[,"rnd"])
    cn[w,treated] <- rpois(length(w)*sum(treated), 
                           as.vector(cn[w,treated])*2)
    assay(se) <- cn
    dd <- .dea(se)
    res <- peakLevelMotifTests(se, dd$dea, pm, design=dd$design, minCount=10)
    expect_equal(res$motif[1], "rnd")
    expect_true(res$difference[1]>0.5)
    expect_lt(res$adj.P.Val[1], 0.01)
    # the other (null) motifs show no effect
    expect_true(all(res$P.Value[res$motif!="rnd"]>0.01))
    expect_true(all(abs(res$difference[res$motif!="rnd"])<0.2))
  })
  
  test_that("peakLevelMotifTests handles the different formats of tables", {
    tb <- d$dea
    tb$bias <- rowData(peakSE)$bias
    ref <- peakLevelMotifTests(peakSE, tb, motifMatches, minCount=10, 
                               vif=FALSE)
    # the correlation requires the design
    expect_message(r0 <- peakLevelMotifTests(peakSE, tb, motifMatches, 
                                             minCount=10), "design")
    expect_equal(r0, ref)
    expect_no_message(peakLevelMotifTests(peakSE, tb, motifMatches, 
                                          minCount=10, vif=FALSE))
    r1 <- peakLevelMotifTests(peakSE, tb, motifMatches, design=d$design,
                              minCount=10)
    r1 <- r1[match(ref$motif, r1$motif),]
    expect_true(all(r1$se>=ref$se-1e-12) && any(r1$se>ref$se))
    # a fixed excess correlation
    r3 <- peakLevelMotifTests(peakSE, tb, motifMatches, minCount=10, vif=0.05)
    r3 <- r3[match(ref$motif, r3$motif),]
    expect_true(all(r3$se>=ref$se-1e-12) && any(r3$se>ref$se))
    expect_error(peakLevelMotifTests(peakSE, tb, motifMatches, vif=2))
    
    # the bias is taken from the table, or else from the object
    tb2 <- tb[,colnames(tb)!="bias"]
    expect_equal(peakLevelMotifTests(peakSE, tb2, motifMatches, minCount=10, 
                                     vif=FALSE), ref)
    expect_equal(unique(ref$contrast), "contrast")
    
    # other formats
    sdv <- abs(tb$logFC/tb$t)
    deseq <- data.frame(log2FoldChange=tb$logFC, lfcSE=sdv, 
                        baseMean=2^tb$AveExpr-1, row.names=row.names(tb))
    expect_equal(peakLevelMotifTests(peakSE, deseq, motifMatches, minCount=10,
                                     vif=FALSE)$difference[1:3], 
                 ref$difference[1:3], tolerance=0.05)
    tstat <- data.frame(logFC=tb$logFC, stat=tb$t, row.names=row.names(tb))
    # (the accessibility is then computed from the counts)
    expect_equal(peakLevelMotifTests(peakSE, tstat, motifMatches, minCount=10,
                                     vif=FALSE), ref)
    fstat <- data.frame(logFC=tb$logFC, F=tb$t^2, logCPM=tb$AveExpr,
                        row.names=row.names(tb))
    expect_equal(peakLevelMotifTests(peakSE, fstat, motifMatches, 
                                     minCount=10, vif=FALSE), ref)
    
    # several contrasts, and matching of the peaks through the row names
    set.seed(1)
    lst <- list(a=tb[sample(nrow(tb)),], b=tb)
    lst$a <- lst$a[order(row.names(lst$a)),]
    lst$b <- lst$a
    lst$b$logFC <- -lst$b$logFC
    lst$b$t <- -lst$b$t
    r4 <- peakLevelMotifTests(peakSE, lst, motifMatches, minCount=10, 
                              vif=FALSE)
    expect_equal(sort(unique(r4$contrast)), c("a","b"))
    gcr <- function(r, co) r$difference[r$contrast==co & r$motif=="GCR"]
    expect_equal(gcr(r4, "a"), -gcr(r4, "b"))
    expect_gt(gcr(r4, "a"), 0)
    # (the order of the peaks of the table doesn't matter)
    expect_equal(r4[r4$contrast=="a",-1], ref[,-1][match(r4$motif[r4$contrast=="a"], ref$motif),])
    
    # a table with missing values or peaks
    tb3 <- tb[-(1:100),]
    tb3$logFC[5] <- NA
    expect_message(r5 <- peakLevelMotifTests(peakSE, tb3, motifMatches, 
                                             minCount=10, vif=FALSE),
                   "101 peaks")
    expect_equal(r5$motif[1], ref$motif[1])
  })
  
  test_that("peakLevelMotifTests checks its input", {
    tb <- d$dea
    tstat <- data.frame(logFC=tb$logFC, stat=tb$t, row.names=row.names(tb))
    suppressMessages({
      expect_error(peakLevelMotifTests(peakSE, tstat[,"stat",drop=FALSE], 
                                       motifMatches), "log-fold-change")
      expect_error(peakLevelMotifTests(peakSE, tstat[,"logFC",drop=FALSE], 
                                       motifMatches), "standard error")
      t2 <- tstat
      row.names(t2) <- paste0("x", seq_len(nrow(t2)))
      expect_error(peakLevelMotifTests(peakSE, t2, motifMatches), "row names")
      expect_error(peakLevelMotifTests(peakSE[1:100,], tstat, motifMatches),
                   "same number of rows")
      expect_error(peakLevelMotifTests(peakSE, tstat, motifMatches, 
                                       design=d$design[1:3,]), "design")
      expect_error(peakLevelMotifTests(tstat, tstat, motifMatches))
    })
  })
  
  test_that("stepwisePeakLevelMotifTests removes the spillover to overlapping motifs", {
    set.seed(42)
    n <- 8000
    A <- runif(n)<0.08
    B <- (A & runif(n)<0.7) | (!A & runif(n)<0.025)  # shares 70% of its peaks
    C <- runif(n)<0.06
    D <- runif(n)<0.07
    ann <- cbind(A=A, B=B, C=C, D=D)
    row.names(ann) <- paste0("p", seq_len(n))
    # (placeholder counts; the table has all that is needed here)
    se0 <- SummarizedExperiment(list(counts=matrix(30L, nrow=n, ncol=6,
                                      dimnames=list(row.names(ann), NULL))))
    se <- exp(rnorm(n, log(0.12), 0.3))
    sim <- function(eff){
      data.frame(logFC=eff+rnorm(n, 0, 0.25)+rnorm(n, 0, se), SE=se, 
                 AveExpr=rnorm(n, 6, 1), row.names=row.names(ann))
    }
    tb <- sim(0.6*A + 0.4*C)
    # B is falsely flagged by the single-motif tests...
    r <- suppressMessages(peakLevelMotifTests(se0, tb, ann, minCount=10))
    expect_lt(r$adj.P.Val[r$motif=="B"], 1e-10)
    expect_gt(r$adj.P.Val[r$motif=="D"], 0.05)
    # ... but not by the stepwise procedure
    st <- suppressMessages(stepwisePeakLevelMotifTests(se0, tb, ann, 
                                                       minCount=10))
    expect_named(st, c("selected","remaining"))
    expect_equal(st$selected$motif, c("A","C"))
    expect_equal(st$selected$step, 1:2)
    expect_true(all(st$selected$adj.P.Val<0.05))
    expect_equal(st$selected$joint.difference, c(0.6, 0.4), tolerance=0.1)
    expect_true(all(c("joint.se","joint.P.Value") %in% colnames(st$selected)))
    expect_setequal(st$remaining$motif, c("B","D"))
    expect_true(all(st$remaining$adj.P.Val>0.05))
    # the first selected motif is the one of the single-motif tests
    expect_equal(st$selected$difference[1], r$difference[r$motif=="A"])
    # maxSteps
    st1 <- suppressMessages(stepwisePeakLevelMotifTests(se0, tb, ann, 
                                                        minCount=10, 
                                                        maxSteps=1))
    expect_equal(st1$selected$motif, "A")
    expect_true(all(c("B","C","D") %in% st1$remaining$motif))
    # nothing selected when there's no effect
    st0 <- suppressMessages(stepwisePeakLevelMotifTests(se0, sim(0*A), ann, 
                                                        minCount=10))
    expect_equal(nrow(st0$selected), 0L)
    expect_true(all(c("contrast","step","motif","joint.difference") %in% 
                      colnames(st0$selected)))
    expect_equal(nrow(st0$remaining), 4L)
    # a fully collinear motif cannot be tested after its twin is selected
    ann2 <- cbind(ann, A2=A)
    st2 <- suppressMessages(stepwisePeakLevelMotifTests(se0, tb, ann2, 
                                                        minCount=10))
    expect_true(st2$selected$motif[1] %in% c("A","A2"))
    expect_equal(st2$selected$motif[2], "C")
    expect_true(all(is.na(st2$remaining$difference[st2$remaining$motif %in% 
                                                     c("A","A2")])))
  })
  
  test_that("stepwisePeakLevelMotifTests works on the example data", {
    st <- stepwisePeakLevelMotifTests(peakSE, d$dea, motifMatches, 
                                      design=d$design, minCount=10)
    expect_true(all(st$selected$motif[1:2] %in% c("GCR","MCR")))
    expect_true(all(st$selected$step==seq_len(nrow(st$selected))))
    r <- peakLevelMotifTests(peakSE, d$dea, motifMatches, design=d$design, 
                             minCount=10)
    expect_equal(st$selected$difference[1], 
                 r$difference[r$motif==st$selected$motif[1]])
    expect_false(any(st$selected$motif %in% st$remaining$motif))
    expect_error(stepwisePeakLevelMotifTests(peakSE, d$dea, motifMatches, 
                                             alpha=2))
  })
  
  test_that("adjustWidth adds a width covariate", {
    r0 <- peakLevelMotifTests(peakSE, d$dea, motifMatches, design=d$design, 
                              minCount=10)
    se2 <- peakSE
    if(methods::is(se2, "RangedSummarizedExperiment")) 
      se2 <- SummarizedExperiment::SummarizedExperiment(
        SummarizedExperiment::assays(se2), 
        colData=SummarizedExperiment::colData(se2),
        rowData=SummarizedExperiment::rowData(se2))
    # no width information: error
    expect_error(peakLevelMotifTests(se2, d$dea, motifMatches, design=d$design,
                                     minCount=10, adjustWidth=TRUE), "width")
    # widths in rowData
    set.seed(1)
    SummarizedExperiment::rowData(se2)$width <- 
      round(exp(rnorm(nrow(se2), 6, 0.4)))
    r1 <- peakLevelMotifTests(se2, d$dea, motifMatches, design=d$design, 
                              minCount=10, adjustWidth=TRUE)
    expect_setequal(r1$motif, r0$motif)
    expect_false(isTRUE(all.equal(r1$difference, 
                                  r0$difference[match(r1$motif, r0$motif)])))
    # constant widths: nothing to adjust for
    SummarizedExperiment::rowData(se2)$width <- 500L
    expect_message(r2 <- peakLevelMotifTests(se2, d$dea, motifMatches, 
                         design=d$design, minCount=10, adjustWidth=TRUE),
                   "constant")
    expect_equal(r2$difference, r0$difference)
    expect_no_error(suppressMessages(stepwisePeakLevelMotifTests(se2, d$dea, 
                     motifMatches, design=d$design, minCount=10, 
                     adjustWidth=TRUE)))
  })
})
