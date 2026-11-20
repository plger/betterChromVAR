#' getInteractionsDeviations
#'
#' Get \code{\link[betterChromVAR]{betterChromVAR}} deviations for the sites
#' of a bait motif, split according to whether they also harbor each of the
#' other motifs.
#'
#' @param se An object inheriting `SummarizedExperiment`, containing a `counts`
#'   assay.
#' @param annotation A (sparse) matrix or object inheriting
#'   `SummarizedExperiment` containing motif annotations for each row of `se`.
#'   Weighted annotations are binarized (any positive value is a match).
#' @param bait The name of the bait motif (corresponding to a column of
#'   `annotation`), for which to search for interacting motifs.
#' @param minCount The minimum number of sites in each of the sets
#'   (see Details) for a motif to be considered.
#' @param ... Passed to `betterChromVAR`. This includes the option to use the
#'   `ridge` argument.
#'
#' @details
#' For each other motif (the 'partner'), three disjoint sets of sites are
#' defined: the sites harboring both the bait and the partner (`both`), only the
#' bait (`baitOnly`), or only the partner (`partnerOnly`). Deviations are
#' computed for each, as well as for the bait motif as a whole. Because the sets
#' are disjoint, their deviations can be compared without one containing the
#' other. In addition to `deviations` and `z`, the returned object has a
#' `sd` assay, which is the standard deviation of the deviations under the null
#' (i.e. `deviations/z`), and is used by \code{\link{discoverMotifInteractions}}
#' to weight the observations.
#'
#' @returns A `SummarizedExperiment` of deviations for the bait and for the sets
#'   of sites defined by the intersection of the bait with each other motif. The
#'   `rowData` indicates the `partner` motif and the type of `set` of each row,
#'   as well as overlap statistics between the bait and the partner.
#' @export
#' @examples
#' attach(getDummyData())
#' dev <- getInteractionsDeviations(counts, motifMatches, bait="motif1",
#'                                  minCount=5)
#' dev
getInteractionsDeviations <- function(se, annotation, bait, minCount=20, ...){
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(nrow(se)==nrow(annotation))
  stopifnot(length(bait)==1 && bait %in% colnames(annotation))
  if(inherits(annotation, "SummarizedExperiment")){
    annotation <- assay(annotation)
  }
  annotation <- drop0(as(annotation>0, "CsparseMatrix"))
  ab <- as.vector(annotation[,bait])
  nBait <- sum(ab)
  if(nBait<minCount) stop("The bait motif has fewer than `minCount` matches.")
  partners <- setdiff(colnames(annotation), bait)
  P <- annotation[,partners,drop=FALSE]*1

  # three disjoint sets of sites for each partner:
  both <- Diagonal(x=as.numeric(ab)) %*% P
  baitOnly <- sparseMatrix(i=rep(which(ab), length(partners)), x=1,
                           j=rep(seq_along(partners), each=nBait),
                           dims=dim(P)) - both
  partnerOnly <- P - both
  sets <- list(both=both, baitOnly=baitOnly, partnerOnly=partnerOnly)
  sets <- lapply(sets, \(x) drop0(as(drop0(x)>0, "CsparseMatrix")))
  sizes <- sapply(sets, MatrixGenerics::colSums)
  keep <- which(rowSums(sizes>=minCount)==3L)
  if(length(keep)==0) stop("No interaction passing filtering.")
  partners <- partners[keep]
  sets <- lapply(sets, \(x) x[,keep,drop=FALSE])

  moi <- cbind(Matrix(ab, ncol=1, sparse=TRUE), sets$both, sets$baitOnly,
               sets$partnerOnly)
  colnames(moi) <- c(bait, paste0(partners, ".", rep(names(sets),
                                                     each=length(partners))))
  dev <- betterChromVAR(se, moi, ...)

  # null SD of the deviations (dev/z is exactly the ratio of the two scales)
  s <- abs(assay(dev, "deviations")/assay(dev, "z"))
  bad <- !is.finite(s) | s==0
  if(any(bad)){
    rmed <- apply(replace(s, bad, NA_real_), 1, stats::median, na.rm=TRUE)
    s[bad] <- rmed[row(s)[bad]]
  }
  assay(dev, "sd") <- s

  # overlap statistics between the bait and each partner
  ov <- sizes[keep,"both"]
  tot <- MatrixGenerics::colSums(P[,keep,drop=FALSE])
  expected <- nBait*tot/nrow(annotation)
  jac <- as.numeric(.colOverlaps(P[,keep,drop=FALSE], ab, metric="jaccard"))
  partner <- c(bait, rep(partners, 3))
  m <- match(partner, partners)
  rowData(dev)$partner <- partner
  rowData(dev)$set <- factor(c("bait", rep(names(sets), each=length(partners))),
                             c("bait", names(sets)))
  # (unnamed, as the NA index of the bait would otherwise give NA names, which
  # break the display of the rowData)
  rowData(dev)$overlap.percentOfBait <- unname(round(100*(ov/nBait),2)[m])
  rowData(dev)$overlap.jaccard <- unname(jac[m])
  rowData(dev)$overlap.log2Enr <- unname(log2(ov/expected)[m])
  rowData(dev)$isBait <- row.names(dev)==bait
  dev
}

#' .colOverlaps
#'
#' Computes pairwise overlap metric between columns (taken from epiwraps)
#'
#' @param x A logical matrix.
#' @param y An optional nother logical matrix or vector. If NULL (default),
#'   overlaps are computed between columns of `x`.
#' @param metric Either 'overlap', 'jaccard', or 'overlapCoef'.
#'
#' @returns A matrix of pairwise metric values.
#' @noRd
.colOverlaps <- function(x, y=NULL, metric=c("overlap","jaccard","overlapCoef")){
  metric <- match.arg(metric)
  x <- as(x, "lMatrix")
  ys <- xs <- colSums(x)
  if(is.null(y)){
    y <- x
  }else{
    if(is.vector(y)) y <- as.matrix(y)
    y <- as(y, "lMatrix")
    ys <- colSums(y)
  }
  intersections <- t(x) %*% as.matrix(y)
  if(metric=="overlap"){
    out <- intersections
  }else if(metric=="overlapCoef"){
    minSize <- outer(xs, ys, "min")
    out <- intersections/minSize
  }else{
    unions <- outer(xs, ys, "+") - intersections
    out <- intersections/unions
  }
  if(is(out, "dgeMatrix")) out <- as.matrix(out)
  out
}

# design matrix from a group variable and covariates; the group coefficients
# are named after the group levels, relative to the first
.interactionDesign <- function(cd, group, covar=c()){
  stopifnot(group %in% colnames(cd))
  stopifnot(all(covar %in% colnames(cd)))
  cd <- as.data.frame(cd[,c(group,covar),drop=FALSE])
  grp <- factor(cd[[group]])
  if(nlevels(grp)<2) stop("`group` should have at least two levels.")
  cd <- cbind(data.frame(.group=grp), cd[,covar,drop=FALSE])
  f <- paste(c(".group", sprintf("`%s`", covar)), collapse="+")
  mm <- model.matrix(as.formula(paste0("~", f)), data=cd)
  gl <- levels(grp)[-1]
  cn <- paste0(gl, "-", levels(grp)[1])
  coefs <- match(paste0(".group", gl), colnames(mm))
  list(design=mm, coefs=setNames(coefs, cn))
}

# overlap statistics columns of a (partner-level) data.frame
.overlapCols <- function(rd, rows){
  cols <- grep("^overlap", colnames(rd), value=TRUE)
  d <- as.data.frame(rd[rows, cols, drop=FALSE])
  d$overlap.N <- rd$N[rows]
  d[,c("overlap.N", setdiff(cols, "overlap.N")),drop=FALSE]
}

#' discoverMotifInteractions
#'
#' Discover motifs interacting with a bait motif, based on the betterChromVAR
#' package.
#'
#' @param dev A deviations SummarizedExperiment produced by
#'   \code{\link{getInteractionsDeviations}}
#' @param group The column of `colData(dev)` containing the experimental groups
#'   of interest. To column will be coerced to a factor, and all comparisons
#'   will be made relative to the first level of the factor.
#' @param covar Optional columns of `colData(dev)` containing covariates to
#'   correct for.
#' @param global Logical; whether to test the contrasts globally, rather than
#'   individually. Has no effect when there is a single contrast
#'   (i.e. two-group comparison).
#' @param weights The weights to use, either 'precision' (default; weights each
#'   sample by the inverse of the variance of the deviations under the null,
#'   which accounts for differences in set size and sequencing depth), or
#'   'none' (no weights; not recommended).
#' @author Pierre-Luc Germain
#'
#' @details
#' For each partner motif, the deviations of the three disjoint sets of sites
#' (bait and partner, bait only, partner only; see
#' \code{\link{getInteractionsDeviations}}) are combined in each sample in an
#' interaction score,
#' `both - baitOnly - partnerOnly`,
#' i.e. the extent to which the deviation at sites with both motifs differs from
#' what would be expected if the bait and the partner acted in an additive 
#' fashion on the accessibility deviations. This score is then regressed on 
#' `group` (and `covar`) with limma.
#'
#' Note that the p-values are meant to be used for ranking rather than 
#' calibrated error rates. Moreover, since the classes of sites are aggregated, 
#' the p-values indicate reflect the reproducibility across replicates, rather 
#' than the generalization across sites. For a test that accounts for the 
#' variability across sites and produces calibrated p-values, see 
#' \code{\link{peakLevelMotifInteractions}}.
#'
#' @returns A data.frame with, for each motif and contrast, the `difference` in
#'   interaction score between the groups, its moderated statistics, and the
#'   differences between groups in each of the sets (`diff.both`,
#'   `diff.baitOnly`, `diff.partnerOnly`.
#' @importFrom limma lmFit eBayes topTable
#' @importFrom stats model.matrix as.formula
#' @export
#' @examples
#' attach(getDummyData())
#' dev <- getInteractionsDeviations(counts, motifMatches, bait="motif1",
#'                                  minCount=5)
#' dev$group <- rep(LETTERS[1:2], each=5)
#' discoverMotifInteractions(dev, "group")
discoverMotifInteractions <- function(dev, group, covar=c(), global=FALSE,
                                      weights=c("precision","none")){
  stopifnot(inherits(dev, "SummarizedExperiment"))
  rd <- rowData(dev)
  if(!all(c("partner","set") %in% colnames(rd)) || !("sd" %in% assayNames(dev)))
    stop("`dev` should be the output of `getInteractionsDeviations`.")
  stopifnot(sum(rd$isBait)==1)
  weights <- match.arg(weights)

  dm <- .interactionDesign(colData(dev), group, covar)
  partners <- unique(rd$partner[!rd$isBait])
  ix <- lapply(setNames(nm=c("both","baitOnly","partnerOnly")), \(st){
    i <- which(rd$set==st)
    i[match(partners, rd$partner[i])]
  })
  if(anyNA(unlist(ix))) stop("Some partner motifs lack one of the sets.")

  d <- assay(dev, "deviations")
  s2 <- assay(dev, "sd")^2
  # the interaction score and its (null) variance
  inter <- d[ix$both,,drop=FALSE] - d[ix$baitOnly,,drop=FALSE] -
    d[ix$partnerOnly,,drop=FALSE]
  v <- s2[ix$both,,drop=FALSE] + s2[ix$baitOnly,,drop=FALSE] +
    s2[ix$partnerOnly,,drop=FALSE]
  row.names(inter) <- row.names(v) <- partners
  w <- switch(weights,
              precision=(1/v)/stats::median(1/v),
              NULL)

  trend <- nrow(inter)>=20
  if(trend){
    N <- rd$N[ix$both]
    trend <- log(N)
  }
  fit <- eBayes(lmFit(inter, dm$design, weights=w), trend=trend)

  # response of each set to the group
  wi <- NULL
  setDiff <- sapply(ix, \(i){
    if(weights!="none") wi <- 1/s2[i,,drop=FALSE]
    lmFit(d[i,,drop=FALSE], dm$design, weights=wi)$coefficients[,dm$coefs,
                                                                drop=FALSE]
  }, simplify=FALSE)

  getRes <- function(co){
    res <- topTable(fit, coef=co, number=Inf, sort.by="none")
    colnames(res)[1] <- "difference"
    res$B <- res$AveExpr <- NULL
    res
  }
  if(!isTRUE(global) || length(dm$coefs)==1){
    res <- do.call(rbind, lapply(seq_along(dm$coefs), \(i){
      res <- getRes(dm$coefs[i])
      for(st in names(ix))
        res[[paste0("diff.",st)]] <- setDiff[[st]][,i]
      cbind(contrast=names(dm$coefs)[i], motif=partners,
            .overlapCols(rd, ix$both), res)
    }))
    res <- res[order(res$contrast, res$P.Value),]
  }else{
    res <- topTable(fit, coef=dm$coefs, number=Inf, sort.by="none")
    res$AveExpr <- NULL
    colnames(res)[seq_along(dm$coefs)] <- paste0("difference.",
                                                 names(dm$coefs))
    res <- cbind(motif=partners, .overlapCols(rd, ix$both), res)
    res <- res[order(res$P.Value),]
  }
  row.names(res) <- NULL
  res
}

# Estimates the mean correlation between peaks of the residuals of a linear
# model, as in limma::camera, and returns the variance inflation factor.
.excessCorrelation <- function(E, design, idx, rhoAll=0){
  if(length(idx)<2) return(0)
  rho <- limma::interGeneCorrelation(E[idx,,drop=FALSE], design)$correlation
  # the component shared by all peaks (e.g. sample-level noise) cancels out of
  # the contrasts between sets, so only the excess is relevant
  min(0.99, max(0, rho-rhoAll))
}

# Weighted regression of logFC across peaks with random-effect heterogeneity
# (per cell) and a correlation-based variance inflation (per cell):
#   y_p = x_p' beta + e_p, var(e_p)=v_p+tau_c^2, cor(e_p,e_q)=rho_c (same cell c)
.peakMetaRegression <- function(y, v, X, cell, rho, nIter=2L){
  ncell <- length(rho)
  tau2 <- rep(0, ncell)
  for(i in seq_len(nIter+1L)){
    w <- 1/(v+tau2[cell])
    XtWX <- crossprod(X*w, X)
    beta <- solve(XtWX, crossprod(X*w, y))
    if(i>nIter) break
    r <- as.numeric(y-X%*%beta)
    tau2 <- vapply(seq_len(ncell), \(k){
      j <- which(cell==k)
      if(length(j)<2) return(0)
      max(0, mean(r[j]^2-v[j])*length(j)/(length(j)-1))
    }, numeric(1))
  }
  # sandwich variance, with the correlation structure within each cell
  meat <- XtWX
  for(k in which(rho>0)){
    j <- which(cell==k)
    a <- X[j,,drop=FALSE]*(w[j]*sqrt(v[j]))
    meat <- meat + rho[k]*(tcrossprod(colSums(a)) - crossprod(a))
  }
  bread <- solve(XtWX)
  list(beta=as.numeric(beta), vcov=bread %*% meat %*% bread, tau2=tau2)
}

# Per-peak differential accessibility (limma-voom) and the elements needed for 
# the peak-level regressions: log-fold-changes and variances for each contrast,
# the adjustment covariates (splines of accessibility and bias), and the 
# correlation between peaks.
.peakLevelSetup <- function(se, dm, splineDF=3, vif=TRUE, useAssay="counts",
                            lib.size=NULL, adjustWidth=FALSE){
  if(is.null(lib.size) && !is.null(se$depth)) lib.size <- se$depth
  v <- voom(as.matrix(assay(se, useAssay)), dm$design, lib.size=lib.size, 
            plot=FALSE)
  fit <- eBayes(lmFit(v, dm$design))
  lfc <- fit$coefficients[,dm$coefs,drop=FALSE]
  vr <- (fit$stdev.unscaled[,dm$coefs,drop=FALSE])^2 * fit$s2.post
  
  covs <- list()
  if(splineDF>0){
    covs$acc <- ns(fit$Amean, df=splineDF)
    bias <- rowData(se)$bias
    if(!is.null(bias) && length(unique(bias))>splineDF+1L)
      covs$bias <- ns(bias, df=splineDF)
  }
  if(adjustWidth) covs$width <- .widthCovariate(se, rep(TRUE, nrow(se)), splineDF)
  covs <- do.call(cbind, covs)
  
  rhoAll <- 0
  if(vif) rhoAll <- max(0, interGeneCorrelation(v$E, dm$design)$correlation)
  list(v=v, lfc=lfc, vr=vr, covs=covs, rhoAll=rhoAll, design=dm$design)
}

#' peakLevelMotifInteractions
#'
#' Tests the interaction between a bait motif and other motifs on the
#' differential accessibility of peaks, using the peaks (rather than the
#' samples) as units of inference.
#'
#' @param se An object inheriting `SummarizedExperiment`, containing a `counts`
#'   assay (or `useAssay`), and optionally a `bias` column in its `rowData`
#'   (see \code{\link[betterChromVAR]{getBackgroundBins}}).
#' @param annotation A (sparse) matrix or object inheriting
#'   `SummarizedExperiment` containing motif annotations for each row of `se`.
#' @param bait The name of the bait motif (corresponding to a column of
#'   `annotation`), for which to search for interacting motifs.
#' @param group The column of `colData(se)` containing the experimental groups
#'   of interest. The column will be coerced to a factor, and all comparisons
#'   will be made relative to the first level of the factor.
#' @param covar Optional columns of `colData(se)` containing covariates to
#'   correct for.
#' @param minCount The minimum number of peaks in each of the sets
#'   (bait and partner, bait only, partner only) for a motif to be considered.
#' @param splineDF The degrees of freedom of the splines used to adjust for
#'   accessibility and bias (set to 0 to disable).
#' @param adjustWidth Logical; whether to additionally adjust for the (log-)
#'   width of the peaks (default FALSE), taken from `rowRanges(se)` (or from a 
#'   `width` column of `rowData(se)`). See \code{\link{peakLevelMotifTests}}.
#' @param vif Logical; whether to inflate the variance to account for the
#'   correlation between the peaks of a set (default TRUE).
#' @param useAssay The assay containing the counts.
#' @param lib.size Optional library sizes (default is the column sums of the
#'   counts), passed to \code{\link[limma]{voom}}.
#'
#' @details
#' The differential accessibility between groups (log2-fold-change and its
#' standard error) is first estimated for each peak using limma-voom, with
#' `group` and `covar` in the model. For each partner motif and contrast, the
#' log-fold-changes of all peaks are then modeled as a function of the presence
#' of the bait motif, of the partner motif, and of their interaction, together
#' with smooth functions (natural splines) of the peak's mean accessibility and
#' bias. The interaction coefficient is the extent to which the response at
#' sites with both motifs differs from the sum of the effects of each motif
#' alone (on the log scale).
#'
#' The regression is weighted by the inverse of the variance of the peak's
#' log-fold-change, plus a between-peak variance (random effect) estimated
#' separately for each of the four types of peaks, so that the heterogeneity of
#' the response across sites is part of the uncertainty. In addition, as in
#' \code{\link[limma]{camera}}, the variance of the coefficients is inflated to
#' account for the average correlation between the residuals of the peaks
#' of a same set (estimated across samples; only the part exceeding the
#' correlation among all peaks is used, as the latter cancels out of the
#' contrasts). P-values use a t-distribution with the number of peaks
#' harboring both motifs minus one as degrees of freedom.
#'
#' The model is fitted separately for each partner, and treats the peaks
#' as independent (beyond the correlation correction), so that very similar
#' motifs matching the same sites cannot be distinguished from a true
#' interaction.
#'
#' This approach is related to motif activity response analysis (MARA, and its 
#' web implementation ISMARA; see the references), which models the signal of 
#' regulatory regions as a function of the motifs they harbor. It differs in 
#' that the bait and each partner are modeled separately from the other motifs 
#' (MARA fits all motifs jointly), that it includes an interaction term, and 
#' that it works on the differential accessibility between groups, with a 
#' random effect for the heterogeneity across sites and a correction for the 
#' correlation between peaks.
#'
#' @references
#'   The FANTOM Consortium and the RIKEN Omics Science Center (2009) The
#'   transcriptional network that controls growth arrest and differentiation
#'   in a human myeloid leukemia cell line, Nature Genetics 41:553-562,
#'   doi: 10.1038/ng.375.
#'
#'   Balwierz P.J., Pachkov M., Arnold P., Gruber A.J., Zavolan M.,
#'   van Nimwegen E. (2014) ISMARA: automated modeling of genomic signals as a
#'   democracy of regulatory motifs, Genome Research 24:869-884,
#'   doi: 10.1101/gr.169508.113
#'
#' @returns A data.frame with, for each motif and contrast, the `difference`
#'   (interaction coefficient, in log2-fold-change), its standard error (`se`),
#'   t-value, p-value and adjusted p-value, the effects of the bait alone
#'   (`bait.diff`) and of the partner alone (`partner.diff`), and the estimated
#'   heterogeneity (`tau.both`, standard deviation of the true response across
#'   sites with both motifs) and excess correlation (`rho.both`) for the
#'   peaks with both motifs.
#' @importFrom limma voom lmFit eBayes interGeneCorrelation
#' @importFrom splines ns
#' @importFrom stats p.adjust pt model.matrix ave
#' @export
#' @examples
#' attach(getDummyData())
#' counts$group <- rep(LETTERS[1:2], each=5)
#' peakLevelMotifInteractions(counts, motifMatches, bait="motif1",
#'                            group="group", minCount=5)
peakLevelMotifInteractions <- function(se, annotation, bait, group, covar=c(),
                                       minCount=20, splineDF=3, adjustWidth=FALSE, 
                                       vif=TRUE, useAssay="counts", lib.size=NULL){
  stopifnot(inherits(se, "SummarizedExperiment"))
  stopifnot(nrow(se)==nrow(annotation))
  stopifnot(length(bait)==1 && bait %in% colnames(annotation))
  if(inherits(annotation, "SummarizedExperiment")){
    annotation <- assay(annotation)
  }
  annotation <- drop0(as(annotation>0, "CsparseMatrix"))
  dm <- .interactionDesign(colData(se), group, covar)

  pl <- .peakLevelSetup(se, dm, splineDF, vif, useAssay, lib.size, adjustWidth)
  v <- pl$v
  lfc <- pl$lfc
  vr <- pl$vr
  covs <- pl$covs
  rhoAll <- pl$rhoAll

  ab <- as.vector(annotation[,bait])
  if(sum(ab)<minCount) stop("The bait motif has fewer than `minCount` matches.")
  partners <- setdiff(colnames(annotation), bait)
  P <- annotation[,partners,drop=FALSE]
  tot <- MatrixGenerics::colSums(P)
  nBoth <- as.numeric(MatrixGenerics::colSums(P & ab))
  keep <- which(nBoth>=minCount & (sum(ab)-nBoth)>=minCount &
                  (tot-nBoth)>=minCount)
  if(length(keep)==0) stop("No interaction passing filtering.")

  res <- lapply(partners[keep], \(m){
    pp <- as.vector(annotation[,m])
    cell <- 1L + ab + 2L*pp  # 1: neither, 2: bait only, 3: partner only, 4: both
    X <- cbind(1, ab, pp, ab*pp, covs)
    rho <- rep(0, 4)
    if(vif) rho <- vapply(1:4, \(k) .excessCorrelation(v$E, dm$design,
                                                        which(cell==k), rhoAll),
                          numeric(1))
    n <- sum(cell==4L)
    out <- lapply(seq_along(dm$coefs), \(i){
      r <- .peakMetaRegression(lfc[,i], vr[,i], X, cell, rho)
      se <- sqrt(diag(r$vcov))
      data.frame(contrast=names(dm$coefs)[i], motif=m,
                 overlap.N=n, N.baitOnly=sum(cell==2L),
                 N.partnerOnly=sum(cell==3L),
                 difference=r$beta[4], se=se[4], t=r$beta[4]/se[4],
                 bait.diff=r$beta[2], partner.diff=r$beta[3],
                 tau.both=sqrt(r$tau2[4]), rho.both=rho[4])
    })
    do.call(rbind, out)
  })
  res <- do.call(rbind, res)
  res$P.Value <- 2*pt(-abs(res$t), df=pmax(1, res$overlap.N-1))
  res$adj.P.Val <- ave(res$P.Value, res$contrast, FUN=\(p) p.adjust(p, "fdr"))

  # overlap statistics
  ov <- nBoth[keep]
  expected <- sum(ab)*tot[keep]/nrow(annotation)
  jac <- as.numeric(.colOverlaps(P[,keep,drop=FALSE], ab, metric="jaccard"))
  m <- match(res$motif, partners[keep])
  res$overlap.percentOfBait <- round(100*(ov/sum(ab)),2)[m]
  res$overlap.jaccard <- jac[m]
  res$overlap.log2Enr <- log2(ov/expected)[m]
  res <- res[order(res$contrast, res$P.Value),]
  row.names(res) <- NULL
  res
}

#' exploreMotifInteraction
#'
#' We compare deviations in sites that harbor a single of the two motifs, or
#' the pair of motifs at given distances from each other.
#'
#' @param counts An object inheriting from `RangedSummarizedExperiment`, with
#'   a `bias` column in its `rowData`.
#' @param motifs A named pair of motifs, in a format recognized by
#'   \code{\link[motifmatchr]{matchMotifs}}.
#' @param genome A BSgenome or FaFile object.
#' @param bg Optional pre-computed background (a `bcvBackground` object, as
#'   produced by \code{\link[betterChromVAR]{computeBackgrounds}}). If given,
#'   the deviations are computed with
#'   \code{\link[betterChromVAR]{computeDeviationsAnalytic}} from this
#'   background; otherwise (default) they are computed with
#'   \code{\link[betterChromVAR]{betterChromVAR}}, which estimates it.
#' @param nDistQuantiles The number of distance quantiles between the motif
#'   matches (up to `maxDist`)
#' @param maxDist The maximum distance between matches.
#' @param minDist The minimum distance between the centers of the matches. This
#'   excludes overlapping matches, which for similar motifs are not independent
#'   observations. To identify competitors (i.e. motifs that bind overlapping
#'   sites), this should be disabled by setting it to 0.
#' @param restrictToRegions Logical; whether both motifs should be within the
#'   peaks (default FALSE, i.e. one of the matches can be up to `maxDist`
#'   outside of the peak, see \code{\link{motifCoOccurence}}).
#' @param ... Any further argument passed to
#'   \code{\link[betterChromVAR]{betterChromVAR}} (or to
#'   \code{\link[betterChromVAR]{computeDeviationsAnalytic}} if `bg` is given).
#' @author Pierre-Luc Germain
#'
#' @returns A SummarizedExperiment object of deviations for different types of
#'   interactions between the motifs.
#'
#' @details
#' Motif pairs are split according to the distance between their centers
#' (between `minDist` and `maxDist`, default 5 to 500) into `nDistQuantiles`
#' (default 5) similarly-populated quantiles. Deviations are reported for these
#' different bins, as well as for sites harboring only a single motif (i.e.
#' having a match for one of the motifs, and no match of the other within the
#' peak or within `maxDist`). Peaks with both motifs but at a distance that
#' falls outside of the range of distances are in neither of these sets. If the
#' distances are heavily tied, fewer than `nDistQuantiles` bins are returned.
#' @export
exploreMotifInteraction <- function(counts, motifs, genome, bg=NULL,
                                    nDistQuantiles=5L, maxDist=500L,
                                    minDist=5L, restrictToRegions=FALSE, ...){
  stopifnot(length(motifs)==2)
  stopifnot(inherits(counts, "RangedSummarizedExperiment"))
  if(!requireNamespace("motifmatchr", quietly=TRUE))
    stop("The `motifmatchr` package is required for this function.")
  moi <- motifmatchr::matchMotifs(motifs, counts, genome=genome)
  moi2 <- motifCoOccurence(motifs, list(names(motifs)), rowRanges(counts),
                           genome = genome, nDistQuantiles=nDistQuantiles,
                           minDist=minDist, maxDist=maxDist,
                           restrictToRegions=restrictToRegions)
  breaks <- attr(moi2, "breaks")[[1]][-1]
  fmt <- function(x) format(round(x,1), trim=TRUE)
  bins <- c(paste0("d<=", fmt(breaks[1])),
            paste0(fmt(breaks[-length(breaks)]), "<d<=", fmt(breaks[-1])))
  # (empty bins can occur with heavily-tied distances)
  moi2 <- moi2[seq_along(bins)]
  moi2flat <- Reduce(cbind, moi2)
  colnames(moi2flat) <- paste(colnames(moi2[[1]]), bins, sep="\n")

  # sites with only one of the motifs, and no partner in the vicinity
  a <- assay(moi)>0
  cooc <- Reduce("|", lapply(moi2, \(x) Matrix::rowSums(x)>0))
  a1 <- as.vector(a[,1])
  a2 <- as.vector(a[,2])
  moiSinglet <- Matrix(cbind(a1 & !a2 & !cooc, a2 & !a1 & !cooc), sparse=TRUE)
  colnames(moiSinglet) <- paste0(colnames(a),"+,", rev(colnames(a)),"-")
  moi <- cbind(a, moiSinglet, moi2flat)
  if(is.null(bg)){
    dev <- betterChromVAR(counts, moi, ...)
  }else{
    dev <- computeDeviationsAnalytic(counts, bg, moi, ...)
  }
  normalizeDevsForSize(dev)
}

