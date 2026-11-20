#' peakLevelMotifTests
#' 
#' Differential motif accessibility between groups based on the peaks' logFC
#' between conditions. This is a peak-level alternative to the analysis of 
#' motif deviations (see \code{\link{fitDeviations}}).
#'
#' @param se An object inheriting `SummarizedExperiment` containing the peak
#'   counts (in the `counts` assay, or `useAssay`), from which the differential
#'   analysis `dea` was obtained. It defines the peaks (rows), and provides the
#'   sample-level data used to estimate the correlation between peaks, as well
#'   as the `bias` column of `rowData` (see \code{\link{addGCBias}}),
#'   if present (which is very much recommended).
#' @param dea The table of the results of a differential accessibility analysis
#'   of the peaks of `se` (or a named list of such tables, one for each 
#'   contrast). Should accept major differential analysis packages (see Details
#'   for the format).
#' @param annotation A (sparse) matrix or object inheriting 
#'   `SummarizedExperiment` containing motif annotations for each row of `se`.
#'   Weighted annotations are binarized (any positive value is a match).
#' @param design Optional design matrix (with a row for each column of `se`) 
#'   of the model used for the differential analysis. It is used to estimate the
#'   correlation between the peaks (see `vif`) from the residuals of the model.
#' @param minCount The minimum number of peaks with the motif (and without it) 
#'   for a motif to be tested.
#' @param splineDF The degrees of freedom of the splines used to adjust for 
#'   accessibility and bias (set to 0 to disable).
#' @param adjustWidth Logical; whether to additionally adjust for the (log-)
#'   width of the peaks (default FALSE), taken from `rowRanges(se)` (or from a 
#'   `width` column of `rowData(se)`). This is useful when peaks have variable 
#'   widths, since wider peaks are both more likely to contain a motif and have
#'   higher counts. Peaks of constant width are not affected. Alternatively, 
#'   the peaks can be resized to a fixed width before counting. 
#' @param vif Controls the inflation of the variance to account for the 
#'   correlation between the peaks of a set. Either `TRUE` (default; the
#'   correlation is estimated from the data of `se`, which requires `design`; 
#'   if `design` is missing the inflation is disabled), `FALSE` (no 
#'   inflation; not recommended for large motif sets), or a number between 0 
#'   and 1, which is the excess correlation assumed between the peaks of each 
#'   motif (see Details).
#' @param useAssay The assay of `se` containing the counts.
#' @param lib.size Optional library sizes (default is the column sums of the 
#'   counts), used to compute the log-counts-per-million from which the 
#'   correlation between peaks is estimated.
#'
#' @details
#' The differential accessibility between groups (log2-fold-change and its 
#' standard error) is first estimated for each peak, using any method, and given
#' as a table (`dea`) with the peaks of `se` as rows (matched to the rows of 
#' `se` through the row names if both have some, and otherwise in the same 
#' order). The table should contain the log-fold-change (column `logFC` or 
#' `log2FoldChange`), its standard error (column `SE`, `se` or `lfcSE`) or a 
#' statistic from which it can be derived (column `t`, `stat`, `F` or `LR`), 
#' and optionally the average accessibility (column `AveExpr`, `logCPM`, 
#' `Amean`, or `baseMean`; if missing, the mean log-counts-per-million of `se`
#' is used), as given for instance by limma's `topTable` (with 
#' `sort.by="none"`), edgeR's `topTags` or DESeq2's `results`. Peaks with 
#' missing values are ignored. A `bias` column can be included to adjust for it
#' (otherwise the `bias` column of the `rowData` of `se` is used, if present). 
#' The standard error should be the one of the log-fold-change itself (and not 
#' of another quantity).
#' 
#' For each motif and contrast, the log-fold-changes of all peaks are then 
#' modeled as a function of the presence of the motif, together with smooth 
#' functions (natural splines) of the peak's mean accessibility and bias (and
#' optionally width, see `adjustWidth`). The coefficient of the motif is thus 
#' the extent to which peaks with the motif change more (or less) than 
#' comparable peaks without it, which is the peak-level analog of a motif 
#' deviation (that is, a competitive test).
#' 
#' The regression is weighted by the inverse of the variance of the peak's 
#' log-fold-change, plus a between-peak variance (random effect) estimated 
#' separately for peaks with and without the motif, so that the heterogeneity of
#' the response across sites is part of the uncertainty. In addition, as in 
#' \code{\link[limma]{camera}} (Wu and Smyth 2012), the variance of the 
#' coefficients can be inflated to account for the average correlation between 
#' the residuals of the peaks of a same set (estimated across samples from the 
#' residuals of the model given by `design`; only the part exceeding the 
#' correlation among all peaks is used, as the latter cancels out of the 
#' contrast). If `design` is not given, an excess correlation can alternatively
#' be assumed through `vif`. (Since the inflation grows with the number of peaks 
#' of the motif, omitting it can make motifs with many peaks appear markedly 
#' more significant than they should.)
#' 
#' Compared to the analysis of deviations, this test asks whether the sites of 
#' a motif change more than other sites (generalizing across sites), rather
#' than whether the set changes reproducibly across samples, and does not 
#' rely on sampling a background. Its power depends on the number of sites of 
#' the motif and on how consistently they respond. Each motif is tested 
#' separately, so motifs that share sites are not distinguished (see 
#' \code{\link{stepwisePeakLevelMotifTests}} for a way to address this, or the 
#' `ridge` argument of \code{\link{betterChromVAR}} for decorrelation in 
#' the deviation-based analysis).
#' 
#' The model is fitted separately for each motif and contrast, which takes 
#' longer than the analysis of deviations when there are many peaks and motifs.
#' 
#' This approach is related in spirit to motif activity response analysis 
#' (MARA, and its web implementation ISMARA), which models the signal of 
#' regulatory regions as a function of the motifs they harbor. It differs in 
#' that motifs are tested separately (MARA fits all motifs jointly), and 
#' that it works on the differential accessibility between groups, with a 
#' random effect for the heterogeneity across sites and a correction for the 
#' correlation between peaks.
#' 
#' @references
#'   The FANTOM Consortium and the RIKEN Omics Science Center (2009) The
#'   transcriptional network that controls growth arrest and differentiation 
#'   in a human myeloid leukemia cell line, Nature Genetics 41:553-562, 
#'   doi: 10.1038/ng.375 (motif activity response analysis, MARA).
#'   
#'   Balwierz P.J., Pachkov M., Arnold P., Gruber A.J., Zavolan M., 
#'   van Nimwegen E. (2014) ISMARA: automated modeling of genomic signals as a 
#'   democracy of regulatory motifs, Genome Research 24:869-884, 
#'   doi: 10.1101/gr.169508.113
#'   
#'   Wu D, Smyth GK (2012). Camera: a competitive gene set test accounting for 
#'   inter-gene correlation. Nucleic Acids Research 40, e133.
#'   doi: 10.1093/nar/gks461
#'
#' @returns A data.frame with, for each motif and contrast, the number of peaks 
#'   with the motif (`N`), the `difference` (log2-fold-change of the peaks with 
#'   the motif relative to the others), its standard error (`se`), t-value, 
#'   p-value and adjusted p-value, and the estimated heterogeneity (`tau`, 
#'   standard deviation of the true response across the sites of the motif) 
#'   and excess correlation (`rho`) of its peaks.
#' @author Pierre-Luc Germain
#' @importFrom stats p.adjust pt model.matrix ave
#' @export
#' @examples
#' attach(getDummyData())
#' counts$group <- rep(LETTERS[1:2], each=5)
#' # differential accessibility analysis of the peaks:
#' library(limma)
#' design <- model.matrix(~group, data=colData(counts))
#' y <- voom(assay(counts), design)
#' dea <- topTable(eBayes(lmFit(y, design)), coef=2, number=Inf, sort.by="none")
#' # peak-level motif tests:
#' peakLevelMotifTests(counts, dea, motifMatches, design=design, minCount=5)
peakLevelMotifTests <- function(se, dea, annotation, design=NULL, 
                                minCount=20, splineDF=3, adjustWidth=FALSE,
                                vif=TRUE, useAssay="counts", lib.size=NULL){
  vif <- .checkVif(vif, !is.null(design))
  inp <- .peakLevelInput(se, dea, annotation, design, splineDF, vif, 
                         useAssay, lib.size, adjustWidth)
  pl <- inp$pl
  annotation <- inp$annotation
  
  n <- as.numeric(Matrix::colSums(annotation))
  keep <- colnames(annotation)[which(n>=minCount & (nrow(annotation)-n)>=minCount)]
  if(length(keep)==0) stop("No motif passing filtering.")
  
  rhoCache <- new.env()
  res <- do.call(rbind, lapply(seq_along(inp$contrasts), \(ci){
    cbind(contrast=inp$contrasts[ci], 
          .peakLevelMotifRegress(pl, annotation, keep, ci, vif=vif, 
                                 rhoCache=rhoCache))
  }))
  .finalizePeakLevel(res)
}

# P-values, adjusted p-values and ordering of peak-level results
.finalizePeakLevel <- function(res){
  res$P.Value <- 2*pt(-abs(res$t), df=pmax(1, res$N-1))
  res$adj.P.Val <- ave(res$P.Value, res$contrast, 
                       FUN=\(p) p.adjust(p, "fdr"))
  res <- res[order(res$contrast, res$P.Value),]
  row.names(res) <- NULL
  res
}

.checkVif <- function(vif, hasDesign){
  if(is.numeric(vif)){
    stopifnot(length(vif)==1 && vif>=0 && vif<1)
    return(vif)
  }
  stopifnot(is.logical(vif) && length(vif)==1 && !is.na(vif))
  if(vif && !hasDesign){
    message("The correlation between peaks is estimated from the residuals of ",
            "the model of the differential analysis, which requires `design`; ",
            "no variance inflation is applied (see `vif`).")
    return(FALSE)
  }
  vif
}

# Reads the table(s) of results of a differential analysis. Returns a list of 
# data.frames with the columns lfc, se, and acc and bias (when available).
.parseDEATable <- function(tab){
  tab <- as.data.frame(tab)
  pick <- function(nms) { 
    w <- intersect(nms, colnames(tab))
    if(length(w)==0) NULL else w[1]
  }
  cl <- pick(c("logFC", "log2FoldChange"))
  if(is.null(cl)) stop("Could not find a log-fold-change column (`logFC` or ",
                       "`log2FoldChange`) in the table.")
  lfc <- as.numeric(tab[[cl]])
  cs <- pick(c("SE", "se", "lfcSE"))
  if(!is.null(cs)){
    se <- as.numeric(tab[[cs]])
  }else if(!is.null(ct <- pick(c("t", "stat")))){
    se <- abs(lfc/as.numeric(tab[[ct]]))
  }else if(!is.null(ct <- pick(c("F", "LR")))){
    se <- abs(lfc)/sqrt(as.numeric(tab[[ct]]))
  }else{
    stop("Could not find the standard error (`SE`, `se` or `lfcSE`) nor a ",
         "statistic to derive it from (`t`, `stat`, `F` or `LR`) in the table.")
  }
  ca <- pick(c("AveExpr", "logCPM", "Amean"))
  acc <- NULL
  if(!is.null(ca)){
    acc <- as.numeric(tab[[ca]])
  }else if("baseMean" %in% colnames(tab)){
    acc <- log2(as.numeric(tab$baseMean)+1)
  }
  bias <- if("bias" %in% colnames(tab)) as.numeric(tab$bias) else NULL
  list(lfc=lfc, se=se, acc=acc, bias=bias, rn=row.names(tab))
}

# Prepares the input of the peak-level analyses: the log-fold-changes and their
# variances (matched to the rows of the object), the adjustment covariates,
# and the data for the correlation between peaks.
.peakLevelInput <- function(se, dea, annotation, design, splineDF, vif, 
                            useAssay, lib.size, adjustWidth=FALSE){
  stopifnot(inherits(se, "SummarizedExperiment"))
  if(inherits(annotation, "SummarizedExperiment")){
    annotation <- assay(annotation)
  }
  annotation <- drop0(as(annotation>0, "CsparseMatrix"))
  if(nrow(se)!=nrow(annotation))
    stop("`se` and `annotation` should have the same number of rows.")
  if(!is.null(design)){
    design <- as.matrix(design)
    if(nrow(design)!=ncol(se))
      stop("`design` should have a row for each column of `se`.")
  }
  
  # table(s) of results
  tabs <- if(is.data.frame(dea) || is(dea, "DataFrame")) list(contrast=dea) else
    as.list(dea)
  if(length(tabs)==0) stop("`dea` is empty.")
  if(is.null(names(tabs))) names(tabs) <- paste0("contrast", seq_along(tabs))
  tabs <- lapply(tabs, .parseDEATable)
  
  # matching of the peaks
  ids <- row.names(se)
  rn <- tabs[[1]]$rn
  if(!all(vapply(tabs, \(x) identical(x$rn, rn), logical(1))))
    stop("All the tables should have the same row names.")
  if(!is.null(ids) && !is.null(rn)){
    idx <- match(ids, rn)
    if(all(is.na(idx))) 
      stop("The row names of `dea` do not match those of `se`.")
  }else{
    if(length(tabs[[1]]$lfc)!=nrow(se))
      stop("`dea` and `se` should have the same number of rows (or ",
           "have row names).")
    idx <- seq_len(nrow(se))
  }
  lfc <- matrix(vapply(tabs, \(x) x$lfc[idx], numeric(length(idx))), 
                nrow=length(idx), dimnames=list(NULL, names(tabs)))
  sdv <- matrix(vapply(tabs, \(x) x$se[idx], numeric(length(idx))), 
                nrow=length(idx), dimnames=list(NULL, names(tabs)))
  
  # log-counts-per-million (as in voom), for the correlation and accessibility
  acc <- tabs[[1]]$acc[idx]
  E <- NULL
  if(isTRUE(vif) || is.null(acc)){
    counts <- as.matrix(assay(se, useAssay))
    if(is.null(lib.size)) lib.size <- colSums(counts)
    E <- log2(t(t(counts+0.5)/(lib.size+1))*1e6)
    if(is.null(acc)) acc <- rowMeans(E)
  }
  bias <- tabs[[1]]$bias[idx]
  if(is.null(bias)) bias <- rowData(se)$bias
  if(is.null(bias))
    message("GC bias not included.",
      "Note that it is highly recommend to include it -- see ?addGCBias")
  ok <- rowSums(!is.finite(lfc) | !is.finite(sdv) | sdv<=0)==0 & is.finite(acc)
  if(!is.null(bias)) ok <- ok & is.finite(bias)
  if(!any(ok)) stop("No peak with valid values.")
  if(sum(!ok)>0) 
    message(sum(!ok), " peaks with missing or invalid values were ignored.")
  
  covs <- list()
  if(splineDF>0){
    covs$acc <- ns(acc[ok], df=splineDF)
    if(!is.null(bias) && length(unique(bias[ok]))>splineDF+1L)
      covs$bias <- ns(bias[ok], df=splineDF)
  }
  if(adjustWidth) covs$width <- .widthCovariate(se, ok, splineDF)
  covs <- do.call(cbind, covs)
  
  rhoAll <- 0
  if(isTRUE(vif)){
    E <- E[ok,,drop=FALSE]
    rhoAll <- max(0, interGeneCorrelation(E, design)$correlation)
  }
  pl <- list(v=if(isTRUE(vif)) list(E=E) else NULL, lfc=lfc[ok,,drop=FALSE], 
             vr=sdv[ok,,drop=FALSE]^2, covs=covs, rhoAll=rhoAll, design=design)
  list(pl=pl, annotation=annotation[ok,,drop=FALSE], contrasts=names(tabs))
}

# Excess correlation of the peaks with (index 2) and without (index 1) a motif
.motifRho <- function(pl, mm, vif){
  if(isFALSE(vif)) return(c(0,0))
  if(is.numeric(vif)) return(c(0,vif))
  vapply(1:2, \(k) .excessCorrelation(pl$v$E, pl$design, which(mm==(k-1L)), 
                                      pl$rhoAll), numeric(1))
}

# Regression of the log-fold-changes of the peaks of contrast `ci` on the 
# presence of each of the motifs, optionally adjusting for other motifs (`adj`, 
# a matrix with the presence of these other motifs in the peaks).
.peakLevelMotifRegress <- function(pl, annotation, motifs, ci, adj=NULL, 
                                   vif=FALSE, rhoCache=new.env()){
  do.call(rbind, lapply(motifs, \(m){
    mm <- as.vector(annotation[,m])
    if(is.null(rhoCache[[m]])) rhoCache[[m]] <- .motifRho(pl, mm, vif)
    rho <- rhoCache[[m]]
    X <- cbind(1, mm, adj, pl$covs)
    r <- tryCatch(.peakMetaRegression(pl$lfc[,ci], pl$vr[,ci], X, 1L+mm, rho),
                  error=function(e) NULL)
    if(is.null(r)){ # (e.g. a motif that is collinear with the adjustment)
      return(data.frame(motif=m, N=sum(mm), difference=NA_real_, se=NA_real_,
                        t=NA_real_, tau=NA_real_, rho=rho[2]))
    }
    sdb <- sqrt(diag(r$vcov))
    data.frame(motif=m, N=sum(mm), difference=r$beta[2], se=sdb[2], 
               t=r$beta[2]/sdb[2], tau=sqrt(r$tau2[2]), rho=rho[2])
  }))
}

#' stepwisePeakLevelMotifTests
#' 
#' Stepwise (forward) selection of the motifs associated with differential 
#' accessibility, using \code{\link{peakLevelMotifTests}}. In each round, the 
#' motifs are tested while adjusting for the motifs already selected; the 
#' strongest significant motif is then added to the selected ones, until no 
#' motif is significant anymore. This addresses the signal that spills over to 
#' motifs that share sites with the truly responsive ones.
#'
#' @inheritParams peakLevelMotifTests
#' @param alpha The threshold on the adjusted p-values (FDR across the motifs
#'   tested in the round) for a motif to be selected.
#' @param maxSteps The maximum number of motifs to select (per contrast).
#'
#' @details
#' In the first round, each motif is tested as in 
#' \code{\link{peakLevelMotifTests}}. If any motif has an adjusted p-value below
#' `alpha`, the one with the smallest p-value is selected, and its presence in 
#' each peak is added to the covariates of the regressions used to test the 
#' other motifs in the next round. The procedure is repeated, adding each time
#' the top motif, until no remaining motif is significant (or `maxSteps` is
#' reached). Each contrast is processed independently.
#' 
#' The effects of the selected motifs are then estimated jointly (all in the 
#' same model). Because the motif of a round is selected as the top among many, 
#' and the rounds are not independent, the p-values at selection are not 
#' calibrated error rates, and the procedure is heuristic; a conservative 
#' `alpha` is advisable. A motif that is (nearly) collinear with those 
#' selected cannot be tested, and is returned with missing values. Motifs are 
#' compared on their presence in peaks only, so motifs in strong linkage with
#' a responsive one can neither be distinguished from it nor from each other
#' if they share all the same sites. 
#' 
#' See \code{\link{peakLevelMotifTests}} for the other details of the model.
#'
#' @returns A list with two data.frames: `selected`, with one row for each 
#'   selected motif and contrast, in order of selection (`step`), with the 
#'   statistics of the motif at the time of its selection (i.e. adjusted for 
#'   the previously-selected motifs; see 
#'   \code{\link{peakLevelMotifTests}} for the columns), as well as the 
#'   estimates when all the selected motifs are fitted jointly 
#'   (`joint.difference`, `joint.se`, `joint.P.Value`); and `remaining`, with 
#'   the statistics of the other motifs in the last round (i.e. adjusted for all
#'   the selected motifs).
#' @author Pierre-Luc Germain
#' @importFrom stats p.adjust pt
#' @export
#' @examples
#' attach(getDummyData())
#' counts$group <- rep(LETTERS[1:2], each=5)
#' library(limma)
#' design <- model.matrix(~group, data=colData(counts))
#' y <- voom(assay(counts), design)
#' dea <- topTable(eBayes(lmFit(y, design)), coef=2, number=Inf, sort.by="none")
#' stepwisePeakLevelMotifTests(counts, dea, motifMatches, design=design, 
#'                             minCount=5)
stepwisePeakLevelMotifTests <- function(se, dea, annotation, design=NULL, 
                                        alpha=0.05, maxSteps=20L, minCount=20,
                                        splineDF=3, adjustWidth=FALSE, vif=TRUE,
                                        useAssay="counts", lib.size=NULL){
  stopifnot(is.numeric(alpha) && length(alpha)==1 && alpha>0 && alpha<1)
  stopifnot(is.numeric(maxSteps) && length(maxSteps)==1 && maxSteps>=1)
  vif <- .checkVif(vif, !is.null(design))
  inp <- .peakLevelInput(se, dea, annotation, design, splineDF, vif, 
                         useAssay, lib.size, adjustWidth)
  pl <- inp$pl
  annotation <- inp$annotation
  
  n <- as.numeric(Matrix::colSums(annotation))
  cand <- colnames(annotation)[which(n>=minCount & (nrow(annotation)-n)>=minCount)]
  if(length(cand)==0) stop("No motif passing filtering.")
  
  rhoCache <- new.env()
  out <- lapply(seq_along(inp$contrasts), \(ci){
    sel <- character()
    selRes <- list()
    remaining <- cand
    repeat{
      adj <- if(length(sel)==0) NULL else 
        as.matrix(annotation[,sel,drop=FALSE])*1
      r <- .peakLevelMotifRegress(pl, annotation, remaining, ci, adj=adj, 
                                  vif=vif, rhoCache=rhoCache)
      r <- .finalizePeakLevel(cbind(contrast=inp$contrasts[ci], r))
      top <- r[1,]
      if(is.na(top$adj.P.Val) || top$adj.P.Val>=alpha ||
         length(sel)>=maxSteps){
        break
      }
      sel <- c(sel, top$motif)
      selRes[[length(sel)]] <- cbind(top[,"contrast",drop=FALSE], 
                                     step=length(sel), top[,-1,drop=FALSE])
      remaining <- setdiff(remaining, top$motif)
      if(length(remaining)==0){
        r <- r[0,]
        break
      }
    }
    selRes <- if(length(selRes)==0) NULL else do.call(rbind, selRes)
    if(!is.null(selRes)) 
      selRes <- cbind(selRes, .jointPeakLevelFit(pl, annotation, sel, ci, vif))
    list(selected=selRes, remaining=r)
  })
  selected <- do.call(rbind, lapply(out, \(x) x$selected))
  if(is.null(selected)){
    selected <- data.frame(contrast=character(), step=integer(), 
                           motif=character(), N=integer(), 
                           difference=numeric(), se=numeric(), t=numeric(), 
                           tau=numeric(), rho=numeric(), P.Value=numeric(),
                           adj.P.Val=numeric(), joint.difference=numeric(), 
                           joint.se=numeric(), joint.P.Value=numeric())
  }
  row.names(selected) <- NULL
  list(selected=selected, remaining=do.call(rbind, lapply(out, \(x) x$remaining)))
}

# Joint fit of the selected motifs
.jointPeakLevelFit <- function(pl, annotation, sel, ci, vif){
  S <- as.matrix(annotation[,sel,drop=FALSE])*1
  anyS <- rowSums(S)>0
  X <- cbind(1, S, pl$covs)
  rho <- .motifRho(pl, anyS, vif)
  r <- tryCatch(.peakMetaRegression(pl$lfc[,ci], pl$vr[,ci], X, 1L+anyS, rho),
                error=function(e) NULL)
  k <- length(sel)
  if(is.null(r)) 
    return(data.frame(joint.difference=rep(NA_real_,k), joint.se=NA_real_, 
                      joint.P.Value=NA_real_))
  b <- r$beta[1+seq_len(k)]
  s <- sqrt(diag(r$vcov))[1+seq_len(k)]
  data.frame(joint.difference=b, joint.se=s, 
             joint.P.Value=2*pt(-abs(b/s), df=pmax(1, colSums(S)-1)))
}

# spline of the log-width of the peaks (NULL if they all have the same width)
.widthCovariate <- function(se, ok, splineDF){
  w <- if(methods::is(se, "RangedSummarizedExperiment")) 
    GenomicRanges::width(SummarizedExperiment::rowRanges(se)) else 
    rowData(se)$width
  if(is.null(w)) 
    stop("`adjustWidth` requires `rowRanges(se)` or a `width` column in ",
         "`rowData(se)`.")
  w <- log(as.numeric(w)[ok])
  if(length(unique(w))<=max(splineDF,1L)+1L){
    message("The peak widths are (nearly) constant; not adjusting for them.")
    return(NULL)
  }
  ns(w, df=max(splineDF, 1L))
}
