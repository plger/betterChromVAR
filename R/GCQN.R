#' GCQN: quantile normalization in GC bins
#'
#' Normalizes peak counts by quantile normalization within bins of GC content,
#' optionally with the group-aware smoothing of the quantile distributions
#' introduced in the `qsmooth` package. This is a reimplementation of the method
#' proposed by Van den Berge et al. (2022), and additionally includes the option
#' not to use GCQN without smoothing.
#'
#' If `grouping` is given, the reference of each sample is, for each quantile,
#' a weighted average of the global reference and of the reference of its
#' group, with a weight reflecting the (smoothed) proportion of the variation
#' across samples that is not explained by the groups, so that group
#' differences in the distribution are partly preserved ('smooth quantile
#' normalization'; Hicks et al., 2018).
#'
#' Unlike [CVnorm()], this approach operates on the whole distribution of
#' counts within a bin and does not use the accessibility of the peaks to define
#' the bins. As a consequence, a shift in the signal of many peaks of the same
#' bin (e.g. a set of differential peaks sharing the same GC content and
#' accessibility) is partially removed as if it were bias. See also
#' [testBiasDimensions()] to evaluate whether accessibility-dependent bias is
#' to be expected in a dataset.
#'
#' @param object A matrix of counts (peaks x samples), or a
#'   SummarizedExperiment-like object with a 'counts' assay.
#' @param bias A vector of length equal to `nrow(object)` specifying the
#'   per-peak GC content. If omitted, will try to get it from
#'   `rowData(object)$bias`.
#' @param grouping Optional grouping of the samples (a vector coercible to
#'   factor of length equal to `ncol(object)`, or the name of a `colData`
#'   column). If given (with at least two groups), the smooth quantile
#'   normalization is applied, preserving part of the differences between
#'   groups.
#' @param g The number of GC bins (default 20).
#' @param summary How to summarize the sorted counts across samples (or across
#'   the samples of a group) to obtain the reference distributions: either
#'   "mean" (default) or "median".
#' @param window The width of the running median used to smooth the weights
#'   across quantiles, as a fraction of the number of peaks in the bin. Ignored
#'   unless `grouping` is given.
#' @param ties How to handle tied counts within a sample: "min" gives them the
#'   reference value of their lowest rank, while "average" gives them the
#'   average of the reference values of the ranks they span. By default,
#'   "min" without `grouping` (as in the original GCQN implementation) and
#'   "average" with it (as in `qsmooth`).
#' @param round Logical; whether to round the normalized values to integers
#'   (default TRUE).
#' @param toAssay The name of the assay in which to store the normalized data
#'   (default 'GCQN'). Ignored unless `object` is a SummarizedExperiment-like
#'   object.
#'
#' @return If `object` is a matrix, then a matrix of normalized counts of the
#'   same dimensions. If `object` is a SummarizedExperiment-like object, then
#'   the object is returned with an extra assay named based on `toAssay`.
#' @references
#'   - Van den Berge K., Chou H.J., Roux de Bezieux H., Street K., Risso D., 
#'     Ngai J., Dudoit S. (2022) Normalization benchmark of ATAC-seq datasets
#'     shows the importance of accounting for GC-content effects, Cell Reports
#'     Methods 2(11):100321, doi: 10.1016/j.crmeth.2022.100321
#'   - Hicks SC, Okrah K, Paulson JN, Quackenbush J, Irizarry RA, Corrado
#'     Bravo H (2018). “Smooth quantile normalization.” Biostatistics 19 (2),
#'     doi: 10.1093/biostatistics/kxx028
#' @importFrom SummarizedExperiment assay assay<- rowData
#' @export
#' @examples
#' counts_se <- getDummyData()$counts
#' counts_se <- GCQN(counts_se)
GCQN <- function(object, bias=NULL, grouping=NULL, g=20L,
                 summary=c("mean","median"), window=0.05, ties=NULL,
                 round=TRUE, toAssay="GCQN"){
  summary <- match.arg(summary)
  if(.isSElike(object)){
    if(is.null(bias)) bias <- rowData(object)$bias
    counts <- assay(object, "counts")
  }else{
    counts <- object
  }
  if(is.null(bias)) stop("`bias` not provided, and not found in the object.")
  stopifnot(length(bias)==nrow(counts), !anyNA(bias))
  stopifnot(is.numeric(g) && length(g)==1 && g>=1)
  stopifnot(is.numeric(window) && length(window)==1 && window>0 && window<=1)
  counts <- as.matrix(counts)
  if(anyNA(counts)) stop("The counts must not contain missing values.")
  grouping <- .groupingInput(grouping, object)
  if(!is.null(grouping)){
    grouping <- factor(grouping)
    stopifnot(length(grouping)==ncol(counts))
    if(nlevels(droplevels(grouping))<2) grouping <- NULL
  }
  if(is.null(ties)) ties <- if(is.null(grouping)) "min" else "average"
  ties <- match.arg(ties, c("min","average"))

  # equal-count GC bins, on the inverse empirical CDF (ties are never split)
  cuts <- unique(stats::quantile(bias, (seq_len(g))/g, type=1, names=FALSE))
  bins <- pmin(length(cuts), findInterval(bias, cuts, left.open=TRUE)+1L)

  out <- matrix(NA_real_, nrow(counts), ncol(counts), dimnames=dimnames(counts))
  for(b in unique(bins)){
    i <- which(bins==b)
    out[i,] <- if(length(i)==1) counts[i,] else
      .qnBin(counts[i,,drop=FALSE], grouping, summary, window, ties)
  }
  if(isTRUE(round)) out <- round(out)
  out[out<0] <- 0

  if(!.isSElike(object)) return(out)
  assay(object, toAssay) <- out
  object
}

# quantile normalization of a (peaks x samples) matrix, with optional
# qsmooth-style smoothing of the reference towards the group references
.qnBin <- function(m, group=NULL, summary="median", window=0.05,
                   ties="min"){
  n <- nrow(m)
  S <- ncol(m)
  smry <- function(x) if(summary=="mean") rowMeans(x) else
    matrixStats::rowMedians(x)
  Q <- apply(m, 2, sort)
  Qref <- smry(Q)
  target <- matrix(Qref, n, S)
  if(!is.null(group)){
    Qg <- vapply(levels(group), function(g){
      i <- which(group==g)
      if(length(i)==0) return(Qref)
      smry(Q[,i,drop=FALSE])
    }, numeric(n))
    Qhat <- Qg[, as.integer(group), drop=FALSE]
    SST <- rowSums((Q-Qref)^2)
    SSB <- rowSums((Qhat-Qref)^2)
    w <- 1 - SSB/SST
    w[SST<1e-6] <- 1
    k <- floor(window*n)
    if(k%%2==0) k <- k+1L
    w <- stats::runmed(w, k=k, endrule="constant")
    target <- w*Qref + (1-w)*Qhat
  }
  out <- matrix(0, n, S, dimnames=dimnames(m))
  for(i in seq_len(S)){
    rmin <- rank(m[,i], ties.method="min")
    if(ties=="min"){
      out[,i] <- target[rmin,i]
    }else{
      rmax <- rank(m[,i], ties.method="max")
      cs <- c(0, cumsum(target[,i]))
      out[,i] <- (cs[rmax+1L]-cs[rmin])/(rmax-rmin+1L)
    }
  }
  out
}
