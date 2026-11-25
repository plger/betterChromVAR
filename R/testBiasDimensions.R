#' testBiasDimensions
#'
#' Tests, in a cross-validated way, which dimensions of peak-level bias (GC
#' content, accessibility and, if available, fragment-length bias) explain the
#' sample-specific biases of a (pseudo)bulk dataset, in the binning framework
#' used by [CVnorm()] and [getBackgroundBins()]. The result is meant to guide
#' the choice of normalization method (see Details).
#'
#' For each subset of the bias dimensions, peaks are binned along these
#' dimensions and, for each bin and sample, a multiplicative bias factor (the
#' ratio of the observed to the expected counts, with a pseudocount, as in
#' [CVnorm()]) is estimated on a training set of peaks. These factors are then
#' used to predict the counts of held-out peaks (the peak-wise expectation
#' being the pooled profile of the peak, multiplied by the library size and by
#' the bin factor), and the models are compared on the negative binomial
#' deviance of the held-out counts, relative to the deviance without
#' correction. Because the evaluation is out-of-sample, models with more bins
#' are not rewarded for overfitting, and when there is no bias the gain is zero
#' or negative. The held-out peaks are never used to estimate the factors that
#' predict them; by default folds are made of whole chromosomes.
#'
#' For each dimension, the *marginal* gain is that of the model using only that
#' dimension, and the *unique* gain is the loss of gain incurred by removing
#' that dimension from the full model (i.e. what it adds to the others). Since
#' bias dimensions are correlated, a dimension can have a large marginal but a
#' negligible unique gain, which means that it is redundant with the others. The
#' results are also given for the most accessible peaks, in which most of the
#' reads lie and where accessibility-dependent bias is often most relevant.
#'
#' Guidance on the choice of normalization (heuristic, see the `recommendation`
#' element of the output): if only GC content has a unique gain, methods based
#' on GC alone (e.g. GC-bin quantile normalization) should be sufficient and
#' the extra accessibility dimension of [CVnorm()] is unlikely to bring
#' anything; if accessibility (or fragment-length bias) has a non-negligible
#' unique gain, then a method that models it (such as [CVnorm()]) is expected
#' to outperform GC-only approaches. Note however that this test only
#' addresses the existence of bias: methods that use the data to define the
#' accessibility dimension can absorb (part of) strong, clustered differential
#' signal, so that `grouping` or `smoothGrouping` should be used in
#' [CVnorm()] when appropriate.
#'
#' The test is designed for (pseudo)bulk data, i.e. a moderate number of
#' well-covered samples. For single-cell data, aggregate the cells first into
#' pseudobulks (e.g. per library or sample, ideally with comparable
#' composition, since differences in cell-type composition are otherwise
#' mistaken for bias unless taken into account through `grouping`).
#'
#' @param object A matrix of counts (peaks x samples), or a
#'   SummarizedExperiment-like object with a 'counts' assay.
#' @param bias A vector of length equal to `nrow(object)` with the per-peak
#'   GC content. If omitted, will try to get it from `rowData(object)$bias`.
#' @param flbias Optional vector of the fragment-length bias of the peaks. If
#'   omitted, will be taken from `rowData(object)$flbias` if present.
#' @param grouping Optional grouping of the samples (vector of length
#'   `ncol(object)`, or the name of a `colData` column). If given, the
#'   expectation of each peak is computed per group (so that biological
#'   differences between groups are not counted as bias, and only the
#'   within-group variation is assessed).
#' @param nbins A named integer vector, giving the number of (approximately
#'   equal-count) bins for each of the dimensions (`accessibility`, `GC`,
#'   `flbias`). Defaults to `c(accessibility=20, GC=20, flbias=5)`. For
#'   accessibility, the top of the distribution is split more finely.
#' @param folds The number of cross-validation folds (default 5). If the object
#'   has genomic ranges for at least as many chromosomes, folds are made of
#'   chromosomes; otherwise peaks are assigned at random.
#' @param dispersion The dispersion of the negative binomial used to compute
#'   the deviance. If NULL (default), it is estimated from the data.
#' @param pseudoCount A pseudocount, in reads, used to compute the bias factors
#'   (see [CVnorm()]).
#' @param useWidthAdj Whether to adjust the accessibility dimension for the
#'   width of the regions (see [CVnorm()]). If NULL, determined as in
#'   [CVnorm()] when possible, and otherwise `FALSE`.
#' @param minShare The minimum share of the gain of the full model that a
#'   dimension should uniquely explain in order to be deemed relevant in the
#'   `recommendation`.
#' @param minGain The minimum gain (as a fraction of the deviance without
#'   correction) of the full model below which the bias is considered
#'   negligible in the `recommendation`.
#' @param seed Random seed used for random folds (the global seed is
#'   preserved).
#'
#' @return An object of class `biasDimTest`, a list with the elements:
#' \itemize{
#' \item{`models`: a data.frame with, for each model (subset of dimensions) and
#' each stratum of peaks (`all`, `low80` (least accessible 80% of peaks),
#' `top20`, `top5`), the `gain` (fraction of the deviance without correction
#' that is explained out-of-sample, pooled over folds) and its standard error
#' across folds (`gain.se`).}
#' \item{`importance`: a data.frame with, for each stratum and dimension, the
#' `marginal` gain, the `unique` gain (and its standard error across folds),
#' and the `share` of the full model's gain that is unique to the dimension.}
#' \item{`recommendation`: a character vector summarizing the interpretation.}
#' \item{`dispersion`, `folds`, `nbins`, `n`: the settings used.}
#' }
#' @author Pierre-Luc Germain
#' @seealso [CVnorm()], [getBackgroundBins()]
#' @importFrom SummarizedExperiment assay rowData rowRanges
#' @importFrom GenomicRanges seqnames width
#' @export
#' @examples
#' counts_se <- getDummyData()$counts
#' # (the dummy data have no real bias, so the gains are expected to be ~0)
#' res <- testBiasDimensions(counts_se, folds=3)
#' res
testBiasDimensions <- function(object, bias=NULL, flbias=NULL, grouping=NULL,
                               nbins=NULL, folds=5L, dispersion=NULL,
                               pseudoCount=1, useWidthAdj=NULL, minShare=0.1,
                               minGain=0.001, seed=1L){
  wi <- NULL
  if(.isSElike(object)){
    if(is.null(bias)) bias <- rowData(object)$bias
    if(is.null(flbias)) flbias <- rowData(object)$flbias
    counts <- assay(object, "counts")
    if(!is.null(rowRanges(object)) && length(rowRanges(object))==nrow(object))
      wi <- width(object)
  }else{
    counts <- object
  }
  if(is.null(bias)) stop("`bias` not provided, and not found in the object.")
  stopifnot(length(bias)==nrow(counts))
  if(!is.null(flbias)) stopifnot(length(flbias)==nrow(counts))
  if(ncol(counts)<3) stop("At least 3 samples are required.")
  grouping <- .groupingInput(grouping, object)
  if(!is.null(grouping)) grouping <- factor(grouping)
  counts <- as.matrix(counts)
  storage.mode(counts) <- "double"

  # chromosomes, for the folds
  chrom <- NULL
  if(.isSElike(object) && !is.null(rowRanges(object)) &&
     length(rowRanges(object))==nrow(counts))
    chrom <- as.character(seqnames(rowRanges(object)))

  keep <- rowSums(counts)>0 & !is.na(bias) &
    (is.null(flbias) | !is.na(flbias))
  if(!all(keep)){
    message(sum(!keep), " peaks with no reads or with missing bias were ",
            "excluded.")
    counts <- counts[keep,,drop=FALSE]
    bias <- bias[keep]
    if(!is.null(flbias)) flbias <- flbias[keep]
    if(!is.null(wi)) wi <- wi[keep]
    if(!is.null(chrom)) chrom <- chrom[keep]
  }
  N <- nrow(counts)
  S <- ncol(counts)
  cs <- colSums(counts)

  # baseline expectation
  if(is.null(grouping)){
    q <- rowSums(counts)/sum(counts)
    mu0 <- outer(q, cs)
    ex <- rowMeans(counts)
  }else{
    mu0 <- matrix(0, N, S)
    for(g in levels(grouping)){
      i <- which(grouping==g)
      if(length(i)==0) next
      qg <- rowSums(counts[,i,drop=FALSE])/sum(counts[,i,drop=FALSE])
      mu0[,i] <- outer(qg, cs[i])
    }
    ex <- getExpectation(counts, grouping)
    q <- ex/sum(ex)
  }
  if(is.null(useWidthAdj))
    useWidthAdj <- !is.null(wi) && (mean(abs(wi-stats::median(wi)))/stats::median(wi)) > 0.1
  if(isTRUE(useWidthAdj)) useWidthAdj <- 200L
  if(!isFALSE(useWidthAdj)){
    if(is.null(wi)) stop("`useWidthAdj` requires an object with rowRanges.")
    ex <- useWidthAdj*ex/pmax(wi, useWidthAdj)
  }

  # dimensions
  dims <- list(accessibility=log10(ex), GC=bias)
  if(!is.null(flbias)) dims$flbias <- flbias
  defNb <- c(accessibility=20L, GC=20L, flbias=5L)
  if(!is.null(nbins)){
    stopifnot(!is.null(names(nbins)) && all(names(nbins) %in% names(defNb)))
    defNb[names(nbins)] <- as.integer(nbins)
  }
  nbins <- defNb[names(dims)]
  bins <- lapply(names(dims), function(d){
    .quantBins(dims[[d]], nbins[[d]], fineTop=(d=="accessibility"))
  })
  names(bins) <- names(dims)

  # strata
  acr <- rank(ex, ties.method="first")/N
  strata <- list(all=rep(TRUE,N), low80=acr<=0.8, top20=acr>0.8, top5=acr>0.95)

  # models: all non-empty subsets of the dimensions
  D <- length(dims)
  subsets <- unlist(lapply(seq_len(D), function(k)
    utils::combn(names(dims), k, simplify=FALSE)), recursive=FALSE)
  mname <- function(s) paste(s, collapse="+")
  cellOf <- function(s){
    id <- rep(1L, N)
    for(d in s) id <- (id-1L)*max(bins[[d]]) + bins[[d]]
    id
  }
  cells <- lapply(subsets, cellOf)
  names(cells) <- vapply(subsets, mname, character(1))

  # dispersion
  if(is.null(dispersion)){
    ok <- which(rowMeans(counts)>=5)
    if(length(ok)<100) ok <- seq_len(N)
    # the expectation is fitted on the same counts: correct the degrees of 
    # freedom (the estimate remains an upper bound, as bias is not modelled)
    G <- if(is.null(grouping)) 1L else nlevels(droplevels(grouping))
    r <- (rowSums((counts[ok,,drop=FALSE]-mu0[ok,,drop=FALSE])^2)*S/(S-G) -
            rowSums(mu0[ok,,drop=FALSE]))/rowSums(mu0[ok,,drop=FALSE]^2)
    dispersion <- min(2, max(0.01, stats::median(r, na.rm=TRUE)))
  }
  stopifnot(is.numeric(dispersion) && length(dispersion)==1 && dispersion>0)

  # folds
  folds <- as.integer(min(folds, N))
  stopifnot(folds>=2)
  fold <- .makeFolds(N, folds, chrom, seed)
  folds <- length(unique(fold))

  # held-out negative binomial deviance
  nbdev <- function(y, mu){
    a <- 1/dispersion
    mu <- pmax(mu, 1e-8)
    t1 <- ifelse(y>0, y*log(y/mu), 0)
    2*(t1 - (y+a)*log((y+a)/(mu+a)))
  }
  Dsum <- function(dv, idx) vapply(strata, function(st)
    sum(rowSums(dv)[st[idx]]), numeric(1))

  # D[fold, model, stratum]
  dev0 <- array(NA_real_, c(folds, length(strata)),
                dimnames=list(NULL, names(strata)))
  devm <- array(NA_real_, c(folds, length(cells), length(strata)),
                dimnames=list(NULL, names(cells), names(strata)))
  for(f in seq_len(folds)){
    te <- which(fold==f)
    tr <- which(fold!=f)
    yte <- counts[te,,drop=FALSE]
    dev0[f,] <- Dsum(nbdev(yte, mu0[te,,drop=FALSE]), te)
    for(m in names(cells)){
      ce <- cells[[m]]
      obs <- rowsum(counts[tr,,drop=FALSE], ce[tr])
      exp <- rowsum(mu0[tr,,drop=FALSE], ce[tr])
      R <- (obs+pseudoCount)/(exp+pseudoCount)
      j <- match(ce[te], as.integer(rownames(R)))
      Rte <- R[j,,drop=FALSE]
      Rte[is.na(j),] <- 1
      devm[f,m,] <- Dsum(nbdev(yte, mu0[te,,drop=FALSE]*Rte), te)
    }
  }
  gainf <- function(f, m, st) 1 - devm[f,m,st]/dev0[f,st]
  pooled <- function(m, st) 1 - sum(devm[,m,st])/sum(dev0[,st])

  models <- do.call(rbind, lapply(names(strata), function(st){
    do.call(rbind, lapply(names(cells), function(m){
      g <- vapply(seq_len(folds), gainf, numeric(1), m=m, st=st)
      data.frame(stratum=st, model=m, gain=pooled(m,st),
                 gain.se=stats::sd(g)/sqrt(folds), stringsAsFactors=FALSE)
    }))
  }))
  full <- names(cells)[length(cells)]
  importance <- do.call(rbind, lapply(names(strata), function(st){
    do.call(rbind, lapply(names(dims), function(d){
      rest <- setdiff(names(dims), d)
      gfull <- pooled(full, st)
      if(length(rest)==0){
        un <- gfull
        uf <- vapply(seq_len(folds), gainf, numeric(1), m=full, st=st)
      }else{
        mr <- mname(rest)
        un <- gfull - pooled(mr, st)
        uf <- vapply(seq_len(folds), function(f)
          gainf(f, full, st) - gainf(f, mr, st), numeric(1))
      }
      data.frame(stratum=st, dimension=d, marginal=pooled(d, st), unique=un,
                 unique.se=stats::sd(uf)/sqrt(folds),
                 share=if(gfull>0) un/gfull else NA_real_,
                 stringsAsFactors=FALSE)
    }))
  }))
  rownames(models) <- rownames(importance) <- NULL

  res <- list(models=models, importance=importance, dispersion=dispersion,
              folds=folds, nbins=nbins, n=c(peaks=N, samples=S))
  res$recommendation <- .biasDimAdvice(res, minShare, minGain)
  class(res) <- "biasDimTest"
  res
}

# approximately equal-count bins (integer vector of the same length as x);
# if fineTop, the top of the distribution is split more finely
.quantBins <- function(x, nb, fineTop=FALSE){
  nb <- as.integer(nb)
  if(fineTop && nb>=8){
    probs <- c(seq(0, 0.9, length.out=nb-3L), 0.97, 0.995, 1)
  }else{
    probs <- seq(0, 1, length.out=nb+1L)
  }
  brk <- unique(stats::quantile(x, probs, na.rm=TRUE))
  if(length(brk)<2) return(rep(1L, length(x)))
  as.integer(cut(x, brk, include.lowest=TRUE, labels=FALSE))
}

# assign peaks to folds, by chromosome when possible
.makeFolds <- function(N, k, chrom=NULL, seed=1L){
  if(!is.null(chrom) && length(unique(chrom))>=k){
    # greedy balancing of the chromosomes
    sz <- sort(table(chrom), decreasing=TRUE)
    load <- rep(0, k)
    asg <- stats::setNames(integer(length(sz)), names(sz))
    for(ch in names(sz)){
      f <- which.min(load)
      asg[ch] <- f
      load[f] <- load[f] + sz[[ch]]
    }
    return(as.integer(asg[chrom]))
  }
  if(!is.null(seed)){
    if(exists(".Random.seed", envir=globalenv())){
      old <- get(".Random.seed", envir=globalenv())
      on.exit(assign(".Random.seed", old, envir=globalenv()))
    }else{
      on.exit(rm(".Random.seed", envir=globalenv()))
    }
    set.seed(seed)
  }
  sample(rep_len(seq_len(k), N))
}

.biasDimAdvice <- function(res, minShare, minGain){
  imp <- res$importance[res$importance$stratum=="all",]
  imt <- res$importance[res$importance$stratum=="top20",]
  fm <- res$models[res$models$model==paste(unique(imp$dimension),
                                           collapse="+") &
                     res$models$stratum=="all",]
  if(fm$gain < minGain)
    return(paste0("Sample-specific bias along these dimensions is negligible ",
                  "(held-out gain of the full model: ",
                  round(100*fm$gain,2), "%): normalization is unlikely to be ",
                  "needed."))
  relevant <- function(i) !is.na(i$share) & i$share>=minShare &
    i$unique > 2*i$unique.se
  rel <- imp$dimension[relevant(imp)]
  relTop <- imt$dimension[relevant(imt)]
  out <- paste0("The full model explains ", round(100*fm$gain,2),
                "% of the deviance of held-out peaks. Dimensions adding to ",
                "the others (unique gain >= ", 100*minShare,
                "% of the total and >2 SE): ",
                if(length(rel)) paste(rel, collapse=", ") else "none", ".")
  if(length(rel)==0){
    out <- c(out, "No single dimension is uniquely informative; the bias may be shared across dimensions (redundant) or weak.")
  }else if(identical(rel, "GC")){
    out <- c(out, paste0("The bias is predominantly GC-dependent: methods ",
             "based on GC alone (e.g. GC-bin quantile normalization) should ",
             "be sufficient, and the accessibility dimension of CVnorm is ",
             "unlikely to bring much."))
  }else if("accessibility" %in% rel){
    out <- c(out, paste0("Accessibility-dependent bias exists beyond what is ",
             "explained by the other dimension(s): approaches correcting only ",
             "GC-dependent shifts would leave it in the data, and a method ",
             "modelling it (e.g. CVnorm) is expected to be preferable."))
  }
  if("flbias" %in% rel)
    out <- c(out, paste0("Fragment-length bias also contributes: ensure that ",
                         "it is accounted for (`flbias` in rowData for ",
                         "CVnorm)."))
  if(length(setdiff(relTop, rel)) || length(setdiff(rel, relTop)))
    out <- c(out, paste0("In the 20% most accessible peaks, the relevant ",
             "dimensions are: ", if(length(relTop)) paste(relTop,
             collapse=", ") else "none", "."))
  out
}

#' @export
print.biasDimTest <- function(x, ...){
  cat("Cross-validated test of the bias dimensions (", x$n[["peaks"]],
      " peaks, ", x$n[["samples"]], " samples, ", x$folds, " folds, NB dispersion ",
      signif(x$dispersion,2), ")\n", sep="")
  cat("\nImportance (gain = fraction of the held-out deviance explained):\n")
  imp <- x$importance
  imp[,c("marginal","unique","unique.se","share")] <-
    round(imp[,c("marginal","unique","unique.se","share")], 4)
  print(imp, row.names=FALSE)
  cat("\n")
  cat(paste(strwrap(x$recommendation, exdent=2), collapse="\n"), "\n")
  invisible(x)
}
