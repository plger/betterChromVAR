#' estimateNoiseKappa
#'
#' Empirically calibrates, for each motif, the analytic estimate of the noise
#' variance of the single-cell deviations, and estimates the reliability of the
#' single-cell deviations of each motif. This is mostly for development and 
#' testing: for practical purposes the analytic estimates are sufficient.
#'
#' @param object A SummarizedExperiment (or SingleCellExperiment) with an assay
#'    'counts', and with a 'bias' column in `rowData(object)` (as for
#'    \code{\link{betterChromVAR}}), or a (sparse) matrix of counts, in which
#'    case `background` is required.
#' @param annotations Peak annotation (sparse) matrix, with motifs as columns,
#'    or a SummarizedExperiment containing this in the first assay.
#' @param expectation Optional vector of expected counts per peak (as used for
#'    the computation of the deviations). If NULL, taken from `background` if
#'    available, otherwise computed from `object` using
#'    \code{\link{getExpectation}} (with `grouping`).
#' @param background Optional `bcvBackground` object as produced by
#'   \code{\link{getBackgroundBins}}. If NULL, created from the expectation and
#'   the bias in `rowData(object)`.
#' @param grouping Optional grouping of the cells, used to compute the
#'   expectation (see \code{\link{betterChromVAR}}); ignored if `expectation`
#'   is given.
#' @param ncells The (maximum) number of randomly-selected cells to use
#'   (default 3000). More cells increase the precision of the estimates (and
#'   the running time).
#' @param pairs Logical; whether counts are Tn5 insertions (e.g. ArchR peak
#'   matrices), in which case a fragment with both ends in a peak contributes 2
#'   and the read-splitting needs to keep such pairs together. If NULL
#'   (default), this is guessed from the data (pairs are assumed if more than
#'   half of all counts are in entries with an even value).
#' @param fpc Logical; whether the finite-population correction was applied to
#'   the deviations whose noise is to be calibrated (see
#'   \code{\link{computeDeviationsAnalytic}}). This should match the setting
#'   used for the deviations passed to \code{\link{smoothDeviationsKNN}}.
#' @param nthreads Either an integer scalar indicating the number of threads to
#'   use, or a `BiocParallelParam` object (default NULL, i.e. no
#'   parallelization). The cells are processed in chunks of `chunkSize`.
#' @param chunkSize Number of cells processed at a time (per worker). Lower
#'   values reduce memory usage.
#' @param cellWeights Logical; whether to weigh the cells by the precision of
#'   their deviations (default TRUE) rather than equally. Precision weighting
#'   prevents the (many) shallow, very noisy cells from dominating the
#'   estimates.
#' @param nBlocks The number of blocks of cells used for the (jackknife)
#'   estimation of the standard error of the reliability.
#' @param verbose Logical; whether to output progress messages.
#' @param ... Passed to \code{\link{getBackgroundBins}} (when `background` is
#'   NULL).
#'
#' @details
#' The variance of the deviations computed analytically overestimates the noise
#' for some motifs. This function first computes the proportion of the variance
#' that is due to noise in each motif, using the difference between the 
#' deviations computed from two independent halves of the reads/fragments. This 
#' provides two pieces of information: first, it allows estimating, for each 
#' motif, the reliability of its single-cell deviations, understood as the 
#' proportion of their variance across cells that is not noise. Second, it 
#' enables computing a factor, `kappa`, that calibrates the analytic variance 
#' computed from the deviations (a value close to 1 means that the analytic 
#' variance is accurate).
#'
#' The selection of the cells and the splitting of their reads are random: use
#' `set.seed` before calling the function for reproducible results (if you
#' provide a `BiocParallelParam` as `nthreads`, set its `RNGseed`).
#'
#' @return A data.frame with one row per motif and the columns `kappa` (to be
#'   given to \code{\link{smoothDeviationsKNN}}), `reliability` and
#'   `reliability.se` (a jackknife estimate of the standard error of the 
#'   reliability).
#'
#' @seealso \code{\link{smoothDeviationsKNN}}
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam bpnworkers
#' @importFrom stats rbinom
#' @export
#' @examples
#' attach(getDummyData())
#' res <- estimateNoiseKappa(counts, motifMatches, pairs=FALSE)
#' res
estimateNoiseKappa <- function(object, annotations, expectation=NULL,
                               background=NULL, grouping=NULL, ncells=3000L,
                               pairs=NULL, fpc=TRUE, nthreads=NULL,
                               chunkSize=500L, cellWeights=TRUE,
                               nBlocks=20L, verbose=FALSE, ...){
  if(.isSElike(object)){
    counts <- assay(object, "counts")
    bias <- rowData(object)$bias
    flbias <- rowData(object)$flbias
  }else{
    counts <- object
    bias <- flbias <- NULL
    if(is.null(background))
      stop("`background` is required when `object` is not a ",
           "SummarizedExperiment.")
  }
  if(!is(counts, "matrix") && !inherits(counts, "Matrix") &&
     !is(counts, "DelayedMatrix"))
    stop("`object` should be a SummarizedExperiment or SingleCellExperiment,",
         " or a (sparse) matrix of counts.")
  if(.isSElike(annotations)) annotations <- assay(annotations)
  .checkAnnotations(annotations)
  stopifnot(nrow(counts)==nrow(annotations))
  stopifnot(length(ncells)==1 && ncells>=2 && length(chunkSize)==1 &&
              chunkSize>=1)
  stopifnot(is.logical(cellWeights) && length(cellWeights)==1 &&
              !is.na(cellWeights))
  stopifnot(length(nBlocks)==1 && nBlocks>=2)
  BPPARAM <- .bpFromThreads(nthreads, progressbar=verbose)

  if(is.null(expectation)){
    if(!is.null(background) && length(background@expectation)==nrow(counts)){
      expectation <- background@expectation
    }else{
      expectation <- getExpectation(counts, .groupingInput(grouping, object))
    }
  }
  stopifnot(length(expectation)==nrow(counts))
  if(is.null(background)){
    if(is.null(bias)) stop("`bias` not found in the object.")
    background <- getBackgroundBins(expectation, bias=bias, flbias=flbias,
                                    verbose=verbose, ...)
  }
  stopifnot(is(background, "bcvBackground"))
  # discard any sample-level information from the background
  background@E <- NULL
  background@V <- NULL
  background@depth <- numeric()

  # choice of the cells
  ncells <- as.integer(min(ncells, ncol(counts)))
  idx <- sort(sample.int(ncol(counts), ncells))
  nBlocks <- as.integer(min(nBlocks, ncells))
  blocks <- (seq_along(idx)-1L) %% nBlocks + 1L   # jackknife blocks
  chunkId <- ceiling(seq_along(idx)/chunkSize)
  chunks <- split(idx, chunkId)
  blocks <- split(blocks, chunkId)

  if(is.null(pairs)){
    x <- .asDgC(counts[, chunks[[1]], drop=FALSE])@x
    pairs <- sum(x[x%%2==0]) > sum(x)/2
    message("Counts assumed to ",
            ifelse(pairs, "be (paired) insertions", "not be paired insertions"),
            " (`pairs=", pairs, "`).")
  }

  if(verbose) message("Estimating noise calibration on ", ncells, " cells")
  res <- bplapply(seq_along(chunks), BPPARAM=BPPARAM, function(i){
    .kappaChunk(counts[, chunks[[i]], drop=FALSE], annotations, background,
                expectation, pairs=pairs, fpc=fpc,
                blocks=blocks[[i]], nBlocks=nBlocks,
                weighted=cellWeights)
  })
  # sums per motif and jackknife block
  st <- lapply(names(res[[1]]), function(n) Reduce("+", lapply(res, `[[`, n)))
  names(st) <- names(res[[1]])
  full <- .kappaStats(lapply(st, rowSums))
  jk <- vapply(seq_len(nBlocks), FUN.VALUE=numeric(nrow(st$sw)), function(b){
    .kappaStats(lapply(st, function(x) rowSums(x[,-b,drop=FALSE])))$rel
  })
  data.frame(kappa=full$kappa, reliability=pmax(0, full$rel),
             reliability.se=sqrt((nBlocks-1)/nBlocks *
                                   rowSums((jk-rowMeans(jk))^2)),
             row.names=colnames(annotations))
}

# kappa and (untruncated) reliability from the (weighted) sums
.kappaStats <- function(s){
  mA <- s$swa/s$sw
  mB <- s$swb/s$sw
  S <- s$swab/s$sw - mA*mB          # variance of the true activity
  N <- s$sdd/s$sw/4                 # noise variance with all the reads
  list(kappa=s$sdd/s$swv, rel=S/(S+N))
}

# NOTE: This function is heavy Claude, the output makes sense, but I haven't 
# checked line by line. Currently used for testing only.
# Processes one chunk of cells and returns, for each motif (rows) and
# jackknife block (columns), the weighted sums from which kappa and the
# reliability are computed (see `.kappaStats`). With dA and dB the deviations
# computed from the two halves of the reads of a cell, v the sum of their
# analytic noise variances, and w the weight of the cell for that motif:
#   sw   = sum(w)               swa = sum(w*dA)          swb = sum(w*dB)
#   swab = sum(w*dA*dB)         sdd = sum(w*(dA-dB)^2)   swv = sum(w*v)
# Returning sums (rather than the final statistics) lets the chunks, which can
# be processed in parallel, be simply added up, and allows leaving out one block
# of cells at a time (jackknife) without recomputing anything.
.kappaChunk <- function(counts, annotations, background, expectation, pairs,
                        fpc, blocks, nBlocks, weighted=TRUE){
  # randomly split the reads of each cell into two halves, A and B
  sp <- .thinCounts(.asDgC(counts), 0.5, pairs)
  rm(counts)
  # deviations of each half (same bins and expectation as for the full data)
  devs <- lapply(sp, function(x){
    # each half has about half of the reads, so the expectation is halved too
    bg <- computeBackgrounds(x, background, expectation=expectation/2)
    r <- computeDeviationsAnalytic(x, bg, annotations, retSE=FALSE, fpc=fpc,
                                   compute=c("deviations","z"))
    # d: deviations; v: their analytic noise variance. Since deviations is
    # (obs-exp)/denominator and z is (obs-exp)/null sd, deviations/z is the
    # null sd divided by the denominator, i.e. the noise sd on the scale of the
    # deviations.
    list(d=r$deviations, v=(r$deviations/r$z)^2)
  })
  rm(sp)
  # both halves are noisy, so the variance of their difference is the sum of
  # their variances
  v <- devs[[1]]$v + devs[[2]]$v
  dA <- devs[[1]]$d
  dB <- devs[[2]]$d
  rm(devs)
  # entries for which the noise variance (or the deviations) is not defined
  # (e.g. zero z-score) are given a weight of 0, i.e. ignored
  ok <- is.finite(v) & v>0 & is.finite(dA) & is.finite(dB)
  # weights: precision (1/variance) of the pair of deviations, or equal
  w <- if(weighted) 1/v else array(1, dim(v))
  w[!ok] <- 0
  v[!ok] <- 0
  dA[!ok] <- 0
  dB[!ok] <- 0
  # indicator matrix (cells x blocks): multiplying a (motifs x cells) matrix
  # by it sums the values of the cells of each jackknife block
  B <- matrix(0, nrow=length(blocks), ncol=nBlocks)
  B[cbind(seq_along(blocks), blocks)] <- 1
  list(sw=w %*% B, swa=(w*dA) %*% B, swb=(w*dB) %*% B, swab=(w*dA*dB) %*% B,
       sdd=(w*(dA-dB)^2) %*% B, swv=(w*v) %*% B)
}

# Binomial thinning of a count matrix (dgCMatrix) in two complementary halves.
# If `pairs`, even counts are thinned as units of two (e.g. both ends of a
# fragment falling in the same peak).
.thinCounts <- function(x, p, pairs=TRUE){
  v <- x@x
  if(pairs){
    u <- floor(v/2)
    a <- 2*rbinom(length(v), u, p) + rbinom(length(v), v-2*u, p)
  }else{
    a <- rbinom(length(v), v, p)
  }
  A <- x
  A@x <- as.numeric(a)
  x@x <- v - A@x
  list(A=drop0(A), B=drop0(x))
}