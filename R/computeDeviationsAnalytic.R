#' computeDeviationsAnalytic
#'
#' @param object A SummarizedExperiment (or SingleCellExperiment) with an assay
#'    'counts', and with a 'bias' column in `rowData(object)`. Note that the 
#'    regions should have similar widths.
#' @param annotations Peak annotation (sparse) matrix, with motifs as columns,
#'    or a SummarizedExperiment containing this in the first assay. Values 
#'    should be either logical or between 0 and 1.
#' @param background A `bcvBackground` object with bins*samples slots filled, 
#'   as produced by \code{\link{computeBackgrounds}}.
#' @param verbose Logical; whether to output progress messages.
#' @param retSE Logical; whether to return a SummarizedExperiment object.
#' @param compute What to compute. Defaults to everything: deviations, z and 
#'   motif variability.
#' @param denominator The type of denominator to use for the deviations.
#'   Either 'global' (default), i.e. the global expectation (same as the 
#'   original chromVAR), 'local' (background expectation of the cell/sample), 
#'   or 'none' (denominator of 1). 'global' (default) is recommended.
#' @param fpc Logical, whether to apply finite population correction (see
#'   details below). Disable (`fpc=FALSE`) to reproduce the behavior of the 
#'   original chromVAR.
#' @param ridge An optional positive number, enabling the decorrelation of 
#'   motifs that share peaks (see details). Smaller values decorrelate more 
#'   strongly; values of 1 to 3 are recommended. If NULL (default) or `Inf`, 
#'   motifs are not decorrelated.
#'
#' @details
#' `fpc` corrects correct the null variance used (as denominator) for the 
#' z-scores for the fact that the expectation of each cell is estimated from 
#' the cell's own counts, including in peaks containing the motif. This 
#' correlation between observed and expected counts reduces the variance of 
#' their difference in a way that depends on the proportion of peaks matching
#' the motif. This is negligilbe for motifs matching few peaks, and especially
#' affects prevalent motifs concentrated in a few background bins (e.g. GC-rich
#' motifs). The correction typically has little effect on cell ranking within a
#' motif.
#' 
#' Since peaks carry multiple motifs and some motif families are very similar,
#' the deviations of motifs tend to be correlated, sometimes highly so. If 
#' `ridge` is given, the excess counts of the motifs are deconvoluted using the 
#' overlap between them (after background subtraction, and weighted by the 
#' expected accessibility of the peaks), with a ridge penalty proportional to 
#' `ridge` and to the size of each motif. Motifs that do not share peaks with 
#' others are unaffected, and the signal of nearly identical motifs is split 
#' between them. Note that this only removes what is due to shared peaks, and 
#' that strong decorrelation (e.g. `ridge<1`) reduces the signal of the motifs. 
#' The decorrelation additionally benefits from annotations that also capture 
#' weak matches (see \code{\link{getGradedMatches}}).
#'
#' @returns A SummarizedExperiment (or a list if `retSE=FALSE`).
#' @importFrom Matrix t crossprod tcrossprod
#' @export
#'
#' @examples
#' attach(getDummyData())
#' # if GC content not already in the object, use:
#' # counts <- addGCBias(counts, genome=YOUR_GENOME)
#' 
#' # we fist get the background bins:
#' bg <- getBackgroundBins(counts)
#' # then we compute the backgrounds for each sample:
#' bg <- computeBackgrounds(counts, bg)
#' # then we can compute the deviations:
#' dev <- computeDeviationsAnalytic(counts, bg, motifMatches)
#' dev
computeDeviationsAnalytic <- function(object, background, annotations,
                                      verbose=FALSE, retSE=TRUE,
                                      compute=c("deviations","z","variability"),
                                      denominator=c("global","local","none"),
                                      fpc=TRUE, ridge=NULL){
  # Check input validity
  stopifnot(nrow(object)==nrow(annotations))
  compute <- match.arg(compute, several.ok=TRUE)
  denominator <- match.arg(denominator)
  stopifnot(is(background, "bcvBackground"))

  if(length(background@depth)==0 || length(background@expectation)==0)
    stop("Incomplete background; please run computeBackgrounds() first.")
  
  if( (length(background@depth)!=ncol(object)) ||
      (length(background@expectation)!=nrow(object)) )
    stop("The `background` object does not match the dimensions of `object`.")

  motifCD <- depth <- NULL
  if(.isSElike(object)){
    counts <- assay(object, "counts")
    if(is.numeric(object$depth)){
      if(verbose) message("Using pre-computed object$depth")
      depth <- object$depth
    }
  }else{
    object <- SummarizedExperiment(list(counts=object))
    counts <- assay(object)
  }
  if(!is(counts, "matrix") && !inherits(counts, "Matrix"))
    stop("`object` should be a SummarizedExperiment or SingleCellExperiment,",
         " or a (sparse) matrix of counts.")
  
  if(.isSElike(annotations)){
    motifCD <- colData(annotations)
    annotations <- assay(annotations)
  } 
  .checkAnnotations(annotations)

  st <- .annotationStats(annotations, background)
  a <- .computeDeviations(counts, background, annotations, st=st, depth=depth,
                          compute=compute, denominator=denominator, fpc=fpc,
                          ridge=ridge, verbose=verbose)

  if(!isTRUE(retSE)) return(a)

  d <- data.frame(N=Matrix::colSums(annotations),
                  Neff=.effectiveN(annotations, background,
                                   background@expectation, fpc, st=st),
                  total=a$total)

  if("variability" %in% compute){
    d <- cbind(d, computeMotifVariability(a$z))
  }

  .packageDevSE(a, object, motifCD, d)
}

.computeDeviations <- function(counts, background, annotations, st,
                               depth=NULL,
                               compute=c("deviations","z","variability"),
                               denominator="global", fpc=TRUE, ridge=NULL,
                               verbose=FALSE){
  if(is.null(depth)) depth <- Matrix::colSums(counts)

  # motif-level background expectation and observed sums (M x S)
  motif_bg_exp <- as.matrix(Matrix::crossprod(st$N, background@E))
  observed_motif_sum <- .annotCrossprod(annotations, counts)

  # deviation = (Obs - bgExpect) / bgExpect; z = (Obs-exp)/sdExpect
  deviations <- observed_motif_sum - motif_bg_exp

  if(!is.null(ridge))
    ridge <- .ridgeOperator(annotations, background, ridge, st=st)
  if(!is.null(ridge)){
    # remove from each motif what is attributable to overlapping motifs
    if(verbose) message("Decorrelating motifs")
    deviations <- matrix(ridge$L %*% deviations, nrow=nrow(deviations),
                         dimnames=dimnames(deviations))
  }

  a <- list()
  if("variability" %in% compute) compute <- union(compute, "z")
  if(any("z" %in% compute)){
    # z = (Obs-exp)/sdExpect
    if(isTRUE(fpc)){
      motif_bg_var <- .fpcVariance(annotations, background, background@V,
                                   st=st)
    }else{
      motif_bg_var <- as.matrix(Matrix::crossprod(st$N, background@V))
    }
    if(!is.null(ridge)) motif_bg_var <- motif_bg_var * ridge$varRatio
    motif_bg_sd <- sqrt(pmax(0, motif_bg_var))
    a$z <- deviations / motif_bg_sd
  }
  if("deviations" %in% compute){
    if(denominator=="local"){
      # use the cell's background as expectation
      deviations <- deviations / motif_bg_exp
    }else if(denominator=="global"){
      # global motif expectation (original CV)
      # denom = motif peak counts scaled by the cell's libsize
      sf <- depth / sum(background@expectation)
      deviations <- deviations/outer(st$globalAvg, sf)
    }
    a$deviations <- deviations
  }
  a$total <- rowSums(observed_motif_sum)
  a
}
