#' computeDeviationsWeighted
#' 
#' A variant of \code{\link{computeDeviationsAnalytic}} enabling the 
#' computation of deviations from weighted foreground counts. Specifically, this
#' functions handles the normalization of the difference in magnitude between
#' the (weighted) foreground and background.
#'
#' @param weightedMotifCounts A matrix of weighted counts per motif (rows) and 
#'   sample (columns), or a `SummarizedExperiment` containing this as first 
#'   assay.
#' @param unweightedPeakCounts A matrix of unweighted counts per peak (rows) and 
#'   sample (columns), or a `SummarizedExperiment` containing this as first 
#'   assay.
#' @param annotations Peak annotation (sparse) matrix, with motifs as columns,
#'    or a SummarizedExperiment containing this in the first assay. Values 
#'    should be either logical or between 0 and 1.
#' @param bg Either a \code{\link[betterChromVAR]{bcvBackground-class}} object
#'   as produced by \code{\link{computeBackgrounds}}, or a 
#'   `SummarizedExperiment` of background peak counts (with bias data in 
#'   `rowData`). If missing, will be created based on `unweightedPeakCounts`.
#' @param retSE Logical; whether to return a SummarizedExperiment object.
#' @param fpc Logical; whether to apply the finite-population correction to the
#'   null variance used for the z-scores (see
#'   \code{\link{computeDeviationsAnalytic}}).
#' @param ... Passed to \link{getBackgroundBins} (can for instance be used to
#'   pass bias info if not contained in the objects). Ignored if `bg` is a
#'   \code{\link[betterChromVAR]{bcvBackground-class}} object.
#'
#' @details
#' The background expectation of the (unweighted) peak counts is rescaled, for
#' each motif, by the ratio between the mean weighted and mean unweighted 
#' counts of the motif, and its variance by the square of that ratio. This 
#' assumes that the weights are similar across the peaks of a motif: if the 
#' weights vary a lot between peaks (coefficient of variation `CV`), the 
#' variance of the weighted counts is larger, by a factor ~`1+CV^2`, than the 
#' variance used here, and the z-scores are overestimated accordingly. The 
#' finite-population correction (`fpc`) also assumes that the weighted counts 
#' are proportional to the unweighted counts of the same peaks, from which the
#' expectation is estimated. Neither the z-scores nor the `deviations` depend 
#' on the overall magnitude of the weights.
#'
#' @returns A SummarizedExperiment (or a list if `retSE=FALSE`). If overall
#'   motif variability and their significance are additionally needed, see 
#'   \code{\link{computeMotifVariability}}.
#' @seealso [computeDeviationsAnalytic()]
#' @export
#' @examples
#' attach(getDummyData())
#' # if GC content not already in the object, use:
#' # counts <- addGCBias(counts, genome=YOUR_GENOME)
#' # For the purpose of this example, we'll use standard (unweighted counts),
#' # although at this step we'd compute counts weighted in the desired fashion:
#' motifCounts <- Matrix::t(motifMatches) %*% assay(counts)
#' dev1 <- computeDeviationsWeighted(motifCounts, counts, motifMatches)
#' dev1
#' # in this case, the results are essentially identical to :
#' dev2 <- betterChromVAR(counts, motifMatches)
#' sapply(assayNames(dev1), function(a){
#'   cor(as.vector(assay(dev1, a)), as.vector(assay(dev2, a)))
#' })
computeDeviationsWeighted <- function(weightedMotifCounts, unweightedPeakCounts,
                                      annotations, bg=NULL, retSE=TRUE, fpc=TRUE, ...){
  stopifnot(ncol(weightedMotifCounts)==ncol(unweightedPeakCounts))
  stopifnot(nrow(weightedMotifCounts)==ncol(annotations))
  stopifnot(nrow(unweightedPeakCounts)==nrow(annotations))
  if(is.null(bg)){
    bg <- computeBackgrounds(unweightedPeakCounts,
                             getBackgroundBins(unweightedPeakCounts, ...))
  }else if(.isSElike(bg)){
    stopifnot(nrow(bg)==nrow(annotations) && 
                ncol(bg)==ncol(weightedMotifCounts))
    bg <- computeBackgrounds(bg, getBackgroundBins(bg, ...))
  }else{
    if(length(bg@depth)==0 || length(bg@expectation)==0)
      stop("Incomplete background; please run computeBackgrounds() first.")
    if( (length(bg@depth)!=ncol(weightedMotifCounts)) ||
        (length(bg@expectation)!=nrow(unweightedPeakCounts)) )
      stop("The `background` object does not match the other objects' dimensions.")
  }
  
  CD <- motifCD <- NULL
  if( inherits(annotations, "SummarizedExperiment") ){
    motifCD <- colData(annotations)
    annotations <- assay(annotations)
  }
  .checkAnnotations(annotations)
  
  if(.isSElike(weightedMotifCounts)){
    CD <- colData(weightedMotifCounts)
    weightedMotifCounts <- assay(weightedMotifCounts)
  }
  if(.isSElike(unweightedPeakCounts)){
    if(is.null(CD)) CD <- colData(unweightedPeakCounts)
    unweightedPeakCounts <- assay(unweightedPeakCounts)
  }
  
  binMap <- bg@peak2bin
  bin2peakMat <- sparseMatrix(i = binMap, j = seq_along(binMap), 
                              dims = c(nrow(bg@binBinProbs), length(binMap)))
  motifBinCounts <- Matrix::t(annotations) %*% Matrix::t(bin2peakMat)
  motif_bg_exp_unweighted <- as.matrix(motifBinCounts %*% bg@E)
  unwMoCounts <- .annotCrossprod(annotations, unweightedPeakCounts)
  unwMoCounts <- Matrix::rowMeans(.fastColNorm(unwMoCounts))*sum(unwMoCounts)
  fg <- sum(weightedMotifCounts)*Matrix::rowMeans(.fastColNorm(weightedMotifCounts))
  motif_sf <- fg / unwMoCounts
  motif_bg_exp <- motif_bg_exp_unweighted * motif_sf
  
  numerator <- weightedMotifCounts - motif_bg_exp
  st <- .fpcStats(annotations, bg)
  
  if(isTRUE(fpc)){
    bg_variance_scaled <- .fpcVariance(annotations, bg, bg@V, st=st) * (motif_sf^2)
  }else{
    bg_variance_scaled <- as.matrix(motifBinCounts %*% bg@V) * (motif_sf^2)
  }
  z <- numerator / sqrt(pmax(0, bg_variance_scaled))
  globalMotifAvg <- as.vector(Matrix::crossprod(annotations, 
                                                bg@expectation))
  sf <- bg@depth/sum(bg@expectation)
  # the global expectation is rescaled to the scale of the weighted counts, as
  # the numerator, so that the deviations are relative to it
  deviations <- numerator/outer(globalMotifAvg*motif_sf, sf)
  a <- lapply(list(deviations=deviations, z=z), \(x){
    if(is(x, "dgeMatrix")) x <- as.matrix(x)
    x
  })
  
  if(!isTRUE(retSE)) return(a)
    
  d <- data.frame(N=Matrix::colSums(annotations),
                  Neff=.effectiveN(annotations, bg, bg@expectation, fpc, st=st))
  
  .packageDevSE(a, SummarizedExperiment(list(), colData=CD),motifCD, d)
}
