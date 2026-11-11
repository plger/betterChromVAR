#' computeBackgrounds
#'
#' @param object A SummarizedExperiment (or SingleCellExperiment) with an assay
#'    'counts', or a (sparse) matrix of counts.
#' @param bins A `bcvBackground` object, as produced by 
#'   \code{\link{getBackgroundBins}}.
#' @param grouping An optional factor or vector coercible to a factor indicating
#'   the groupings of the columns of `object`. This is optionally used to 
#'   compute the base expectation such that rare cell types are given as much 
#'   weight as abundant ones. In single-cell data, the grouping can for 
#'   instance be the interaction of samples and cell types. (The name of a 
#'   colData column of `object` can also be provided.)
#' @param expectation Optional vector of length equal to `nrow(object)` 
#'   giving the expected counts. If NULL, defaults to mean counts (eventually
#'   grouped, see `grouping` and \code{\link{getExpectation}}).
#' @param verbose Logical; whether to output progress messages.
#'
#' @returns A `bcvBackground` object with bins*samples slots filled, for use
#'   with \code{\link{computeDeviationsAnalytic}}.
#' @export
#'
#' @examples
#' attach(getDummyData())
#' # if GC content not already in the object, use:
#' # counts <- addGCBias(counts, genome=YOUR_GENOME)
#' 
#' # we fist get the background bins:
#' bg <- getBackgroundBins(counts)
#' # then we can compute the backgrounds for each sample:
#' bg <- computeBackgrounds(counts, bg)
#' # for use in computeDeviationsAnalytic...
computeBackgrounds <- function(object, bins, grouping=NULL, expectation=NULL,
                               verbose=FALSE){
  # Check input validity
  stopifnot(is.null(expectation) || length(expectation)==nrow(object))
  depth <- NULL
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
  stopifnot(length(dim(counts))==2)
  if(!is(counts, "matrix") && !inherits(counts, "Matrix"))
    stop("`object` should be a SummarizedExperiment or SingleCellExperiment,",
         " or a (sparse) matrix of counts.")
  
  stopifnot(is(bins, "bcvBackground"))
  if(length(bins@peak2bin)!=nrow(object))
    stop("The number of peaks in the `bins` object does not match the rows",
         "of `object`.")
  
  if(length(bins@depth)>0) depth <- bins@depth
  if(is.null(depth)) depth <- colSums(counts)

  grouping <- .groupingInput(grouping, object)
  if(is.null(expectation)){
    if(length(bins@expectation)==nrow(counts)){
      if(verbose) message("Using pre-computed expectation")
      expectation <- bins@expectation
    }else{
      expectation <- getExpectation(counts, grouping)
    }
  }
  if(any(expectation==0)){
    stop("Some peaks have an expectation of zero, most likely because they ",
         "have zero counts. Please remove them.")
  }
  
  bin_map <- bins@peak2bin
  binBinProbs <- bins@binBinProbs
  
  # sparse mapping from peaks to bins
  bin2peakMat <- sparseMatrix(i=bin_map, j=seq_along(expectation), 
                              dims=c(nrow(binBinProbs), length(expectation)))
  
  binCounts <- bin2peakMat %*% counts
  
  # mean and variance of the counts of each bin's background peaks (B x S)
  E <- binBinProbs %*% binCounts
  V <- binBinProbs %*% (bin2peakMat %*% (counts^2))
  if(is(V,"CsparseMatrix") && identical(V@p, E@p) && identical(V@i, E@i)){
    # same sparsity pattern: subtract on the values to save memory
    V@x <- V@x - E@x^2
  }else{
    V <- V - E^2
  }
  V@x[which(V@x<0)] <- 0
  V <- drop0(V)
  
  if(is(E,"dgeMatrix")) E <- as.matrix(E)
  if(is(V,"dgeMatrix")) V <- as.matrix(V)
  
  bins@E <- E
  bins@V <- V
  bins@expectation <- expectation
  bins@depth <- depth
  validObject(bins)
  
  bins
}
