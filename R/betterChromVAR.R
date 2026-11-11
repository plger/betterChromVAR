#' betterChromVAR
#' 
#' A fast, analytic implementation of `chromVAR`.
#' This is a wrapper around the \code{\link{getBackgroundBins}},
#' \code{\link{computeBackgrounds}}, and \code{\link{computeDeviationsAnalytic}}
#' steps. It additionally allows for multithreading. For more control or 
#' optimization, see the individual steps.
#' 
#' @param object A SummarizedExperiment (or SingleCellExperiment) with an assay
#'    'counts', and with a 'bias' column in `rowData(object)`. Note that the 
#'    regions should have similar widths.
#' @param annotations Peak annotation (sparse) matrix, with motifs as columns,
#'    or a SummarizedExperiment containing this in the first assay. Values 
#'    should be either logical or between 0 and 1.
#' @param grouping An optional factor or vector coercible to a factor indicating
#'   the groupings of the columns of `object`. This is optionally used to 
#'   compute the base expectation such that rare cell types are given as much 
#'   weight as abundant ones. In single-cell data, the grouping can for 
#'   instance be the interaction of samples and cell types. This should either 
#'   be a vector coercible to factor of length equal to `ncol(object)`, or a 
#'   character of length 1 specifying a column of `colData(object)`.
#' @param nthreads Either an integer scalar indicating the number of threads to
#'   use, or a `BiocParallelParam` object.
#' @param verbose Logical; whether to output progress messages (default FALSE).
#' @param fpc Logical; whether to apply the finite-population correction to the
#'   null variance used for the z-scores (default TRUE, see
#'   \code{\link{computeDeviationsAnalytic}}).
#' @param ridge An optional positive number, enabling the decorrelation of 
#'   motifs that share peaks (smaller values decorrelate more strongly; values
#'   of 1 to 3 are recommended). By default (NULL), motifs are not 
#'   decorrelated. See \code{\link{computeDeviationsAnalytic}} for details.
#' @param ... Passed to \code{\link{getBackgroundBins}}.
#' @author Pierre-Luc Germain
#' 
#' @details
#' Contrarily to the original chromVAR, this function is entirely deterministic,
#' and achieves higher precision and much higher efficiency through two changes:
#' 1) working with expected background sampling mean and variances, rather than
#' actual permutations, and 2) computing expectations and variance at the level
#' of bias bins, instead of in the peak-space. The function can additionally 
#' handles non-binary annotations (e.g. probability scores), a third 
#' bias dimension (fragment length bias, which should be stored in 
#' `rowData(object)$flbias` see \code{\link{getBackgroundBins}} for details),
#' an adaptive kernel width, and ridge-based decorrelation of the motifs.
#' 
#' @references
#'   Schep A.N., Wu B., Buenrostro J.D., Greenleaf W.J. (2017) chromVAR: 
#'   inferring transcription-factor-associated accessibility from 
#'   single-cell epigenomic data, Nature Methods, doi: 10.1038/nmeth.4401
#' 
#' @return A SummarizedExperiment containing the adjusted deviations and 
#'   z-scores for each motif/sample. The rowData additionally contains the 
#'   number of motif matches and their variability.
#' @importFrom SummarizedExperiment SummarizedExperiment assay colData rowData
#' @importFrom S4Vectors metadata
#' @importFrom Matrix crossprod sparseMatrix kronecker Diagonal cbind2 colSums
#' @importFrom BiocParallel bplapply SerialParam MulticoreParam bpnworkers
#' @importFrom stats p.adjust pchisq
#' @importFrom matrixStats rowSds
#' @export
#' @examples
#' attach(getDummyData())
#' # if GC content not already in the object, use:
#' # counts <- addGCBias(counts, genome=YOUR_GENOME)
#' dev <- betterChromVAR(counts, motifMatches)
#' dev
#' # note that this is the exact equivalent of doing:
#' # bg <- getBackgroundBins(counts)
#' # bg <- computeBackgrounds(counts, bg)
#' # dev <- computeDeviationsAnalytic(counts, bg, motifMatches)
betterChromVAR <- function(object, annotations, grouping=NULL, nthreads=NULL,
                           verbose=FALSE, fpc=TRUE, ridge=NULL, ...){
  
  stopifnot(.isSElike(object))
  stopifnot(nrow(object) == nrow(annotations))
  stopifnot(!is.null(rowData(object)$bias))
  bias <- rowData(object)$bias
  motifCD <- flbias <- NULL
  if(!is.null(rowData(object)$flbias)) flbias <- rowData(object)$flbias
  
  if(.isSElike(annotations)){
    motifCD <- colData(annotations)
    annotations <- assay(annotations)
  } 
  .checkAnnotations(annotations)
  
  grouping <- .groupingInput(grouping, object, TRUE)
  
  if(verbose) message("Preparing bias bins")
  expectation <- getExpectation(object, grouping)
  bg <- getBackgroundBins(expectation, bias=bias, flbias=flbias, 
                          verbose=verbose, ...)
  
  # what depends only on the annotations is computed once, rather than for
  # each chunk of cells
  st <- .annotationStats(annotations, bg, expectation)
  if(!is.null(ridge)){
    if(verbose) message("Computing the overlap between motifs")
    ridge <- .ridgeOperator(annotations, bg, ridge, expectation=expectation,
                            st=st)
  }

  ngroups <- length(levels(grouping))
                    
  BPPARAM <- .bpFromThreads(nthreads, progressbar=(verbose && ngroups>1))

  i <- seq_len(ncol(object))
  counts <- assay(object, "counts")
  depth <- NULL
  if(is.numeric(object$depth)) depth <- object$depth
  if((nW <- BiocParallel::bpnworkers(BPPARAM))>1 && ncol(counts) > 100){
    if(verbose) message("Computing backgrounds and deviations")
    chunks <- split(i, cut(i, nW, labels=FALSE))
    res <- bplapply(chunks, BPPARAM=BPPARAM, function(i){
      bg <- computeBackgrounds(counts[,i], bg[,i], expectation=expectation,
                               verbose=FALSE)
      .computeDeviations(counts[,i], bg, annotations, st=st, depth=depth[i],
                         compute=c("deviations","z"), fpc=fpc, ridge=ridge)
    })
    res <- list(deviations=Reduce(cbind2,
                                  lapply(res, function(x) x$deviations)),
                z=Reduce(cbind2, lapply(res, function(x) x$z)),
                total=Reduce("+", lapply(res, function(x) x$total)))
  }else{
    if(verbose) message("Computing backgrounds")
    bg <- computeBackgrounds(counts, bg, expectation=expectation,
                             verbose=verbose)
    if(verbose) message("Computing deviations")
    res <- .computeDeviations(counts, bg, annotations, st=st, depth=depth,
                              compute=c("deviations","z"), fpc=fpc,
                              ridge=ridge, verbose=verbose)
  }
  res <- res[c("deviations","z","total")]

  d <- cbind(data.frame(N=colSums(annotations),
                        Neff=.effectiveN(annotations, bg, expectation, fpc, st=st),
                        total=res$total),
             computeMotifVariability(res$z))

  .packageDevSE(res[1:2], object, motifCD, d)
}
