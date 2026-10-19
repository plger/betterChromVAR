#' getBackgroundBins
#' 
#' Computes chromVAR-like background (i.e. bias) bins, as well as bin-to-bin 
#' selection probabilities needed for \code{\link{betterChromVAR}}.
#' 
#' @param x A SummarizedExperiment containing a 'counts' assay, or a matrix of
#'   counts, or a vector of expected (e.g. mean) counts.
#' @param bias A vector of length equal to `ncol(object)` specifying the 
#'   per-peak bias (i.e. GC content). If omitted, will try to get it
#'   from `rowData(x)$bias`.
#' @param flbias A vector of length equal to `ncol(object)` specifying the 
#'   per-peak fragment length bias (e.g. log10-transformed median length of 
#'   fragments overlapping each region). If omitted, will try to get it
#'   from `rowData(object)$flbias`. If absent, will use 2-dimensional bias bins
#'   (which is already pretty good).
#' @param w Standard deviation of the Gaussian kernel.
#' @param bs Number of bins per dimension. This can be a single integer (total
#'   bins = `bs^2`), or an integer vector of length 2 (if `flbias=NULL`) or 3 
#'   (in which case there are `prod(bs)` total bins). The values specify the 
#'   number of bins for, in order: enrichment, GC and fragment length. By 
#'   default, `bs=50` if `flbias` is not provided (mimicking chromVAR), and
#'   `bs=c(30, 30, 6)` if it is.
#' @param maxShare Optional maximum share that a peak can have in the 
#'   expectation of its own background. If not NULL, the width of the kernel is
#'   increased (up to `maxWidth` times `w`) for the bins that would otherwise 
#'   exceed it, so that sparsely-populated bins draw more from their neighbors.
#' @param maxWidthFactor The maximum factor by which `w` can be multiplied when 
#'   `maxShare` is used.
#' @param pseudo Optional pseudocount to be added. This should not be needed 
#'   with standard workflows.
#' @param verbose Whether to print processing info.
#' 
#' @details
#' The procedure underlying this function is the same as in 
#' `chromVAR::getBackgroundPeaks`, with the following differences:
#' * Rather than producing a set of background peaks for each input peak, the 
#'   function returns peak-to-bin mappings and bin-to-bin background selection
#'   probabilities, which enables an analytic background computation. It is, as
#'   such, entirely deterministic.
#' * The function supports the optional use of a third bias dimension, provided
#'   through the `flbias` argument, meant for fragment length bias. This is 
#'   still an experimental feature.
#' 
#' @return A `bcvBackground` object, to be used with 
#'   \code{\link{computeBackgrounds}}.
#' @references
#'   Schep A.N., Wu B., Buenrostro J.D., Greenleaf W.J. (2017) chromVAR: 
#'   inferring transcription-factor-associated accessibility from 
#'   single-cell epigenomic data, Nature Methods, doi: 10.1038/nmeth.4401
#'   
#' @importFrom stats cov dnorm dist
#' @importFrom SummarizedExperiment assay rowData
#' @importFrom matrixStats colMins colMaxs
#' @importFrom methods as is new
#' @importFrom Matrix drop0
#' @export
#' @examples
#' counts_se <- getDummyData()$counts
#' background <- getBackgroundBins(counts_se)
getBackgroundBins <- function(x, bias=NULL, flbias=NULL, w=0.1, bs=NULL,
                              maxShare=0.1, maxWidthFactor=4, pseudo=0, 
                              verbose=TRUE){
  if(.isSElike(x)){
    if(is.null(bias)) bias <- rowData(x)$bias
    if(is.null(flbias)) flbias <- rowData(x)$flbias
    x <- Matrix::rowMeans(assay(x, "counts"))
  }else if(is.matrix(x) || is(x, "Matrix")){
    x <- Matrix::rowMeans(x)
  }
  if(is.null(bias)) stop("`bias` not provided, and not found in the object.")
  stopifnot(length(bias)==nrow(x))
  
  stopifnot("Some peaks have zero counts; please remove them."=sum(x==0)==0)
  
  if(!is.null(flbias)){
    if(is.null(bs)) bs <- c(30L, 30L, 6L)
    stopifnot(length(bs)==3)
    stopifnot(length(flbias)==length(bias))
  }else{
    if(is.null(bs)) bs <- 50
    if(length(bs)==1) bs <- c(bs,bs)
    stopifnot(length(bs)==2)
  }
  bs <- as.integer(bs)
  stopifnot(all(bs>=1))
  if(verbose) 
    message("Creating ", paste(bs,collapse="*"),"=",prod(bs)," bias bins ",
            "and computing their sampling distances")
  
  # Mahalanobis transformation
  norm_mat <- cbind(log10(x+pseudo), bias)
  if(!is.null(flbias)) norm_mat <- cbind(norm_mat, flbias)
  cv <- cov(norm_mat)
  diag(cv) <- diag(cv) + pseudo/1000
  transMat <- t(forwardsolve(t(chol(cv)), t(norm_mat)))
  
  # Calculate min/max for the bin boundaries
  minCoords <- colMins(transMat)
  maxCoords <- colMaxs(transMat)
  range_coords <- maxCoords - minCoords
  
  # Map peak coordinates to bins
  idx1 <- round((transMat[, 1] - minCoords[1]) / 
                  range_coords[1] * (bs[1] - 1)) + 1
  idx2 <- round((transMat[, 2] - minCoords[2]) /
                  range_coords[2] * (bs[2] - 1)) + 1
  # Ensure indices fall within [1, bs]
  idx1 <- pmax(1, pmin(bs[1], idx1))
  idx2 <- pmax(1, pmin(bs[2], idx2))
  
  if(!is.null(flbias)){
    idx3 <- round((transMat[, 3] - minCoords[3]) /
                    range_coords[3] * (bs[3] - 1)) + 1
    idx3 <- pmax(1, pmin(bs[3], idx3))
  }
  
  # Linearize index
  peak2bin <- idx1 + (idx2 - 1) * bs[1]
  if(!is.null(flbias)){
    peak2bin <- peak2bin  + (idx3 - 1) * (bs[1] * bs[2])
  }
  
  binDensity <- tabulate(peak2bin, nbins = prod(bs))
  
  # bin center grid (for distance calculation)
  grid_args <- lapply(seq_along(minCoords), function(i) {
    seq(minCoords[i], maxCoords[i], length.out = bs[i])
  })
  bin_data <- do.call(expand.grid, grid_args)

  # bin-to-bin probability matrix
  D <- as.matrix(dist(bin_data))
  W <- dnorm(D, 0, w)
  
  if(!is.null(maxShare)){
    stopifnot(length(maxShare)==1 && maxShare>0 && maxShare<1)
    stopifnot(length(maxWidthFactor)==1 && maxWidthFactor>=1)
    W <- .adaptiveKernel(W, D, w, maxShare, maxWidthFactor, x, peak2bin)
  }
  
  normalizer <- as.vector(W %*% binDensity)
  # Avoid division by zero for empty regions
  normalizer[normalizer < 1e-9] <- 1
  binBinProbs <- W/normalizer
  
  binBinProbs[binBinProbs < 1e-6] <- 0
  binBinProbs <- as(binBinProbs, "sparseMatrix")
    
  new("bcvBackground",
      dims = bs,
      peak2bin = as.integer(peak2bin), 
      binDensity = binDensity, 
      binBinProbs = binBinProbs,
      E = NULL,
      V = NULL,
      expectation = numeric(), 
      depth = integer())
}


#' sampleBackgroundPeaks
#' 
#' Given a background generated by \code{\link{getBackgroundBins}}, samples 
#' background peaks for each input peak.
#'
#' @param background A `bcvBackground` produced by 
#'   \code{\link{getBackgroundBins}}.
#' @param niterations Number of background peaks to sample for each target peak. 
#'
#' @return A peaks x niterations matrix of integers representing the indices of
#'   the sampled background peaks.
#'   
#' @details
#' This function is not used by \code{\link{betterChromVAR}}, which is 
#' deterministic, but for other applications requiring an outputs similar to 
#' that of the original \code{\link[chromVAR]{getBackgroundPeaks}}.
#' 
#' @export
#'
#' @examples
#' counts_se <- getDummyData()$counts
#' background <- getBackgroundBins(counts_se)
#' bg_peaks <- sampleBackgroundPeaks(background, niterations=20)
sampleBackgroundPeaks <- function(background, niterations=50){
  stopifnot(is(background, "bcvBackground"))
  peak2bin <- background@peak2bin
  binBinProbs <- background@binBinProbs
  density <- background@binDensity
  n_bins <- nrow(binBinProbs)

  out <- matrix(0L, nrow=length(peak2bin), ncol=niterations)
  
  peaks_in_bin <- split(seq_along(peak2bin), peak2bin)
  
  for (b in names(peaks_in_bin)) {
    target_peaks <- peaks_in_bin[[b]]
    b <- as.integer(b)
    n_targets <- length(target_peaks)
    n_samples <- n_targets * niterations
    
    # Prob of selecting candidate bin `j` for target bin `b` is proportional to:
    # (Probability of selecting a specific peak in `j`) * (Total peaks in `j`)
    b_probs <- as.numeric(binBinProbs[b, ]) * density
    if(sum(b_probs) == 0) b_probs[b] <- 1
    
    # sample bins
    sampled_bins <- sample(n_bins, size=n_samples, replace=TRUE, prob=b_probs)
    
    # convert sampled bins to sampled peaks
    sampled_counts <- tabulate(sampled_bins, nbins = n_bins)
    active_candidate_bins <- which(sampled_counts > 0)
    
    sampled_peaks <- integer(n_samples)
    for (j in active_candidate_bins) {
      cand_peaks <- peaks_in_bin[[as.character(j)]]
      # handle the single-peak bin edge case to avoid sampling from 1:x.
      if(length(cand_peaks) == 1){
        sampled_peaks[sampled_bins == j] <- cand_peaks
      }else{
        sampled_peaks[sampled_bins == j] <- 
          sample(cand_peaks, size=sampled_counts[j], replace=TRUE)
      }
    }
    
    # assign to the target rows in the output matrix
    out[target_peaks, ] <- matrix(sampled_peaks, ncol=niterations)
  }
  
  return(out)
}




# Widens the gaussian kernel of the bins where a peak would otherwise make up
# more than `maxShare` of the expectation of its own background.
.adaptiveKernel <- function(W, D, w, maxShare, maxWidth, x, peak2bin){
  B <- nrow(W)
  fb <- factor(peak2bin, levels=seq_len(B))
  ebin <- as.numeric(tapply(x, fb, sum))
  emax <- as.numeric(tapply(x, fb, max))
  ebin[is.na(ebin)] <- 0
  emax[is.na(emax)] <- 0
  share <- (diag(W) * emax) / as.vector(W %*% ebin)
  todo <- which(emax>0 & share > maxShare)
  for(u in todo){
    m <- 1
    Wu <- W[u,]
    while(Wu[u]*emax[u]/sum(Wu*ebin) > maxShare && m < maxWidth){
      m <- min(m*1.25, maxWidth)
      Wu <- dnorm(D[u,], 0, w*m)
    }
    W[u,] <- Wu
  }
  W
}
