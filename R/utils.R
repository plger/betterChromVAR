.fastColNorm <- function(x, cs=Matrix::colSums(x)){
  x %*% Diagonal(x = 1/cs)
}

.fastColAgg <- function(x, by){
  by <- factor(by)
  stopifnot(length(by)==ncol(x))
  mapMat <- sparseMatrix(
    i = as.integer(by),
    j = seq_along(by),
    dims = c(length(levels(by)), ncol(x)),
    dimnames = list(levels(by), colnames(x))
  )
  tcrossprod(x, mapMat)
}

#' getExpectation
#' 
#' Computes expected counts (a glorified rowMeans)
#'
#' @param counts A count matrix, or object inheriting SummarizedExperiment with
#'   a 'counts' assay.
#' @param grouping An optional vector of length equal to `ncol(counts)`,
#'   indicating the grouping of the cells. If provided, cells will be averaged
#'   by group before averaging across groups.
#' @param normalize Logical; whether to normalize data between averaging (but
#'   after grouping). Default TRUE and highly recommended if providing 
#'   `grouping`.
#'
#' @returns A vector of expectation for each row of `counts`
#' @export
#'
#' @importFrom DelayedMatrixStats rowMeans2
#' @examples
#' attach(getDummyData())
#' e <- getExpectation(counts)
getExpectation <- function(counts, grouping=NULL, normalize=TRUE){
  if(.isSElike(counts)) counts <- assay(counts, "counts")
  if(is.null(grouping) || length(unique(grouping))==1)
    return(DelayedMatrixStats::rowMeans2(counts))
  grouping <- factor(grouping)
  stopifnot(length(grouping)==ncol(counts))
  # compute expectation based on an average of group averages
  agcnt <- .fastColAgg(counts, grouping)
  if(is(agcnt, "dgeMatrix")) agcnt <- as.matrix(agcnt)
  if(normalize){
    cs <- Matrix::colSums(agcnt)
    agcnt <- .fastColNorm(agcnt, cs=cs)*median(cs)
    if(is(agcnt, "dgeMatrix")) agcnt <- as.matrix(agcnt)
  }
  DelayedMatrixStats::rowMeans2(agcnt)
}


#' Dummy data for testing purposes
#'
#' @param nRegions Number of regions to generate
#' @param nSamples Number of samples to generate
#' @param nMotifs Number of motifs to generate
#'
#' @returns A list with the slots `counts` (a peak counts SummarizedExperiment)
#'   and `matches` (a sparse matrix of binary motif matches per peaks)
#' @importFrom Matrix Matrix
#' @importFrom stats rnorm rnbinom
#' @importFrom GenomicRanges GRanges
#' @importFrom IRanges IRanges
#' @importFrom SummarizedExperiment rowData<-
#' @export
#'
#' @examples
#' out <- getDummyData()
#' (counts <- out$counts)
#' matches <- out$motifMatches
getDummyData <- function(nRegions=500, nSamples=10, nMotifs=5){
  mu <- sample.int(100, nRegions, replace=TRUE)
  counts <- matrix(rnbinom(nRegions * nSamples, mu=rep(mu,nSamples), size=2),
                   nrow=nRegions, ncol=nSamples)
  counts <- counts[which(rowSums(counts)>0),]
  nRegions <- nrow(counts)
  colnames(counts) <- paste0("sample", seq_len(nSamples))
  gr <- GRanges("chr1", IRanges(seq_len(nRegions)*100, width=20))
  counts <- SummarizedExperiment(list(counts=counts), rowRanges=gr)
  rowData(counts)$bias <- pmin(1,pmax(0,rnorm(nRegions, mean=0.5, sd=0.05)))
  matches <- Matrix(
    data=sample(c(0L, 1L), nRegions*nMotifs, replace=TRUE, prob=c(0.85, 0.15)), 
    nrow=nRegions, ncol=nMotifs, sparse=TRUE)
  colnames(matches) <- paste0("motif",seq_len(nMotifs))
  list(counts=counts, motifMatches=matches)
}

#' addGCBias
#' 
#' Add the `bias` column to the object's rowData, containing the regions' 
#' proportion of Gs and Cs.
#'
#' @param object An object inheriting RangedSummarizedExperiment or GRanges.
#' @param genome A BSgenome object or any other genome object supported by 
#'   \code{\link[Biostrings]{getSeq}}.
#'
#' @returns `object` with the GC content in `mcols(object)$bias` (if GRanges) 
#'   or `rowData(object)$bias`.
#' @importFrom Biostrings getSeq letterFrequency
#' @importFrom SummarizedExperiment rowRanges rowRanges<- mcols mcols<-
#' @export
#'
#' @examples
#' # not run:
#' # se <- addGCBias(se, genome)
addGCBias <- function(object, genome){
  if(.isSElike(object)){
    stopifnot(!is.null(rowRanges(object)))
    rowRanges(object) <- addGCBias(rowRanges(object), genome)
    return(object)
  }
  seqs <- Biostrings::getSeq(x=genome, object)
  # same as chromVAR:
  nucfreqs <- letterFrequency(seqs, c("A", "C", "G", "T"))
  gc <- rowSums(nucfreqs[, 2:3]) / rowSums(nucfreqs)
  mcols(object)$bias <- gc
  object
}

# Finite-population correction for the null variance of a motif's summed counts.
#
# For a motif m, the deviation of a cell is S = sum_q a_qm (x_q - E_b(q)), where 
# x_q is the count of peak q, a_qm its annotation (0/1 or a weight), and E_b the
# expected count per peak of its bin b, which is estimated from the cell's own 
# counts in the bins around b (E_b = sum_b' P_bb' T_b', with T_b' the counts in
# bin b' and P the bin-to-bin probabilities). The peaks of the motif are thus
# part of their own expectation. Writing S as sum_q c_q x_q, with
#   c_q = a_qm - g_b(q),   g_b' = sum_b N_b P_bb'   (N_b: sum of a_qm in bin b),
# and assuming independent peaks with variance V_b, the variance of S is
#   sum_q c_q^2 V_b(q) = sum_b V_b * coef_b,   
#   coef_b = A2_b - 2 g_b N_b + n_b g_b^2   (A2_b: sum of a_qm^2 in bin b; n_b: 
#   number of peaks in bin b)
# instead of sum_b V_b N_b. `.fpcStats` returns these quantities (bins x motifs).
.fpcStats <- function(annotations, background){
  nb <- background@binDensity
  bin2peak <- sparseMatrix(i=background@peak2bin, j=seq_along(background@peak2bin),
                           dims=c(length(nb), nrow(annotations)))
  N <- as.matrix(bin2peak %*% annotations)
  A2 <- as.matrix(bin2peak %*% (annotations * annotations))
  g <- as.matrix(Matrix::crossprod(background@binBinProbs, N))
  list(N=N, A2=A2, g=g, coef=A2 - 2*g*N + nb*g^2)
}

# Null variance of the motif sums (motifs x cells) with the finite-population 
# correction; V is the bins x cells matrix of background variances.
.fpcVariance <- function(annotations, background, V, st=NULL){
  if(is.null(st)) st <- .fpcStats(annotations, background)
  as.matrix(Matrix::crossprod(st$coef, V))
}

# Quantities that depend only on the annotations and the bins (and not on the
# cells): the output of `.fpcStats` and the global expectation of each motif.
# They are computed once and passed on (as `st`), rather than for each chunk
# of cells.
.annotationStats <- function(annotations, background, 
                             expectation=background@expectation){
  st <- .fpcStats(annotations, background)
  st$globalAvg <- as.vector(Matrix::crossprod(annotations, expectation))
  st
}

# Whether annotations are sparse enough for sparse algebra to be worth it
.isSparseAnnotation <- function(x, maxDensity=0.2){
  inherits(x, "sparseMatrix") && Matrix::nnzero(x)/prod(dim(x)) < maxDensity
}

# Splits the rows of x in blocks which, once dense, have at most `blockCells`
# values
.rowBlocks <- function(x, blockCells=2.5e7){
  i <- seq_len(nrow(x))
  unname(split(i, ceiling(i/max(1L, floor(blockCells/ncol(x))))))
}

# Motif sums, i.e. t(annotations) %*% counts, as a dense matrix (motifs x 
# cells). Sparse annotations (e.g. binary matches) use a sparse cross-product.
# For annotations that are not sparse (e.g. graded matches), this is slow and
# memory-hungry, and the annotations are instead processed as dense blocks of 
# rows.
.annotCrossprod <- function(annotations, counts, maxDensity=0.2, 
                            blockCells=2.5e7){
  if(!inherits(annotations, "sparseMatrix") || 
     .isSparseAnnotation(annotations, maxDensity))
    return(as.matrix(Matrix::crossprod(annotations, counts)))
  annotations <- .asDgC(annotations)
  res <- Reduce("+", lapply(.rowBlocks(annotations, blockCells), function(i){
    as.matrix(Matrix::crossprod(as.matrix(annotations[i,,drop=FALSE]),
                                counts[i,,drop=FALSE]))
  }))
  dimnames(res) <- list(colnames(annotations), colnames(counts))
  res
}

# Effective number of peaks of each motif: the number of peaks (N) of a motif 
# with the same sensitivity of its z-score to a same relative change in all of 
# its peaks. A change in the counts of the motif's peaks also changes the 
# expectation (see above), which attenuates the signal as well as the noise; 
# the z-score of the motif is proportional to the square root of
#   Neff = N * signal^2 / (variance * X),  
# with signal = sum_b E_b (A2_b - g_b N_b), variance = sum_b E_b coef_b (or 
# sum_b E_b N_b without finite-population correction) and X = sum_b E_b N_b, 
# where E_b is the expected count per peak around bin b (the variance of the 
# counts being taken as proportional to their mean). Without overlap between 
# a motif's peaks and their expectation Neff = N; with an even spread of the
# motif's peaks across bins, Neff = N(1-N/P).
.effectiveN <- function(annotations, background, expectation, fpc=TRUE,
                        st=NULL){
  if(is.null(st)) st <- .fpcStats(annotations, background)
  bin2peak <- sparseMatrix(i=background@peak2bin, j=seq_along(background@peak2bin),
                           dims=c(length(background@binDensity), nrow(annotations)))
  Eb <- as.numeric(background@binBinProbs %*% (bin2peak %*% expectation))
  X <- colSums(Eb * st$N)
  signal <- colSums(Eb * (st$A2 - st$g*st$N))
  variance <- colSums(Eb * (if(isTRUE(fpc)) st$coef else st$N))
  colSums(st$N) * signal^2 / (variance * X)
}

.groupingInput <- function(grouping, object, fillNULL=FALSE){
  if(is.null(grouping)){
    if(fillNULL){
      grouping <- rep(factor("all"), ncol(object))
    }else{
      return(NULL)
    }
  }
  if(is.character(grouping) && length(grouping)==1 && 
     .isSElike(object) &&
     grouping %in% colnames(colData(object))){
    grouping <- colData(object)[[grouping]]
  }
  stopifnot(length(grouping)==ncol(object))
  factor(grouping)
}

#' normalizeDevsForSize
#' 
#' Normalizes the z-scores assay of a deviations object to make the scores 
#' comparable across motifs with different number of matches.
#' 
#' @details
#' If a TF changes the accessibility of all its motif matches by the same 
#' relative amount, the z-score of its motif increases with the square root of
#' the effective number of matches, `Neff` (see the `rowData` of the object).
#' `Neff` is smaller than the number of matches `N` for motifs matching a large 
#' share of the peaks in the bins in which they are found, since the 
#' expectation of the cells is estimated from the same peaks. The z-scores are 
#' multiplied by `sqrt(median(Neff)/Neff)`. For objects without `Neff` (e.g. 
#' produced by earlier versions) `N` is used instead.
#'
#' @param dev A SummarizedExperiment object as produced by 
#'   \code{\link{betterChromVAR}} or \code{\link{computeDeviationsAnalytic}}.
#'
#' @returns The `dev` object with an additional assay named 'norm'.
#' @export
#'
#' @examples
#' attach(getDummyData())
#' dev <- betterChromVAR(counts, motifMatches)
#' dev <- normalizeDevsForSize(dev)
#' dev
normalizeDevsForSize <- function(dev){
  stopifnot(inherits(dev, "SummarizedExperiment"))
  stopifnot("z" %in% assayNames(dev))
  N <- rowData(dev)$Neff
  if(is.null(N)) N <- rowData(dev)$N
  stopifnot(!is.null(N))
  assay(dev, "norm") <- assay(dev, "z")*sqrt(median(N)/N)
  dev
}


.packageDevSE <- function(a, object, motifCD, d){
  if(!is.null(motifCD)) d <- cbind(motifCD, d)
  a <- a[intersect(c("deviations","z"),names(a))]
  SummarizedExperiment(
    assays = a,
    colData = colData(object),
    rowData = d,
    metadata = metadata(object)
  )  
}

.checkAnnotations <- function(annotations){
  stopifnot(length(dim(annotations))==2)
  if(max(annotations) > 1 || min(annotations)<0)
    warning("`annotations` should be either binary or weights from 0 to 1.")
}

.isSElike <- function(object){
  inherits(object, "SummarizedExperiment") ||
    inherits(object, "SingleCellExperiment") ||
    inherits(object, "RangedSummarizedExperiment")
}


.bpFromThreads <- function(nthreads, progressbar=FALSE, default=1){
  if(is.null(nthreads)) nthreads <- default
  if(inherits(nthreads, "BiocParallelParam")) return(nthreads)
  if(is.numeric(nthreads) && length(nthreads)==1 && nthreads>0 &&
     nthreads==round(nthreads)){
    # make the global seed impact the BiocParallel one:
    seed <- sample.int(1000, 1L) 
    if(nthreads==1) return(SerialParam(progressbar=progressbar, RNGseed=seed))
    return(MulticoreParam(as.integer(nthreads), progressbar=progressbar, 
                          RNGseed=seed))
  }
  stop("`nthreads` should either be a positive integer, or a ",
       "BiocParallelParam object.")
}

# Validates and expands the `kappa` argument (scalar, vector with one value per
# motif, or output of estimateNoiseKappa)
.checkKappa <- function(kappa, n){
  if(is.data.frame(kappa)) kappa <- kappa$kappa
  if(!is.numeric(kappa) || !(length(kappa) %in% c(1L, n)) || anyNA(kappa) ||
     any(kappa<=0))
    stop("`kappa` should be a positive number, or a vector of positive ",
         "numbers with one value per row of `dev`.")
  rep(kappa, length.out=n)
}

.asDgC <- function(x){
  if(is(x, "dgCMatrix")) return(x)
  if(is(x, "DelayedMatrix")) x <- as.matrix(x)
  as(as(x, "dMatrix"), "CsparseMatrix")
}
