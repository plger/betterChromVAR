# Decorrelation of motifs sharing peaks (ridge deconvolution)
#
# The motif sums of a cell are S = T x, with x the counts of the peaks and
# T = A' - N' P B the annotation after subtraction of the background (A: peaks
# x motifs annotation, B: bins x peaks membership, N = B A, P: bin-to-bin
# probabilities; i.e. column m of T' is c_qm = a_qm - g_b(q),m, see .fpcStats).
# If the accessibility of the peaks of each motif m changes by a relative
# amount beta_m, x = sf * e * (1 + A beta), and (the expectation e being itself
# removed by T) S is approximately sf * O beta, with
#   O = T diag(e) T'
# the overlap between motifs, after background subtraction and weighted by the
# expected accessibility of the peaks. With counts whose variance is
# proportional to their expectation, O is also proportional to the covariance
# of the motif sums under the null.
# `.motifOverlap` returns O (motifs x motifs).
#' @importFrom Matrix crossprod sparseMatrix
.motifOverlap <- function(annotations, background,
                          expectation=background@expectation, st=NULL){
  if(is.null(st)) st <- .fpcStats(annotations, background)
  bin2peak <- sparseMatrix(i=background@peak2bin,
                           j=seq_along(background@peak2bin),
                           dims=c(length(background@binDensity),
                                  nrow(annotations)))
  Q <- as.matrix(bin2peak %*% (annotations * expectation))
  binExp <- as.numeric(bin2peak %*% expectation)
  Qg <- crossprod(Q, st$g)
  O <- .weightedCrossprod(annotations, expectation) - Qg - t(Qg) +
    crossprod(st$g, binExp * st$g)
  O <- (O + t(O))/2
  dimnames(O) <- list(colnames(annotations), colnames(annotations))
  O
}

# t(x) %*% diag(w) %*% x for non-negative weights w, computed as the
# (symmetric) cross-product of x*sqrt(w). Annotations that are not sparse (e.g.
# graded matches) are processed as dense blocks of rows (see .annotCrossprod).
.weightedCrossprod <- function(x, w, maxDensity=0.2, blockCells=2.5e7){
  stopifnot(all(w>=0))
  w <- sqrt(w)
  if(!inherits(x, "sparseMatrix")) return(crossprod(x * w))
  if(.isSparseAnnotation(x, maxDensity))
    return(as.matrix(Matrix::crossprod(x * w)))
  x <- .asDgC(x)
  Reduce("+", lapply(.rowBlocks(x, blockCells), function(i){
    crossprod(as.matrix(x[i,,drop=FALSE]) * w[i])
  }))
}

# Ridge deconvolution of the motif sums:
#   L = (1 + ridge) D (O + ridge*D)^-1,   D = diag(O)
# so that L %*% S are motif sums in the original units, from which what can be
# attributed to the other motifs was removed.
# `varRatio` is the factor by which the null variance of the motif sums changes
# (taking their covariance as proportional to O).
.ridgeOperator <- function(annotations, background, ridge,
                           expectation=background@expectation, st=NULL){
  if(is.list(ridge) && !is.null(ridge$L)) return(ridge) # already computed
  if(!is.numeric(ridge) || length(ridge)!=1L || is.na(ridge) || ridge<=0)
    stop("`ridge` should be NULL or a single positive number.")
  if(is.infinite(ridge)) return(NULL)
  O <- .motifOverlap(annotations, background, expectation, st=st)
  D <- diag(O)
  ok <- which(D > 0)
  L <- diag(nrow(O))
  varRatio <- rep(1, nrow(O))
  if(length(ok)>0){
    Ook <- O[ok,ok,drop=FALSE]
    Lok <- (1 + ridge) * D[ok] * solve(Ook + diag(ridge*D[ok], length(ok)))
    L[ok,ok] <- Lok
    varRatio[ok] <- rowSums((Lok %*% Ook) * Lok)/D[ok]
  }
  dimnames(L) <- dimnames(O)
  list(L=L, varRatio=varRatio, ridge=ridge)
}
