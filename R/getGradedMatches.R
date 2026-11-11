#' getGradedMatches
#'
#' Scans regions for motif matches and returns graded (rather than binary)
#' matches, in which weaker matches are given a weight between 0 and 1
#' according to their score. This is a wrapper around
#' \code{\link[motifmatchr]{matchMotifs}}, which needs to be installed.
#'
#' @param pwms The motifs, in any of the formats accepted by
#'   \code{\link[motifmatchr]{matchMotifs}} (e.g. a `PWMatrixList`).
#' @param subject The regions in which to search for matches, in any of the
#'   formats accepted by \code{\link[motifmatchr]{matchMotifs}} (e.g. a
#'   `GRanges`, a `RangedSummarizedExperiment`, or a `DNAStringSet`).
#' @param genome The genome (e.g. a BSgenome object), needed if `subject` does
#'   not contain the sequences.
#' @param pCutoffs A numeric vector of length 2 indicating the p-value
#'   thresholds respectively below which a match starts to receive a weight
#'   (lenient threshold), and below which it receives a weight of 1 (stringent
#'   threshold).
#' @param minWeight Weights below this value are set to 0 (keeps the
#'   annotation sparser).
#' @param chunkSize The number of motifs processed at a time (this affects
#'   memory usage and not the results).
#' @param verbose Logical; whether to output progress messages.
#' @param ... Passed to \code{\link[motifmatchr]{matchMotifs}}.
#'
#' @details
#' For each motif and region, the weight is based on the score of the best match
#'  in the region: it is 0 at the score corresponding to the lenient threshold
#' (`pCutoffs[1]`), and increases linearly with the (log-odds) score up to 1 at
#' the score corresponding to the stringent threshold (`pCutoffs[2]`, by
#' default the default threshold of `motifmatchr`). Regions with a match
#' passing the stringent threshold therefore have the same value as with binary
#' matches.
#'
#' The two score thresholds of a motif are taken as the lowest scores observed
#' at each cutoff, so that `subject` should not be too small.
#'
#' The result is considerably denser than binary matches (e.g. 30-50% non-zero).
#' Graded matches are mostly useful in combination with the decorrelation of
#' the motifs (see the `ridge` argument of \code{\link{betterChromVAR}}).
#'
#' @returns A sparse matrix of weights between 0 and 1, with regions as rows
#'   and motifs as columns, which can be used as the `annotations` of
#'   \code{\link{betterChromVAR}} or \code{\link{computeDeviationsAnalytic}}.
#' @importFrom Matrix cbind2 drop0
#' @importFrom methods as is
#' @export
#'
#' @examples
#' # not run:
#' # matches <- getGradedMatches(pwms, se, genome=genome)
#' # dev <- betterChromVAR(se, matches, ridge=3)
getGradedMatches <- function(pwms, subject, genome=NULL,
                             pCutoffs=c(5e-4, 5e-5), minWeight=0.05,
                             chunkSize=100L, verbose=FALSE, ...){
  if(!requireNamespace("motifmatchr", quietly=TRUE))
    stop("This function requires the 'motifmatchr' package.")
  stopifnot(is.numeric(pCutoffs) && length(pCutoffs)==2L && all(pCutoffs>0))
  stopifnot(is.numeric(minWeight) && length(minWeight)==1L && minWeight>=0 &&
              minWeight<1)
  pCutoffs <- sort(pCutoffs, decreasing=TRUE)
  if(pCutoffs[1]==pCutoffs[2]) stop("The two `pCutoffs` should differ.")

  getScores <- function(i, p){
    if(is.null(genome)){
      m <- motifmatchr::matchMotifs(pwms[i], subject, out="scores", 
                                    p.cutoff=p, ...)
    }else{
      m <- motifmatchr::matchMotifs(pwms[i], subject, genome=genome, 
                                    out="scores", p.cutoff=p, ...)
    }
    .asDgC(motifmatchr::motifScores(m))
  }
  
  i <- seq_along(pwms)
  chunks <- split(i, ceiling(i/chunkSize))
  res <- lapply(seq_along(chunks), function(ci){
    if(verbose && length(chunks)>1)
      message("Scanning motifs (chunk ", ci, "/", length(chunks), ")")
    i <- chunks[[ci]]
    .gradeScores(getScores(i, pCutoffs[1]), getScores(i, pCutoffs[2]),
                 minWeight=minWeight)
  })
  res <- Reduce(cbind2, res)
  if(!is.null(names(pwms)) && !any(duplicated(names(pwms))))
    colnames(res) <- names(pwms)
  res
}

# Converts sparse matrices of the best match scores (regions x motifs) obtained
# at a lenient (`lenient`) and at a stringent (`stringent`) threshold into
# weights going from 0 at the lowest lenient score of the motif to 1 at its
# lowest stringent score.
.gradeScores <- function(lenient, stringent, minWeight=0.05){
  stopifnot(identical(dim(lenient), dim(stringent)))
  lenient <- .asDgC(lenient)
  stringent <- .asDgC(stringent)
  colMin <- function(x){
    vapply(seq_len(ncol(x)), FUN.VALUE=numeric(1), function(j){
      v <- x@x[seq.int(x@p[j]+1L, length.out=x@p[j+1L]-x@p[j])]
      if(length(v)==0) return(NA_real_)
      min(v)
    })
  }
  lo <- colMin(lenient)
  hi <- colMin(stringent)
  j <- rep(seq_len(ncol(lenient)), diff(lenient@p))
  w <- (lenient@x - lo[j])/(hi[j] - lo[j])
  # motifs without any stringent match get no weight; if the two thresholds
  # coincide the matches are binary
  w[!is.finite(w)] <- as.numeric(is.finite(hi[j]))[!is.finite(w)]
  w <- pmin(1, pmax(0, w))
  w[w < minWeight] <- 0
  lenient@x <- w
  drop0(lenient)
}
