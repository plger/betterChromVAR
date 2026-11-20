#' motifCoOccurence
#' 
#' Returns regions that have matches for given pairs of motifs within certain 
#' distances of each other.
#'
#' @param motifs A named `PFMatrixList` or `PWMatrixList` object containing 
#'   motifs (only those specified in `pairs` will be used). If you're not 
#'   familiar with these objects, see the `TFBSTools` package, and the 
#'   `universalmotif` package on how to convert motifs.
#' @param pairs A list of pairs of motifs for which to compute co-occurence.
#'   Specifically, this should be a (optionally named) list of character vectors
#'   of length 2, corresponding to names in `motifs`.
#' @param regions The regions in which to search for motif matches.
#' @param genome A genome object or path to a genome fasta file.
#' @param centerDist Logical; whether to compute distances between the centers
#'   of the matches (default). If FALSE, the distance is that between the edges
#'   of the matches closest to each other (as in 
#'   `GenomicRanges::distanceToNearest`), i.e. 0 for overlapping matches, and
#'   so independent of the width of the motifs.
#' @param minDist The minimum distance between matches for a co-occurrence. This
#'   is primarily used to exclude interactions between highly-similar, 
#'   overlapping motifs, and to focus on putative cooperative interactions 
#'   between TFs. If interested in competing TFs, set this to 0. Also used
#'   (as the lower bound of the first bin) when using `nDistQuantiles`.
#' @param maxDist The maximum distance between the matches for a co-occurence.
#' @param exclusiveDist Logical; whether the co-occurences should be counted 
#'  only for the smallest of the distance thresholds. Ignored if `minDist` and
#'  `maxDist` have a length of 1.
#' @param restrictToRegions Logical; whether to restrict the search for motifs 
#'  to `regions`. If FALSE (default), one of the motifs still needs to be within
#'  the regions, but the other can be +/- `maxDist` around the region.
#' @param nDistQuantiles The number of distance quantiles to use. Disabled by
#'   default and requires `exclusiveDist=TRUE`. When a positive integer, the 
#'   span from `minDist` to `maxDist` is split into `nDistQuantiles` number of 
#'   quantiles for each pair of motifs (in each region, the shortest distance
#'   between matches of the two motifs is used).
#' @param ignore.strand Logical; whether to ignore the strand for co-occurence
#'   (default TRUE).
#' 
#' @details
#' Note that both `minDist` and `maxDist`, rather than specifying a single 
#' threshold, can compute co-occurence for multiple thresholds (which is must
#' faster than running the function multiple times). `minDist` and `maxDist` 
#' should have the same length, and the corresponding entries will be used 
#' together to produce multiple matrices.
#' If `exclusiveDist=TRUE`, the distance threshold are exclusive, meaning that 
#' the higher threshold does not include co-occurences of the lower thresholds.
#' When using multiple exclusive distance bins, it is also possible to have the
#' bin boundaries adjusted for each pair of motifs so that each bin is equally
#' populated by setting `nDistQuantiles` (the breaks used can be retrieved from
#' `attr(results, "breaks")` ). If the distances are heavily tied, fewer bins 
#' than `nDistQuantiles` can result (the corresponding matrices are empty).
#' Distances are between the centers of the matches (when `centerDist=TRUE`) or
#' between their closest edges (otherwise), and a co-occurence is attributed 
#' only to the regions that contain one of the two matches (not to those lying 
#' in between).
#' Also note that all matches are stored in memory, so using this function 
#' across the entire genome is not advisable (unless for very few motifs).
#' 
#'
#' @returns A list of sparse logical matrices, with one matrix for each value 
#'   of `minDist`/`maxDist` (or each quantile bin).
#' @author Pierre-Luc Germain
#' @export
#' @importFrom IRanges findOverlapPairs IRanges overlapsAny
#' @importFrom S4Vectors queryHits subjectHits from to
#' @importFrom GenomicRanges findOverlaps GRanges granges resize seqnames start end width
#' @importFrom stats quantile
motifCoOccurence <- function(motifs, pairs, regions, genome, centerDist=TRUE,
                             minDist=5, maxDist=50, exclusiveDist=TRUE,
                             restrictToRegions=FALSE,
                             nDistQuantiles=NULL, ignore.strand=TRUE){
  
  if(!requireNamespace("motifmatchr", quietly=TRUE))
    stop("The `motifmatchr` package is required for this function.")
  stopifnot(is(motifs, "XMatrixList"))
  stopifnot(is.list(pairs) && all(lengths(pairs)==2))
  stopifnot(is(regions, "GRanges"))
  stopifnot(length(minDist)==length(maxDist))
  stopifnot(is.numeric(maxDist) && all(maxDist>1))
  stopifnot(is.numeric(minDist) && all(minDist>=0))
  
  if(!is.null(nDistQuantiles)){
    stopifnot(is.numeric(nDistQuantiles) && length(nDistQuantiles)==1 &&
                nDistQuantiles>1 && (nDistQuantiles %%1==0))
    if(!exclusiveDist)
      stop("`nDistQuantiles` only makes sense with `exclusiveDist=TRUE`")
    if(length(maxDist)>1)
      stop("`maxDist` and `minDist` should each have a length of 1 when using",
           " `nDistQuantiles`.")
  }
  if(!all(sapply(pairs, \(x) x %in% names(motifs)))){
    stop("Some of the motifs specific in `pairs` appear not to be in `motifs`.")
  }
  if(is.null(names(pairs)))
    pairs <- setNames(pairs, sapply(pairs, paste, collapse="+"))
  if(is.character(genome) && length(genome)==1)
    genome <- .getFaFile(genome)
  
  if(restrictToRegions){
    r2 <- regions
  }else{
    r2 <- resize(regions, width=width(regions)+2*max(maxDist), fix="center")
  }
  matches <- motifmatchr::matchMotifs(motifs[unique(unlist(pairs))], r2,
                         genome=genome, out="positions")
  
  if(centerDist) matches <- lapply(matches, resize, fix="center", width=1)
  
  if(is.null(nDistQuantiles)){
    distCrit <- setNames(seq_along(minDist), paste0(minDist,"<= d <=",maxDist))
    return(lapply(distCrit, \(i){
      as(sapply(pairs, \(x){
        o <- findOverlapPairs(matches[[x[1]]], matches[[x[2]]],
                              maxgap=maxDist[i]-centerDist, 
                              ignore.strand=ignore.strand)
        m1 <- S4Vectors::first(o)
        m2 <- S4Vectors::second(o)
        d <- .matchDist(m1, m2, centerDist)
        keep <- which(d<=maxDist[i] & (!exclusiveDist | d>=minDist[i]))
        # regions containing either of the two matches (not the span between)
        gr <- c(granges(m1[keep]), granges(m2[keep]))
        return(overlapsAny(regions, gr, ignore.strand=TRUE))
      }), "sparseMatrix")
    }))
  }

  # Using adaptative distance thresholds
  res <- lapply(pairs, \(x){
    
    out <- list(bins=vector(mode="integer", length=length(regions)), q=NULL)
    
    # all pairs of matches within maxDist of each other
    o <- findOverlaps(matches[[x[1]]], matches[[x[2]]], 
                      maxgap=maxDist-centerDist, ignore.strand=ignore.strand)
    if(length(o)==0) return(out)
    m1 <- matches[[x[1]]][queryHits(o)]
    m2 <- matches[[x[2]]][subjectHits(o)]
    d <- .matchDist(m1, m2, centerDist)
    keep <- which(d>=minDist & d<=maxDist)
    if(length(keep)==0) return(out)
    m1 <- m1[keep]
    m2 <- m2[keep]
    d <- d[keep]
    
    # attribute each pair to the regions containing one of the two matches,
    # and find for each region the shortest distance between the motifs
    o1 <- findOverlaps(m1, regions, ignore.strand=TRUE)
    o2 <- findOverlaps(m2, regions, ignore.strand=TRUE)
    rg <- c(to(o1), to(o2))
    if(length(rg)==0) return(out)
    dd <- c(d[from(o1)], d[from(o2)])
    ag <- tapply(dd, rg, min)
    if(length(ag)<2) return(out)
    
    # define the distance quantiles
    q <- unique(quantile(ag, c(0,seq_len(nDistQuantiles))/nDistQuantiles))
    if(length(q)==1) return(out)
    out$bins[as.integer(names(ag))] <- cut(as.numeric(ag), breaks=q, 
                                           labels=FALSE, include.lowest=TRUE)
    out$q <- as.numeric(q)
    out
  })
  res <- res[!sapply(res, \(x) is.null(x$q))]
  if(length(res)==0) stop("No co-occurrence found for the given pairs.")
  bins <- seq_len(nDistQuantiles);
  names(bins) <- paste0("bin", bins)
  a <- lapply(bins, \(i){
    as(sapply(res, \(x) x$bins==i), "sparseMatrix")
  })
  attr(a, "breaks") <- lapply(res, \(x) x$q)
  a
}

# genome fasta file (requires Rsamtools)
.getFaFile <- function(path){
  if(!requireNamespace("Rsamtools", quietly=TRUE))
    stop("The `Rsamtools` package is required to use a fasta file as genome.")
  Rsamtools::FaFile(path)
}

# distance between pairs of matches: between their centers if `centerDist` (in 
# which case the matches have been reduced to their center), otherwise between 
# their closest edges (0 if they overlap or are adjacent, as distanceToNearest)
.matchDist <- function(m1, m2, centerDist){
  if(centerDist) return(abs(start(m1)-start(m2)))
  pmax(0L, pmax(start(m1), start(m2)) - pmin(end(m1), end(m2)) - 1L)
}
