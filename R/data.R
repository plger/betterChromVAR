#' NR3C1example
#'
#' A small example dataset of ATAC-seq data upon ligand-based activation of the 
#' glucocorticoid receptor (GCR, encoded by the NR3C1 gene) in A549 cells 
#' (a subset of the data from McDowell et al., 2018: three control and three 
#' treated samples). Loading it (`data(NR3C1example)`) provides four objects 
#' (see Format).
#'
#' @format Four objects:
#' * `peakSE`: a `RangedSummarizedExperiment` of 19,143 peaks by 6 samples, with
#'   peak counts (assay `counts`), a `bias` column in the `rowData` and the 
#'   `group` and `depth` columns in the `colData`.
#' * `motifMatches`: a `RangedSummarizedExperiment` of 19,143 peaks by 29 
#'   motifs, with a logical assay (`motifMatches`) indicating the matches of 
#'   each motif in each peak of `peakSE`.
#' * `motifs`: a `PFMatrixList` of three motifs (GCR, JUND and CTCFL), for use 
#'   with functions that need the motifs themselves 
#'   (e.g. \code{\link{exploreMotifInteraction}}).
#' * `didevs`: a list with the output of 
#'   \code{\link{exploreMotifInteraction}} for the pairs GCR-JUND and GCR-CTCFL
#'   (precomputed on the hg38 genome, which is needed to run it).
#'
#' @name NR3C1example
#' @aliases peakSE motifMatches motifs didevs
#' @docType data
#' @keywords datasets
#' @examples
#' data(NR3C1example)
#' peakSE
#' motifMatches
NULL
