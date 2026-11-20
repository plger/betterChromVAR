#' NR3C1example
#'
#' A small example dataset of ATAC-seq data upon ligand-based activation of the 
#' glucocorticoid receptor (GCR, encoded by the NR3C1 gene) in A549 cells 
#' (a subset of the data from McDowell et al., 2018: three control and three 
#' treated samples). Loading it (`data(NR3C1example)`) provides the following objects:
#' 
#' * `peakSE`: a `RangedSummarizedExperiment` of peak counts (assay `counts`), 
#'   with a `bias` column in the `rowData` and the `group` and `depth` columns in
#'   the `colData`.
#' * `motifs`: a `PFMatrixList` of three motifs (GCR, JUND and CTCFL), for use 
#'   with functions that need the motif themselves 
#'   (e.g. \code{\link{exploreMotifInteraction}}).
#' * `didevs`: a list with the output of 
#'   \code{\link{exploreMotifInteraction}} for the pairs GCR-JUND and GCR-CTCFL
#'   (precomputed on the hg38 genome, which is needed to run it).
#' * `motifMatches`: a `RangedSummarizedExperiment` with a logical assay 
#'   (`motifMatches`) of the matches of 29 motifs in each peak of `peakSE`.
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
