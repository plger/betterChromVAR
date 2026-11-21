#' fitDeviations
#' 
#' A wrapper to fit a model to motif deviations, using limma to estimate 
#' posterior variability from a trend with the (log) effective number of 
#' matches.
#'
#' @param dev A SummarizedExperiment of motif deviations, as produced for
#'   instance by \code{\link{betterChromVAR}}.
#' @param design The design matrix of the experiment, as produced for instance 
#'   by \code{\link[stats]{model.matrix}}, or a formula (e.g. `~group`)
#'   using `colData` columsn of `dev`.
#' @param type The type of assay to use, either 'deviations' (default, 
#'   background-adjusted motif deviations), 'z' (z-scores), or 'norm' 
#'   (set size-normalized z-scores).
#' @param contrasts Optional contrasts, passed to 
#'   \code{\link[limma]{contrasts.fit}}.
#' @param robust Logical, passed to \code{\link[limma]{eBayes}}; whether the 
#'   fit of the priors should be robust (not recommended).
#' @param ... Other parameters passed to \code{\link[limma]{lmFit}}.
#'
#' @details
#' The prior variance of the empirical Bayes moderation is estimated as a 
#' function of the log of the effective number of matches of the motifs (`Neff`,
#' or `N` for objects that lack it), consistent with the scaling used by
#' \code{\link{normalizeDevsForSize}}. If there are fewer than 20 motifs, or 
#' no variation in set sizes, no trend is used.
#' 
#' @returns A fitted `MArrayLM` object, which can for instance be used with
#'   limma's \code{\link[limma]{topTable}} function to extract differential
#'   analysis results.
#' @export
#' @importFrom limma lmFit contrasts.fit eBayes
#' @importFrom stats model.matrix
#'
#' @examples
#' dd <- getDummyData()
#' dev <- betterChromVAR(dd$counts, dd$motifMatches)
#' dev$group <- rep(LETTERS[1:2], each=5)
#' fit <- fitDeviations(dev, ~group)
#' limma::topTable(fit)
#' 
#' # this is equivalent to:
#' design <- stats::model.matrix(~group, data=SummarizedExperiment::colData(dev))
#' fit <- fitDeviations(dev, design)
#' limma::topTable(fit)
fitDeviations <- function(dev, design, type=c("deviations","norm","z"),
                          contrasts=NULL, robust=FALSE, ...){
  stopifnot(inherits(dev, "SummarizedExperiment"))
  stopifnot(!is.null(rowData(dev)$N) || !is.null(rowData(dev)$Neff))
  type <- match.arg(type)
  if(inherits(design, "formula")){
    design <- model.matrix(design, data=as.data.frame(colData(dev)))
  }
  if(!("norm" %in% assayNames(dev)) && type=="norm"){
    dev <- normalizeDevsForSize(dev)
  }
  stopifnot(type %in% assayNames(dev))
  fit <- lmFit(assay(dev, type), design, ...)
  if(!is.null(contrasts)){
    fit <- contrasts.fit(fit, contrasts)
  }
  N <- rowData(dev)$Neff
  if(is.null(N)) N <- rowData(dev)$N
  logN <- log(pmax(N, 1))
  trend <- nrow(dev)>=20 && length(unique(logN))>1
  if(trend) trend <- logN
  eBayes(fit, trend=trend, robust=robust)
}
