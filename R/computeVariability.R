#' computeMotifVariability
#'
#' @param z A matrix of z-scores (with motifs as rows), or a 
#'   `SummaziedExperiment` object with such a matrix in assay `z` (as produced
#'   e.g. by \code{\link{computeDeviationsAnalytic}}).
#' @param confInt The the confidence interval (a numeric scalar between 0 and 
#'   1, default 0.95).
#' @param method The method used to compute the confidence interval. 'normal'
#'   computes it analytically, assuming that the z-scores are normally 
#'   distributed. 'bonett' (default) adjusts this analytic estimate for 
#'   kurtosis (based on Bonett, Computational Statistics & Data Analysis, 2006).
#'   'bootstrap' uses bootstrapping, which does not scale very well.
#' @param n The number of bootstrap samples. Ignored unless 
#'   `method="bootstrap"`.
#'
#' @returns A data.frame containing the variability and confidence interval 
#'   around it, as well as significance, for each motif. If `z` is a 
#'   `SummarizedExperiment`, the data.frame will be stored in `rowData(z)`.
#' @seealso \code{\link{estimateNoiseKappa}}.
#' @export
#' @importFrom stats qchisq
#' @importFrom matrixStats rowSds rowQuantiles
#' @importFrom MatrixGenerics rowMeans2
#' @examples
#' # we generate random z-scores:
#' z <- matrix(rnorm(mean=rnorm(10), sd=runif(10, max=2), 200), nrow=10)
#' var <- computeMotifVariability
computeMotifVariability <- function(z, confInt=0.95, n=100,
                                    method=c("bonett", "normal", "bootstrap")){
  method <- match.arg(method)
  stopifnot(length(confInt)==1 && confInt>0 && confInt<1)

  if(.isSElike(z)){
    out <- computeMotifVariability(assay(z, "z"), confInt=confInt, n=n, method=method)
    rowData(z) <- cbind(
      rowData(z)[,setdiff(colnames(rowData(z)),colnames(out)),drop=FALSE], out)
    return(z)
  }
  
  out <- if(method=="bonett"){
    computeBonettVariability(z, confInt)
  }else if(method=="bootstrap"){
    computeBootstrapVariability(z, confInt, n=n)
  }else{
    computeNormalVariability(z, confInt)
  }
  out
}

computeNormalVariability <- function(z, confInt=0.95){
  df <- ncol(z)-1L
  variability <- matrixStats::rowSds(as.matrix(z))
  vars <- variability^2
  
  alpha2 <- (1-confInt)/2
  chi_ci <- stats::qchisq(c(1-alpha2, alpha2), df=df)
  
  pval <- stats::pchisq(df*vars, df=df, lower.tail=FALSE)
  
  data.frame(
    variability = variability,
    var.lower = sqrt(df*vars/chi_ci[1]),
    var.upper = sqrt(df*vars/chi_ci[2]),
    var.pval=pval,
    var.adjPval=p.adjust(p=pval, method="BH"),
    row.names = rownames(z)
  )
}


# Compute variability using Bonett's Method (Kurtosis-adjusted variability)
# 
# @param z_scores A dense or sparse M x C matrix of analytical Z-scores.
# @param confInt The confidence interval (default 0.95).
# @return A data.frame with variability (SD), kurtosis, and confidence bounds.
computeBonettVariability <- function(z, confInt=0.95){
  z <- as.matrix(z)
  nC <- ncol(z)
  SD <- matrixStats::rowSds(z)
  vars <- SD^2
  z <- z - rowMeans2(z)
  
  # Calculate kurtosis
  m2 <- rowMeans2(z^2)
  kurtosis <- ifelse(m2 == 0, 3, rowMeans2(z^4) / (m2^2))
  
  # Bonett's CI calculation
  # Get the standard normal critical value
  alpha <- 1 - confInt
  z_crit <- stats::qnorm(1 - alpha / 2)
  
  # Bonett's standard error of the log-variance
  se_ln_var <- sqrt( (kurtosis - (nC-3)/nC) / (nC-1) )
  # Bonett's small-sample adjustment constant 'c'
  # (approaches 1 for large nC)
  c <- nC / (nC - z_crit)
  ln_var <- log(c * pmax(vars, 1e-16))
  
  # Compute CI for the Standard Deviation 
  ci_lower <- exp(0.5 * (ln_var - z_crit * se_ln_var))
  ci_upper <- exp(0.5 * (ln_var + z_crit * se_ln_var))
  
  # Correct the bounds for motifs with absolutely zero variance
  ci_lower[which(vars==0)] <- 0
  ci_upper[which(vars==0)] <- 0
  
  pval <- stats::pnorm(ln_var/se_ln_var, lower.tail=FALSE)
  
  data.frame(
    variability=SD,
    var.lower=ci_lower,
    var.upper=ci_upper,
    var.pval=pval,
    var.adjPval=p.adjust(p=pval, method="BH"),
    row.names=rownames(z)
  )
}

computeBootstrapVariability <- function(z, confInt=0.95, n=100){
  z <- as.matrix(z)
  bsds <- vapply(seq_len(n), FUN.VALUE=numeric(nrow(z)), FUN=\(i){
    idx <- sample.int(ncol(z), size=ncol(z), replace=TRUE)
    matrixStats::rowSds(z[, idx])
  })
  alpha2 <- (1-confInt)/2
  CI <- matrixStats::rowQuantiles(bsds, probs=c(alpha2, 1-alpha2))
  pval <- (rowSums(bsds <= 1)+1) / (ncol(bsds)+1)
  data.frame(
    variability=matrixStats::rowSds(z),
    var.lower=CI[,1],
    var.upper=CI[,2],
    var.pval=pval,
    var.adjPval=p.adjust(p=pval, method="BH"),
    row.names=rownames(z)
  )
}
