# betterChromVAR 1.99.4

## Changes
* Simplified the adaptive kernel width of `getBackgroundBins()` using a minimum
  peak counts (the `maxShare` argument is deprecated, and `minCount` produces 
  virtually the same results).

## New features
* `CVnorm()` now warns when the experimental groups explain much of the bin-level
  bias factors.
* Reimplementation of `GCQN()` for convenience.
* `testBiasDimensions()`: cross-validated test of the importance of the different
  bias dimensions (accessibility, GC, fragment-length) for (pseudo)bulk data, to
  guide the choice of normalization (GC-only vs `CVnorm()`).


# betterChromVAR 1.99.1

## New features
* `fitDeviations()` (moved from `weightedMotifAccess`)
* `peakLevelMotifTests()` and `stepwisePeakLevelMotifTests()` as peak-logFC-based
  alternative to the analysis of deviations.
* Motif interaction analysis, including functions moved from `weightedMotifAccess`
  (and improved), as well as peak-level functions and a vignette.

# betterChromVAR 1.1.14

## New features
* `ridge` argument of `betterChromVAR()` and `computeDeviationsAnalytic()`: 
  optional decorrelation of the motifs that share peaks (ridge deconvolution 
  based on the overlap between motifs after background subtraction).
* `getGradedMatches()`: graded (instead of binary) motif matches, giving
  partial weights to weak matches according to their score.

## Changes
* Removed the bias shrinkage options -- see downstream smoothing instead.
* The bin-to-bin background probabilities of `getBackgroundBins()` are now 
  balanced (`tol` argument; `tol=NULL` restores the earlier behavior).
  This slightly changes the results, removing a bias in mean z.
* Minor gains in efficiency.

# betterChromVAR 1.1.13
* added pseudoCount and `loo` to make CVnorm more robust.
* added adaptive kernel width

# betterChromVAR 1.1.12

## New features
* `smoothDeviationsKNN()`: adaptive EB smoothing of single-cell deviations
  across nearest-neighbor cells.

## Bug fixes
* Fixed wrong denominator in `computeDeviationsWeighted`
* Fixed mistake in background variance (affecting the scale of the z-scores):
  it is now the variance of the counts of each bin's background peaks (as
  drawn by `sampleBackgroundPeaks()`), computed before any shrinkage

## Changes
* The z-scores (and things based on it, e.g. motif variability) now include a
  finite-population correction (`fpc=TRUE`) of the null variance.
* functions now also return an effective `N` per motif that is used in
  `normalizeDevsForSize` for more accurate adjustment.

# betterChromVAR 1.1.9
* Fixed bug that some SE-like objects were not recognized as such

# betterChromVAR 1.0.8
* Added methods for kNN-based background

# betterChromVAR 1.0.2

## Bug fixes
* Fixed namespace error in getBackgroundBins when using dgCMatrix counts
* Fixed bug in the binning when using frag length bias
* Fixed delay in betterChromVAR

# betterChromVAR 1.0.0

Initial release.
