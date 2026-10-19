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
