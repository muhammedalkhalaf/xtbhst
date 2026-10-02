## xtbhst 1.1.0

This release corrects the computations below; the version on CRAN is 1.0.2.

- Unit-specific error variance. The previous versions divided the pooled
  fixed-effects residual sum of squares of each unit by T minus the number
  of partialled-out variables, which is neither the estimator of Blomquist
  and Westerlund (2016) nor that of Pesaran and Yamagata (2008). A new
  argument `variance = c("bw", "py")` selects the estimator: `"bw"`
  (default) is T^-1 times the residual sum of squares of the unit-specific
  least-squares regression (Blomquist and Westerlund, 2016, Section 3.1);
  `"py"` is the pooled fixed-effects residual sum of squares divided by
  T - K - 1 (Pesaran and Yamagata, 2008). The chosen estimator is used in
  the weights of the weighted fixed-effects estimator, in S and in every
  bootstrap statistic.
- Adjusted Delta. The variance in the denominator was
  2K(T - K - Kpartial - 1)/(T - Kpartial + 1), where Kpartial counted the
  constant; it is now 2K(T - K - 1)/(T + 1) as in Pesaran and Yamagata
  (2008). The silent fallback to 2K when this was non-positive has been
  removed; `delta_adj` is `NA` with a warning when T - K - 1 <= 0.
- p-values. The bootstrap p-value `pval` is now based on the statistic S
  that the paper bootstraps (proportion of S* at least as large as S).
  Delta and adjusted Delta are fixed rescalings of S, so their bootstrap
  p-values coincided with each other; the duplicated `pval_adj` and
  `delta_adj_stars` elements have been removed and `pval_delta_asy` and
  `pval_delta_adj_asy` (asymptotic standard normal p-values) are reported
  instead. Delta and adjusted Delta are always computed with the Pesaran
  and Yamagata (2008) variance estimator (from the new element `S_py`),
  whatever the value of `variance`, so that they are the statistics of
  that paper; `variance` governs only `S` and the bootstrap. The
  conclusion in `summary()` uses the bootstrap p-value only. New elements
  `S`, `S_py`, `S_stars`, `sigma2`, `variance` and `n_dropped`.
- Lagged cross-sectional averages (`csa_lags > 0`). Rows with missing
  lagged averages were previously set to zero (including the constant)
  without reducing T, which distorted the unit estimates. The first
  `csa_lags` periods are now dropped for all units, T is reduced
  accordingly, and when the dependent variable enters `csa` its
  cross-sectional averages are recomputed from the pseudo-data in each
  bootstrap replication. The documentation states that CSA augmentation is
  an extension not covered by the paper.
- Block length. The default is now `round(2 * T^(1/3))`, which gives the
  block lengths 6, 7 and 9 for T = 25, 50 and 100 stated in Section 4.2 of
  the paper (previously `floor()`, giving 5, 7 and 9). A `blocklength`
  larger than T is now an error.
- Missing values. Rows with missing values in any variable used are
  dropped before the balance check, and the check now also detects
  duplicated (id, time) pairs; previously missing values in the model
  variables caused a "subscript out of bounds" error.
- Singular unit regressions now stop with an error naming the unit
  instead of silently applying a ridge correction.
- Citation corrected everywhere to Blomquist and Westerlund (2016),
  Empirical Economics, 50(4), 1359-1381; the reference to Pesaran and
  Yamagata (2008) has been added.
- `plot()`: plot 1 now shows the bootstrap distribution of S and plot 2
  that of Delta.
- Removed the unused suggested dependency on 'plm'.

## Test environments

* Ubuntu 24.04, R 4.3.3 and R-devel, R CMD check --as-cran

## R CMD check results

0 errors | 0 warnings | 0 notes
