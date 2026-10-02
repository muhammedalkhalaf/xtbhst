#' xtbhst: Bootstrap Slope Heterogeneity Test for Panel Data
#'
#' Implements the bootstrap slope heterogeneity test for panel data based on
#' Blomquist and Westerlund (2016). The test examines whether slope coefficients
#' are homogeneous across cross-sectional units in a panel regression.
#'
#' @section Main Function:
#' \describe{
#'   \item{\code{\link{xtbhst}}}{Performs the bootstrap slope heterogeneity test.}
#' }
#'
#' @section Methods:
#' \describe{
#'   \item{\code{\link{print.xtbhst}}}{Print test results.}
#'   \item{\code{\link{summary.xtbhst}}}{Detailed summary of test results.}
#'   \item{\code{\link{plot.xtbhst}}}{Diagnostic plots for the bootstrap test.}
#' }
#'
#' @section Test Details:
#' The null hypothesis is that all cross-sectional units share the same
#' slope coefficients (H0: homogeneous slopes). Rejection indicates
#' significant slope heterogeneity, suggesting that pooled or fixed effects
#' estimators may be inappropriate.
#'
#' The bootstrapped statistic is the Swamy-type statistic \eqn{S} of
#' Blomquist and Westerlund (2016), a weighted sum of squared deviations
#' of the individual slope estimates from the weighted fixed-effects
#' estimate. The bootstrap procedure resamples the unit-specific residuals
#' in blocks along the time dimension, which preserves serial correlation
#' and cross-sectional dependence. The standardized statistics
#' \eqn{\Delta} and \eqn{\Delta_{adj}} of Pesaran and Yamagata (2008)
#' are reported with asymptotic normal p-values for reference. See
#' \code{\link{xtbhst}} for the formulas.
#'
#' @section Cross-Sectional Dependence:
#' The bootstrap itself is robust to cross-sectional dependence. As an
#' extension not covered by Blomquist and Westerlund (2016), the package
#' can also partial out cross-sectional averages (CSA) of user-specified
#' variables, in the spirit of Pesaran (2006).
#'
#' @references
#' Blomquist, J. and Westerlund, J. (2016).
#' Panel bootstrap tests of slope homogeneity.
#' \emph{Empirical Economics}, 50(4), 1359-1381.
#' \doi{10.1007/s00181-015-0978-z}
#'
#' Pesaran, M. H. and Yamagata, T. (2008).
#' Testing slope homogeneity in large panels.
#' \emph{Journal of Econometrics}, 142(1), 50-93.
#' \doi{10.1016/j.jeconom.2007.05.010}
#'
#' Pesaran, M. H. (2006).
#' Estimation and inference in large heterogeneous panels with a
#' multifactor error structure.
#' \emph{Econometrica}, 74(4), 967-1012.
#' \doi{10.1111/j.1468-0262.2006.00692.x}
#'
#' Swamy, P. A. V. B. (1970).
#' Efficient inference in a random coefficient regression model.
#' \emph{Econometrica}, 38(2), 311-323.
#' \doi{10.2307/1913012}
#'
#' @keywords internal
"_PACKAGE"
