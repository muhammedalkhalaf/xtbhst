#' xtbhst: Bootstrap Slope Heterogeneity Test for Panel Data
#'
#' Implements the bootstrap slope heterogeneity test for panel data based on
#' Blomquist and Westerlund (2015). The test examines whether slope coefficients
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
#' The test is based on the Delta statistic:
#' \deqn{\Delta = \sqrt{N} \frac{\tilde{S}/N - K}{\sqrt{2K}}}
#'
#' where \eqn{\tilde{S}} measures the weighted sum of squared deviations
#' of individual slope estimates from the weighted pooled estimate.
#'
#' The bootstrap procedure uses a block bootstrap to preserve serial
#' correlation in the residuals.
#'
#' @section Cross-Sectional Dependence:
#' The package supports handling cross-sectional dependence through
#' cross-sectional averages (CSA), following the approach of Pesaran (2006).
#' Users can specify variables for which CSA should be computed and
#' partialled out.
#'
#' @references
#' Blomquist, J., & Westerlund, J. (2015).
#' Panel bootstrap tests of slope homogeneity.
#' \emph{Empirical Economics}, 48(1), 1191-1204.
#' \doi{10.1007/s00181-015-0978-z}
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
