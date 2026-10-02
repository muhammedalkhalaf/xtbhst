#' Bootstrap Slope Heterogeneity Test for Panel Data
#'
#' Implements the panel bootstrap test of slope homogeneity of
#' Blomquist and Westerlund (2016). The null hypothesis is that the slope
#' coefficients are equal across all cross-sectional units.
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...} specifying
#'   the dependent variable and regressors.
#' @param data A data frame containing the panel data.
#' @param id A character string specifying the name of the cross-sectional
#'   identifier variable.
#' @param time A character string specifying the name of the time variable.
#' @param reps Integer. Number of bootstrap replications (default: 999).
#' @param blocklength Integer. Block length for the block bootstrap. If
#'   \code{NULL} (default), set to \code{round(2 * T^(1/3))}, the
#'   deterministic rule of Blomquist and Westerlund (2016, Section 4.2),
#'   where \code{T} is the number of time periods in the estimation sample.
#'   Must lie between 1 and \code{T}.
#' @param variance Character. Which unit-specific error variance estimate
#'   to use in the weights of the weighted fixed-effects estimator and in
#'   the test statistic. \code{"bw"} (default) is the estimator of
#'   Blomquist and Westerlund (2016, Section 3.1); \code{"py"} is the
#'   estimator of Pesaran and Yamagata (2008). See Details.
#' @param partial Optional formula specifying variables to be partialled out
#'   unit by unit. For example, \code{~ z1 + z2}.
#' @param csa Optional formula specifying variables for which cross-sectional
#'   averages should be computed and partialled out, in the spirit of
#'   Pesaran (2006). See Details.
#' @param csa_lags Integer. Number of lags of the cross-sectional averages
#'   to include (default: 0).
#' @param constant Logical. If \code{TRUE} (default), a unit-specific
#'   constant is partialled out.
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class \code{"xtbhst"} containing:
#' \describe{
#'   \item{S}{The Swamy-type statistic \eqn{S} that is bootstrapped,
#'     computed with the variance estimator selected by \code{variance}.}
#'   \item{S_py}{The same statistic computed with the variance estimator
#'     of Pesaran and Yamagata (2008); it underlies \code{delta} and
#'     \code{delta_adj}. Equal to \code{S} when \code{variance = "py"}.}
#'   \item{pval}{Bootstrap p-value: the proportion of bootstrap statistics
#'     \eqn{S^*} that are at least as large as \eqn{S}.}
#'   \item{delta}{The Delta statistic of Pesaran and Yamagata (2008),
#'     computed from \code{S_py} whatever the value of \code{variance}.}
#'   \item{delta_adj}{The adjusted Delta statistic of Pesaran and
#'     Yamagata (2008), or \code{NA} if \eqn{T - K - 1 \le 0}.}
#'   \item{pval_delta_asy}{Asymptotic (standard normal, upper tail)
#'     p-value of \code{delta}.}
#'   \item{pval_delta_adj_asy}{Asymptotic (standard normal, upper tail)
#'     p-value of \code{delta_adj}.}
#'   \item{S_stars}{Vector of bootstrap statistics \eqn{S^*}.}
#'   \item{delta_stars}{\code{S_stars} on the Delta scale, that is
#'     \code{sqrt(N) (S_stars/N - K)/sqrt(2K)}.}
#'   \item{blocklength}{The block length used in the bootstrap.}
#'   \item{reps}{Number of bootstrap replications.}
#'   \item{variance}{The variance estimator used (\code{"bw"} or
#'     \code{"py"}).}
#'   \item{N}{Number of cross-sectional units.}
#'   \item{T}{Number of time periods in the estimation sample (after
#'     dropping the first \code{csa_lags} periods).}
#'   \item{K}{Number of regressors.}
#'   \item{Kpartial}{Number of partialled-out variables (including the
#'     constant).}
#'   \item{sigma2}{Vector of unit-specific error variance estimates.}
#'   \item{beta_i}{Matrix of individual slope estimates (N x K).}
#'   \item{beta_fe}{Vector of weighted fixed-effects (pooled) estimates.}
#'   \item{n_dropped}{Number of rows of \code{data} dropped because of
#'     missing values.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details
#' Let \eqn{\hat\beta_i} be the least-squares estimate of the slopes of
#' unit \eqn{i} after partialling out the constant (and any further
#' variables given in \code{partial} and \code{csa}), let \eqn{M} denote
#' the corresponding projection matrix, and let
#' \deqn{\hat\beta_{WFE} = \left(\sum_{i=1}^N \frac{x_i' M x_i}{\hat\sigma_i^2}\right)^{-1}
#'   \sum_{i=1}^N \frac{x_i' M y_i}{\hat\sigma_i^2}}{
#'   b_WFE = (sum_i x_i' M x_i / s_i^2)^-1 sum_i x_i' M y_i / s_i^2}
#' be the weighted fixed-effects estimator. The test statistic of
#' Blomquist and Westerlund (2016, Section 3.1) is
#' \deqn{S = \sum_{i=1}^N (\hat\beta_i - \hat\beta_{WFE})'
#'   \frac{x_i' M x_i}{\hat\sigma_i^2} (\hat\beta_i - \hat\beta_{WFE}).}{
#'   S = sum_i (b_i - b_WFE)' (x_i' M x_i / s_i^2) (b_i - b_WFE).}
#'
#' Two estimators of \eqn{\sigma_i^2} are available:
#' \describe{
#'   \item{\code{variance = "bw"}}{\eqn{\hat\sigma_i^2 = T^{-1} \sum_t
#'     \hat\varepsilon_{i,t}^2}, where \eqn{\hat\varepsilon_{i,t}} are
#'     the residuals of the unit-specific least-squares regression. This
#'     is the estimator used by Blomquist and Westerlund (2016).}
#'   \item{\code{variance = "py"}}{\eqn{\tilde\sigma_i^2 = (T - K - 1)^{-1}
#'     (y_i - x_i \hat\beta_{FE})' M (y_i - x_i \hat\beta_{FE})}, where
#'     \eqn{\hat\beta_{FE}} is the pooled fixed-effects estimator. This is
#'     the estimator used by Pesaran and Yamagata (2008). When further
#'     variables are partialled out, the degrees of freedom are
#'     \eqn{T - K - K_{partial}}, which reduces to \eqn{T - K - 1} in the
#'     standard case.}
#' }
#' The same choice is used for the weights of \eqn{\hat\beta_{WFE}}, in
#' \eqn{S} and in every bootstrap statistic \eqn{S^*}.
#'
#' The bootstrap follows Algorithm BOOT of Blomquist and Westerlund
#' (2016): the unit-specific least-squares residuals are resampled in
#' blocks of length \code{blocklength} along the time dimension (keeping
#' the cross-section intact), pseudo-data are generated under the null
#' hypothesis with slopes \eqn{\hat\beta_{WFE}}, and \eqn{S^*} is
#' computed exactly as \eqn{S}. The reported bootstrap p-value is the
#' proportion of \eqn{S^*} that are at least as large as \eqn{S}.
#'
#' For reference, the standardized statistics of Pesaran and Yamagata
#' (2008) are also reported with their asymptotic standard normal p-values:
#' \deqn{\Delta = \sqrt{N} \frac{N^{-1} S - K}{\sqrt{2K}}, \qquad
#'   \Delta_{adj} = \sqrt{N} \frac{N^{-1} S - K}{\sqrt{2K (T - K - 1)/(T + 1)}}.}{
#'   Delta = sqrt(N) (S/N - K) / sqrt(2K),
#'   Delta_adj = sqrt(N) (S/N - K) / sqrt(2K (T - K - 1)/(T + 1)).}
#' These are always computed from \eqn{S} evaluated with the Pesaran and
#' Yamagata (2008) variance estimator \eqn{\tilde\sigma_i^2} (returned as
#' \code{S_py}), whatever the value of \code{variance}, so that they are
#' the statistics of that paper with their asymptotic null distribution;
#' \code{variance} governs only \eqn{S} and the bootstrap. Because
#' \eqn{\Delta} and \eqn{\Delta_{adj}} are fixed rescalings of
#' \code{S_py}, a bootstrap p-value for them would coincide with the
#' bootstrap p-value of \eqn{S} when \code{variance = "py"}; only the
#' asymptotic p-values are therefore reported for them. \eqn{\Delta_{adj}} is \code{NA} when
#' \eqn{T - K - 1 \le 0}.
#'
#' Rows of \code{data} with missing values in any of the variables used
#' (the model, \code{partial}, \code{csa}, \code{id} and \code{time}) are
#' dropped before the panel structure is checked. The remaining panel must
#' be strongly balanced: every unit must be observed in exactly the same
#' set of time periods, and every (id, time) pair must occur once.
#'
#' Partialling out cross-sectional averages (\code{csa}, \code{csa_lags})
#' is an extension not covered by Blomquist and Westerlund (2016). With
#' \code{csa_lags > 0} the first \code{csa_lags} periods are dropped for
#' all units, so that the estimation sample is a common sample of
#' \code{T - csa_lags} periods. If the dependent variable appears
#' in \code{csa}, its cross-sectional averages are recomputed from the
#' bootstrap pseudo-data in each replication (observed values are used
#' for lagged averages that fall before the estimation sample). The
#' dependent variable must then be a plain column of \code{data}.
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
#' @examples
#' \donttest{
#' # Generate example panel data
#' set.seed(123)
#' N <- 20  # cross-sectional units
#' T_periods <- 30  # time periods
#'
#' # Homogeneous slopes (H0 is true)
#' data_hom <- data.frame(
#'   id = rep(1:N, each = T_periods),
#'   time = rep(1:T_periods, N),
#'   x = rnorm(N * T_periods)
#' )
#' data_hom$y <- 1 + 0.5 * data_hom$x + rnorm(N * T_periods)
#'
#' # Test for slope heterogeneity
#' result <- xtbhst(y ~ x, data = data_hom, id = "id", time = "time",
#'                  reps = 199, seed = 42)
#' print(result)
#' summary(result)
#'
#' # Pesaran and Yamagata (2008) variance estimator
#' result_py <- xtbhst(y ~ x, data = data_hom, id = "id", time = "time",
#'                     reps = 199, seed = 42, variance = "py")
#' }
#'
#' @export
xtbhst <- function(formula, data, id, time, reps = 999L, blocklength = NULL,
                   variance = c("bw", "py"),
                   partial = NULL, csa = NULL, csa_lags = 0L,
                   constant = TRUE, seed = NULL) {
  cl <- match.call()
  variance <- match.arg(variance)

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate inputs
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.")
  }
  if (length(formula) != 3L) {
    stop("'formula' must have a response, e.g. y ~ x1 + x2.")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  if (!is.character(id) || length(id) != 1L || !id %in% names(data)) {
    stop("Variable '", id, "' not found in data.")
  }
  if (!is.character(time) || length(time) != 1L || !time %in% names(data)) {
    stop("Variable '", time, "' not found in data.")
  }
  if (!is.null(partial) && !inherits(partial, "formula")) {
    stop("'partial' must be a formula (e.g., ~ z1 + z2).")
  }
  if (!is.null(csa) && !inherits(csa, "formula")) {
    stop("'csa' must be a formula (e.g., ~ x1 + x2).")
  }
  csa_lags <- as.integer(csa_lags)
  if (length(csa_lags) != 1L || is.na(csa_lags) || csa_lags < 0L) {
    stop("'csa_lags' must be a non-negative integer.")
  }
  if (is.null(csa) && csa_lags > 0L) {
    stop("'csa_lags' requires 'csa'.")
  }

  reps <- as.integer(reps)
  if (length(reps) != 1L || is.na(reps) || reps < 1L) {
    stop("'reps' must be a positive integer.")
  }

  # Variables used and missing-value handling: drop incomplete rows first
  vars_needed <- unique(c(all.vars(formula),
                          if (!is.null(partial)) all.vars(partial),
                          if (!is.null(csa)) all.vars(csa),
                          id, time))
  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0L) {
    stop("Variable(s) not found in data: ",
         paste(missing_vars, collapse = ", "), ".")
  }
  data <- as.data.frame(data)[, vars_needed, drop = FALSE]
  complete <- stats::complete.cases(data)
  n_dropped <- sum(!complete)
  data <- data[complete, , drop = FALSE]
  if (nrow(data) == 0L) {
    stop("No complete observations in 'data'.")
  }

  # Sort data by id and time
  data <- data[order(data[[id]], data[[time]]), , drop = FALSE]
  rownames(data) <- NULL
  id_var <- data[[id]]
  time_var <- data[[time]]

  # Check for a strongly balanced panel on the complete data
  units <- unique(id_var)
  N_g <- length(units)
  periods <- sort(unique(time_var))
  T_full <- length(periods)

  if (anyDuplicated(data[, c(id, time), drop = FALSE]) > 0L) {
    stop("xtbhst requires a strongly balanced panel. ",
         "Some (id, time) pairs occur more than once.")
  }
  unit_counts <- table(id_var)
  if (nrow(data) != N_g * T_full || any(unit_counts != T_full)) {
    stop("xtbhst requires a strongly balanced panel (every unit observed ",
         "in the same time periods). After dropping ", n_dropped,
         " row(s) with missing values, ", N_g, " units and ", T_full,
         " time periods were found, but ", nrow(data),
         " rather than ", N_g * T_full, " complete observations.")
  }
  if (N_g < 2L) {
    stop("At least two cross-sectional units are required.")
  }

  # Extract model matrices
  mf <- stats::model.frame(formula, data = data)
  Y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(formula, data = mf)

  # Remove intercept from X (partialled out if constant = TRUE)
  intercept_col <- which(colnames(X) == "(Intercept)")
  if (length(intercept_col) > 0L) {
    X <- X[, -intercept_col, drop = FALSE]
  }
  if (ncol(X) < 1L) {
    stop("At least one regressor is required.")
  }
  K <- ncol(X)
  xnames <- colnames(X)

  # Build partial variables matrix (Z), full sample
  Z <- NULL
  if (constant) {
    Z <- matrix(1, nrow = nrow(X), ncol = 1L)
    colnames(Z) <- "constant"
  }
  if (!is.null(partial)) {
    partial_X <- stats::model.matrix(partial, data = data)
    int_col <- which(colnames(partial_X) == "(Intercept)")
    if (length(int_col) > 0L) {
      partial_X <- partial_X[, -int_col, drop = FALSE]
    }
    if (ncol(partial_X) > 0L) {
      Z <- if (is.null(Z)) partial_X else cbind(Z, partial_X)
    }
  }

  # Cross-sectional averages (extension, not covered by the paper)
  csa_depends_on_y <- FALSE
  yname <- NULL
  if (!is.null(csa)) {
    csa_vars <- all.vars(csa)
    resp_vars <- all.vars(formula[[2L]])
    if (any(resp_vars %in% csa_vars)) {
      csa_depends_on_y <- TRUE
      if (!is.name(formula[[2L]])) {
        stop("When the dependent variable enters 'csa', it must be a plain ",
             "column of 'data' (not a transformation); add the transformed ",
             "variable to 'data' first.")
      }
      yname <- as.character(formula[[2L]])
    }
    csa_X <- stats::model.matrix(csa, data = data)
    int_col <- which(colnames(csa_X) == "(Intercept)")
    if (length(int_col) > 0L) {
      csa_X <- csa_X[, -int_col, drop = FALSE]
    }
    if (ncol(csa_X) > 0L) {
      csa_mat <- .compute_csa(csa_X, id_var, time_var, csa_lags)
      Z <- if (is.null(Z)) csa_mat else cbind(Z, csa_mat)
    } else {
      csa_depends_on_y <- FALSE
    }
  }
  Kpartial <- if (is.null(Z)) 0L else ncol(Z)

  # Common estimation sample: drop the first csa_lags periods for all units
  time_rank <- match(time_var, periods)
  keep <- time_rank > csa_lags
  T_g <- T_full - csa_lags
  if (T_g <= K + Kpartial) {
    stop("Too few time periods: T = ", T_g, " (after dropping ", csa_lags,
         " period(s) for lagged cross-sectional averages) but ",
         K + Kpartial, " coefficients are estimated per unit.")
  }

  # Function rebuilding Z from a full-length y vector (only needed when the
  # cross-sectional averages depend on the dependent variable)
  Z_fun <- NULL
  if (csa_depends_on_y) {
    Z_fixed <- Z[, seq_len(Kpartial - ncol(csa_mat)), drop = FALSE]
    Z_fun <- function(y_full) {
      d2 <- data
      d2[[yname]] <- y_full
      cx <- stats::model.matrix(csa, data = d2)
      ic <- which(colnames(cx) == "(Intercept)")
      if (length(ic) > 0L) cx <- cx[, -ic, drop = FALSE]
      cm <- .compute_csa(cx, id_var, time_var, csa_lags)
      Zf <- if (ncol(Z_fixed) > 0L) cbind(Z_fixed, cm) else cm
      Zf[keep, , drop = FALSE]
    }
  }

  # Set and validate block length
  if (is.null(blocklength)) {
    blocklength <- round(2 * T_g^(1 / 3))
  }
  blocklength <- as.integer(blocklength)
  if (length(blocklength) != 1L || is.na(blocklength) || blocklength < 1L) {
    stop("'blocklength' must be a positive integer.")
  }
  if (blocklength > T_g) {
    stop("'blocklength' (", blocklength, ") must not exceed the number of ",
         "time periods in the estimation sample (", T_g, ").")
  }

  # Run the bootstrap test on the estimation sample
  result <- .xtbhst_bootstrap(
    Y = Y[keep], X = X[keep, , drop = FALSE],
    Z = if (is.null(Z)) NULL else Z[keep, , drop = FALSE],
    id_var = id_var[keep], units = units,
    N_g = N_g, T_g = T_g, K = K, Kpartial = Kpartial,
    variance = variance, reps = reps, blocklength = blocklength,
    Z_fun = Z_fun, y_full = Y, keep = keep
  )

  beta_i <- result$beta_i
  colnames(beta_i) <- xnames
  rownames(beta_i) <- as.character(units)
  beta_fe <- as.vector(result$beta_fe)
  names(beta_fe) <- xnames

  out <- list(
    S = result$S,
    S_py = result$S_py,
    pval = result$pval,
    delta = result$delta,
    delta_adj = result$delta_adj,
    pval_delta_asy = result$pval_delta_asy,
    pval_delta_adj_asy = result$pval_delta_adj_asy,
    S_stars = result$S_stars,
    delta_stars = result$delta_stars,
    blocklength = blocklength,
    reps = reps,
    variance = variance,
    N = N_g,
    T = T_g,
    K = K,
    Kpartial = Kpartial,
    sigma2 = result$sigma2,
    beta_i = beta_i,
    beta_fe = beta_fe,
    n_dropped = n_dropped,
    formula = formula,
    partial = partial,
    csa = csa,
    csa_lags = csa_lags,
    constant = constant,
    call = cl
  )

  class(out) <- "xtbhst"
  out
}


#' Compute cross-sectional averages
#'
#' @param X Matrix of variables (rows sorted by id, then time; balanced).
#' @param id_var Vector of cross-sectional identifiers.
#' @param time_var Vector of time identifiers.
#' @param lags Number of lags to include.
#' @return Matrix of cross-sectional averages (with lags if requested);
#'   lagged values are \code{NA} in the first \code{lags} periods.
#' @keywords internal
#' @noRd
.compute_csa <- function(X, id_var, time_var, lags = 0L) {
  n_vars <- ncol(X)
  csa_mat <- matrix(NA_real_, nrow = nrow(X), ncol = n_vars)
  colnames(csa_mat) <- paste0("csa_", colnames(X))
  for (v in seq_len(n_vars)) {
    csa_mat[, v] <- stats::ave(X[, v], time_var)
  }

  if (lags > 0L) {
    lagged_list <- list(csa_mat)
    for (lag in seq_len(lags)) {
      lagged_mat <- matrix(NA_real_, nrow = nrow(X), ncol = n_vars)
      colnames(lagged_mat) <- paste0("csa_", colnames(X), "_L", lag)
      for (v in seq_len(n_vars)) {
        lagged_mat[, v] <- stats::ave(csa_mat[, v], id_var,
                                      FUN = function(z) {
                                        n <- length(z)
                                        if (n <= lag) return(rep(NA_real_, n))
                                        c(rep(NA_real_, lag), z[seq_len(n - lag)])
                                      })
      }
      lagged_list[[lag + 1L]] <- lagged_mat
    }
    csa_mat <- do.call(cbind, lagged_list)
  }

  csa_mat
}


#' Core computation of the statistic S and its ingredients
#'
#' @param Y Dependent variable vector (estimation sample, raw scale).
#' @param X Regressor matrix (raw scale).
#' @param Z Matrix of variables to partial out unit by unit, or NULL.
#' @param index List of row indices per unit.
#' @param units Unit labels (for error messages).
#' @param N_g,T_g,K,Kpartial Dimensions.
#' @param variance "bw" or "py".
#' @return List with S, beta_i, beta_wfe, beta_fe, sigma2, E (T x N
#'   matrix of unit-specific residuals), and the projections of Y and X.
#' @keywords internal
#' @noRd
.bhst_stats <- function(Y, X, Z, index, units, N_g, T_g, K, Kpartial,
                        variance) {
  Yt <- Y
  Xt <- X
  ZZinv_list <- vector("list", N_g)

  # Partial out Z unit by unit
  if (Kpartial > 0L) {
    for (i in seq_len(N_g)) {
      idx <- index[[i]]
      Zi <- Z[idx, , drop = FALSE]
      ZZinv <- .solve_unit(crossprod(Zi), units[i], "partialled-out variables")
      ZZinv_list[[i]] <- ZZinv
      Yt[idx] <- Y[idx] - Zi %*% (ZZinv %*% crossprod(Zi, Y[idx]))
      Xt[idx, ] <- X[idx, , drop = FALSE] -
        Zi %*% (ZZinv %*% crossprod(Zi, X[idx, , drop = FALSE]))
    }
  }

  XX_list <- vector("list", N_g)
  XXinv_list <- vector("list", N_g)
  XY_list <- vector("list", N_g)
  beta_i <- matrix(NA_real_, nrow = N_g, ncol = K)
  E <- matrix(NA_real_, nrow = T_g, ncol = N_g)
  sigma2_bw <- numeric(N_g)
  XX_sum <- matrix(0, K, K)
  XY_sum <- numeric(K)

  for (i in seq_len(N_g)) {
    idx <- index[[i]]
    Xi <- Xt[idx, , drop = FALSE]
    Yi <- Yt[idx]
    XX <- crossprod(Xi)
    XXinv <- .solve_unit(XX, units[i], "regressors")
    XY <- crossprod(Xi, Yi)
    XX_list[[i]] <- XX
    XXinv_list[[i]] <- XXinv
    XY_list[[i]] <- XY
    b_i <- XXinv %*% XY
    beta_i[i, ] <- as.vector(b_i)
    e_i <- Yi - Xi %*% b_i
    E[, i] <- e_i
    sigma2_bw[i] <- sum(e_i^2) / T_g
    XX_sum <- XX_sum + XX
    XY_sum <- XY_sum + as.vector(XY)
  }

  # Pesaran and Yamagata (2008): pooled fixed-effects residuals
  beta_fe <- solve(XX_sum, XY_sum)
  sigma2_py <- numeric(N_g)
  df_py <- T_g - K - Kpartial
  for (i in seq_len(N_g)) {
    idx <- index[[i]]
    r_i <- Yt[idx] - Xt[idx, , drop = FALSE] %*% beta_fe
    sigma2_py[i] <- sum(r_i^2) / df_py
  }
  sigma2 <- if (variance == "bw") sigma2_bw else sigma2_py
  for (s2 in list(sigma2, sigma2_py)) {
    if (any(!is.finite(s2)) || any(s2 <= 0)) {
      bad <- which(!is.finite(s2) | s2 <= 0)
      stop("The error variance estimate is zero or not finite for unit(s) ",
           paste(units[bad], collapse = ", "),
           "; the unit regression may fit the data exactly.")
    }
  }

  # Weighted fixed-effects estimator and Swamy-type statistic S for a
  # given vector of unit variances
  swamy <- function(s2) {
    A <- matrix(0, K, K)
    b <- numeric(K)
    for (i in seq_len(N_g)) {
      A <- A + XX_list[[i]] / s2[i]
      b <- b + as.vector(XY_list[[i]]) / s2[i]
    }
    bw <- solve(A, b)
    S <- 0
    for (i in seq_len(N_g)) {
      d <- beta_i[i, ] - bw
      S <- S + as.numeric(t(d) %*% XX_list[[i]] %*% d) / s2[i]
    }
    list(S = S, beta_wfe = bw)
  }
  main <- swamy(sigma2)
  S_py <- if (variance == "py") main$S else swamy(sigma2_py)$S

  list(S = main$S, S_py = S_py, beta_i = beta_i, beta_wfe = main$beta_wfe,
       sigma2 = sigma2, sigma2_py = sigma2_py, E = E, ZZinv_list = ZZinv_list)
}


#' Bootstrap core computation
#'
#' @param Y Dependent variable vector (estimation sample).
#' @param X Regressor matrix (estimation sample).
#' @param Z Matrix of variables to partial out (or NULL), estimation sample.
#' @param id_var Cross-sectional identifier vector (estimation sample).
#' @param units Unit labels.
#' @param N_g Number of cross-sectional units.
#' @param T_g Number of time periods in the estimation sample.
#' @param K Number of regressors.
#' @param Kpartial Number of partialled-out variables.
#' @param variance "bw" or "py".
#' @param reps Number of bootstrap replications.
#' @param blocklength Block length for block bootstrap.
#' @param Z_fun Function returning Z (estimation sample) from a full-length
#'   dependent variable, or NULL when Z does not depend on y.
#' @param y_full Full-length dependent variable (all periods).
#' @param keep Logical vector selecting the estimation sample in y_full.
#' @return List with test statistics and bootstrap results.
#' @keywords internal
#' @noRd
.xtbhst_bootstrap <- function(Y, X, Z, id_var, units,
                               N_g, T_g, K, Kpartial, variance,
                               reps, blocklength,
                               Z_fun = NULL, y_full = NULL, keep = NULL) {
  index <- lapply(units, function(u) which(id_var == u))

  st <- .bhst_stats(Y, X, Z, index, units, N_g, T_g, K, Kpartial, variance)
  S <- st$S
  S_py <- st$S_py
  beta_wfe <- st$beta_wfe
  E <- st$E

  # Delta statistics of Pesaran and Yamagata (2008), always computed from
  # S_py (their variance estimator), with asymptotic p-values
  delta <- sqrt(N_g) * (S_py / N_g - K) / sqrt(2 * K)
  pval_delta_asy <- stats::pnorm(delta, lower.tail = FALSE)
  if (T_g - K - 1 > 0) {
    var_adj <- 2 * K * (T_g - K - 1) / (T_g + 1)
    delta_adj <- sqrt(N_g) * (S_py / N_g - K) / sqrt(var_adj)
    pval_delta_adj_asy <- stats::pnorm(delta_adj, lower.tail = FALSE)
  } else {
    warning("The adjusted Delta statistic is not defined because ",
            "T - K - 1 = ", T_g - K - 1, " <= 0; it is reported as NA.")
    delta_adj <- NA_real_
    pval_delta_adj_asy <- NA_real_
  }

  # Unit-specific coefficients on Z under H0 (raw-scale bootstrap DGP):
  # y*_i = Z_i gamma_i + X_i beta_wfe + e*_i
  Zg <- numeric(length(Y))
  if (Kpartial > 0L) {
    for (i in seq_len(N_g)) {
      idx <- index[[i]]
      Zi <- Z[idx, , drop = FALSE]
      r <- Y[idx] - X[idx, , drop = FALSE] %*% beta_wfe
      gamma_i <- st$ZZinv_list[[i]] %*% crossprod(Zi, r)
      Zg[idx] <- Zi %*% gamma_i
    }
  }
  Xb <- as.vector(X %*% beta_wfe)

  # Bootstrap
  S_stars <- numeric(reps)
  n_blocks_start <- T_g - blocklength + 1L

  for (b in seq_len(reps)) {
    # Block bootstrap resampling of the T x N residual matrix
    E_star <- matrix(NA_real_, nrow = T_g, ncol = N_g)
    t_start <- 1L
    while (t_start <= T_g) {
      rand_idx <- sample.int(n_blocks_start, 1L)
      len <- min(blocklength, T_g - t_start + 1L)
      E_star[t_start:(t_start + len - 1L), ] <-
        E[rand_idx:(rand_idx + len - 1L), , drop = FALSE]
      t_start <- t_start + len
    }

    # Pseudo-data under H0
    Y_star <- Zg + Xb
    for (i in seq_len(N_g)) {
      idx <- index[[i]]
      Y_star[idx] <- Y_star[idx] + E_star[, i]
    }

    Z_star <- Z
    if (!is.null(Z_fun)) {
      y_full_star <- y_full
      y_full_star[keep] <- Y_star
      Z_star <- Z_fun(y_full_star)
    }

    S_stars[b] <- .bhst_stats(Y_star, X, Z_star, index, units,
                              N_g, T_g, K, Kpartial, variance)$S
  }

  pval <- mean(S_stars >= S)
  delta_stars <- sqrt(N_g) * (S_stars / N_g - K) / sqrt(2 * K)

  list(
    S = S,
    S_py = S_py,
    pval = pval,
    delta = delta,
    delta_adj = delta_adj,
    pval_delta_asy = pval_delta_asy,
    pval_delta_adj_asy = pval_delta_adj_asy,
    S_stars = S_stars,
    delta_stars = delta_stars,
    sigma2 = st$sigma2,
    beta_i = st$beta_i,
    beta_fe = beta_wfe
  )
}


#' Matrix inversion with an informative error for singular unit matrices
#'
#' @param x A square cross-product matrix to invert.
#' @param unit Label of the cross-sectional unit.
#' @param what Description of the variables (for the error message).
#' @return The inverse of x.
#' @keywords internal
#' @noRd
.solve_unit <- function(x, unit, what) {
  tryCatch(
    solve(x),
    error = function(e) {
      stop("The cross-product matrix of the ", what, " is singular for unit '",
           unit, "'; the ", what, " are collinear (or constant) within this ",
           "unit. ", conditionMessage(e), call. = FALSE)
    }
  )
}
