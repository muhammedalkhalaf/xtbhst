#' Bootstrap Slope Heterogeneity Test for Panel Data
#'
#' Implements the bootstrap slope heterogeneity test based on
#' Blomquist and Westerlund (2015). Tests the null hypothesis that slope
#' coefficients are homogeneous across cross-sectional units.
#'
#' @param formula A formula of the form \code{y ~ x1 + x2 + ...} specifying
#'   the dependent variable and regressors.
#' @param data A data frame containing the panel data.
#' @param id A character string specifying the name of the cross-sectional
#'   identifier variable.
#' @param time A character string specifying the name of the time variable.
#' @param reps Integer. Number of bootstrap replications (default: 999).
#' @param blocklength Integer. Block length for the block bootstrap. If
#'   \code{NULL} (default), set automatically to \code{floor(2 * T^(1/3))}
#'   where T is the number of time periods.
#' @param partial Optional formula specifying variables to be partialled out.
#'   For example, \code{~ z1 + z2}.
#' @param csa Optional formula specifying variables for which cross-sectional
#'   averages should be computed and included. Used to handle cross-sectional
#'   dependence following Pesaran (2006).
#' @param csa_lags Integer. Number of lags of the cross-sectional averages
#'   to include (default: 0).
#' @param constant Logical. If \code{TRUE} (default), include a constant
#'   (partialled out).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return An object of class \code{"xtbhst"} containing:
#' \describe{
#'   \item{delta}{The Delta test statistic.}
#'   \item{delta_adj}{The adjusted Delta test statistic.}
#'   \item{pval}{Bootstrap p-value for the Delta statistic.}
#'   \item{pval_adj}{Bootstrap p-value for the adjusted Delta statistic.}
#'   \item{blocklength}{The block length used in the bootstrap.}
#'   \item{reps}{Number of bootstrap replications.}
#'   \item{N}{Number of cross-sectional units.}
#'   \item{T}{Number of time periods.}
#'   \item{K}{Number of regressors.}
#'   \item{beta_i}{Matrix of individual slope estimates (N x K).}
#'   \item{beta_fe}{Vector of weighted fixed effects (pooled) estimates.}
#'   \item{delta_stars}{Vector of bootstrap Delta statistics.}
#'   \item{delta_adj_stars}{Vector of bootstrap adjusted Delta statistics.}
#'   \item{call}{The matched call.}
#' }
#'
#' @details
#' The test is based on the following Delta statistic:
#' \deqn{\Delta = \sqrt{N} \frac{\tilde{S}/N - K}{\sqrt{2K}}}
#' where \eqn{\tilde{S}} is a weighted sum of squared deviations of
#' individual slope estimates from the weighted pooled estimate.
#'
#' The adjusted Delta statistic corrects for small sample bias:
#' \deqn{\Delta_{adj} = \sqrt{N} \frac{\tilde{S}/N - K}{\sqrt{V}}}
#' where \eqn{V = 2K(T-K-K_{partial}-1)/(T-K_{partial}+1)}.
#'
#' Under the null hypothesis of slope homogeneity, both statistics
#' are asymptotically standard normal. The bootstrap procedure provides
#' finite-sample p-values.
#'
#' The function requires a strongly balanced panel (all units observed
#' for all time periods).
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
#' }
#'
#' @export
xtbhst <- function(formula, data, id, time, reps = 999L, blocklength = NULL,
                   partial = NULL, csa = NULL, csa_lags = 0L,
                   constant = TRUE, seed = NULL) {
  cl <- match.call()

  # Set seed if provided
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Validate inputs
  if (!inherits(formula, "formula")) {
    stop("'formula' must be a formula object.")
  }
  if (!is.data.frame(data)) {
    stop("'data' must be a data frame.")
  }
  if (!id %in% names(data)) {
    stop("Variable '", id, "' not found in data.")
  }
  if (!time %in% names(data)) {
    stop("Variable '", time, "' not found in data.")
  }

  reps <- as.integer(reps)
  if (reps < 1L) {
    stop("'reps' must be a positive integer.")
  }

  # Sort data by id and time
  data <- data[order(data[[id]], data[[time]]), , drop = FALSE]

  # Extract model matrices
  mf <- stats::model.frame(formula, data = data)
  Y <- stats::model.response(mf)
  X <- stats::model.matrix(formula, data = mf)

 # Remove intercept from X (will be partialled out if requested)
  intercept_col <- which(colnames(X) == "(Intercept)")
  if (length(intercept_col) > 0L) {
    X <- X[, -intercept_col, drop = FALSE]
  }

  if (ncol(X) < 1L) {
    stop("At least one regressor is required.")
  }

  # Extract id and time variables
  id_var <- data[[id]]
  time_var <- data[[time]]

  # Check for balanced panel
  units <- unique(id_var)
  N_g <- length(units)
  periods <- unique(time_var)
  T_g <- length(periods)

  if (nrow(data) != N_g * T_g) {
    stop("xtbhst requires a strongly balanced panel. ",
         "Expected ", N_g * T_g, " observations but found ", nrow(data), ".")
  }

  # Check that each unit has exactly T observations
  unit_counts <- table(id_var)
  if (any(unit_counts != T_g)) {
    stop("xtbhst requires a strongly balanced panel. ",
         "Not all units have the same number of observations.")
  }

  K <- ncol(X)

  # Build partial variables matrix (Z)
  Z <- NULL
  Kpartial <- 0L

  # Add constant to partial if requested
  if (constant) {
    Z <- matrix(1, nrow = nrow(X), ncol = 1L)
    colnames(Z) <- "constant"
    Kpartial <- 1L
  }

  # Add user-specified partial variables
  if (!is.null(partial)) {
    if (!inherits(partial, "formula")) {
      stop("'partial' must be a formula (e.g., ~ z1 + z2).")
    }
    partial_mf <- stats::model.frame(partial, data = data)
    partial_X <- stats::model.matrix(partial, data = data)
    # Remove intercept if present
    int_col <- which(colnames(partial_X) == "(Intercept)")
    if (length(int_col) > 0L) {
      partial_X <- partial_X[, -int_col, drop = FALSE]
    }
    if (ncol(partial_X) > 0L) {
      if (is.null(Z)) {
        Z <- partial_X
      } else {
        Z <- cbind(Z, partial_X)
      }
      Kpartial <- ncol(Z)
    }
  }

  # Add cross-sectional averages if requested
  if (!is.null(csa)) {
    if (!inherits(csa, "formula")) {
      stop("'csa' must be a formula (e.g., ~ x1 + x2).")
    }
    csa_mf <- stats::model.frame(csa, data = data)
    csa_X <- stats::model.matrix(csa, data = data)
    int_col <- which(colnames(csa_X) == "(Intercept)")
    if (length(int_col) > 0L) {
      csa_X <- csa_X[, -int_col, drop = FALSE]
    }

    if (ncol(csa_X) > 0L) {
      # Compute cross-sectional averages for each time period
      csa_vars <- .compute_csa(csa_X, id_var, time_var, csa_lags)
      if (is.null(Z)) {
        Z <- csa_vars
      } else {
        Z <- cbind(Z, csa_vars)
      }
      Kpartial <- ncol(Z)
    }
  }

  # Set block length if not provided
  if (is.null(blocklength)) {
    blocklength <- floor(2 * T_g^(1 / 3))
  }
  blocklength <- max(1L, as.integer(blocklength))

  # Run the bootstrap test
  result <- .xtbhst_bootstrap(
    Y = Y, X = X, Z = Z,
    id_var = id_var, time_var = time_var,
    N_g = N_g, T_g = T_g, K = K, Kpartial = Kpartial,
    reps = reps, blocklength = blocklength
 )

  # Build output
  out <- list(
    delta = result$delta,
    delta_adj = result$delta_adj,
    pval = result$pval,
    pval_adj = result$pval_adj,
    blocklength = blocklength,
    reps = reps,
    N = N_g,
    T = T_g,
    K = K,
    Kpartial = Kpartial,
    beta_i = result$beta_i,
    beta_fe = result$beta_fe,
    delta_stars = result$delta_stars,
    delta_adj_stars = result$delta_adj_stars,
    formula = formula,
    partial = partial,
    csa = csa,
    constant = constant,
    call = cl
  )

  class(out) <- "xtbhst"
  out
}


#' Compute cross-sectional averages
#'
#' @param X Matrix of variables.
#' @param id_var Vector of cross-sectional identifiers.
#' @param time_var Vector of time identifiers.
#' @param lags Number of lags to include.
#' @return Matrix of cross-sectional averages (with lags if requested).
#' @keywords internal
#' @noRd
.compute_csa <- function(X, id_var, time_var, lags = 0L) {
  times <- unique(time_var)
  n_times <- length(times)
  n_obs <- length(id_var)
  n_vars <- ncol(X)

  # Compute contemporaneous cross-sectional averages
  csa_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_vars)
  colnames(csa_mat) <- paste0("csa_", colnames(X))

  for (v in seq_len(n_vars)) {
    for (t_idx in seq_along(times)) {
      t_val <- times[t_idx]
      idx <- which(time_var == t_val)
      csa_mat[idx, v] <- mean(X[idx, v], na.rm = TRUE)
    }
  }

  # Add lags if requested
  if (lags > 0L) {
    units <- unique(id_var)
    lagged_list <- list(csa_mat)

    for (lag in seq_len(lags)) {
      lagged_mat <- matrix(NA_real_, nrow = n_obs, ncol = n_vars)
      colnames(lagged_mat) <- paste0("csa_", colnames(X), "_L", lag)

      for (u in units) {
        u_idx <- which(id_var == u)
        u_times <- time_var[u_idx]
        ord <- order(u_times)
        u_idx_sorted <- u_idx[ord]

        for (v in seq_len(n_vars)) {
          vals <- csa_mat[u_idx_sorted, v]
          if (length(vals) > lag) {
            lagged_mat[u_idx_sorted[(lag + 1):length(vals)], v] <-
              vals[1:(length(vals) - lag)]
          }
        }
      }
      lagged_list[[lag + 1L]] <- lagged_mat
    }

    csa_mat <- do.call(cbind, lagged_list)
  }

  csa_mat
}


#' Bootstrap core computation
#'
#' @param Y Dependent variable vector.
#' @param X Regressor matrix.
#' @param Z Matrix of variables to partial out (or NULL).
#' @param id_var Cross-sectional identifier vector.
#' @param time_var Time identifier vector.
#' @param N_g Number of cross-sectional units.
#' @param T_g Number of time periods.
#' @param K Number of regressors.
#' @param Kpartial Number of partialled-out variables.
#' @param reps Number of bootstrap replications.
#' @param blocklength Block length for block bootstrap.
#' @return List with test statistics and bootstrap results.
#' @keywords internal
#' @noRd
.xtbhst_bootstrap <- function(Y, X, Z, id_var, time_var,
                               N_g, T_g, K, Kpartial,
                               reps, blocklength) {
  units <- unique(id_var)

  # Create index mapping: for each unit, store row indices
  index <- vector("list", N_g)
  for (i in seq_along(units)) {
    index[[i]] <- which(id_var == units[i])
  }

  # Storage for precomputed matrices
  tmp_xx_list <- vector("list", N_g)
  tmp_xx1_list <- vector("list", N_g)
  tmp_zz1_list <- if (Kpartial > 0L) vector("list", N_g) else NULL

  # Partial out Z from Y and X if needed
  Y_orig <- Y
  X_orig <- X

  if (Kpartial > 0L && !is.null(Z)) {
    for (i in seq_len(N_g)) {
      idx <- index[[i]]
      Yi <- Y[idx]
      Xi <- X[idx, , drop = FALSE]
      Zi <- Z[idx, , drop = FALSE]

      # Remove rows with NA in Z
      complete <- stats::complete.cases(Zi)
      if (!all(complete)) {
        # For CSA with lags, we may have NAs at the beginning
        # This is a simplified handling - ideally would adjust T_g
        Zi[!complete, ] <- 0
      }

      tmp_zz <- crossprod(Zi)
      tmp_zz1 <- .safe_solve(tmp_zz)
      tmp_zz1_list[[i]] <- tmp_zz1

      # Partial out Z
      Y[idx] <- Yi - Zi %*% tmp_zz1 %*% crossprod(Zi, Yi)
      X[idx, ] <- Xi - Zi %*% tmp_zz1 %*% crossprod(Zi, Xi)
    }
  }

  # Compute pooled FE estimates and residuals
  XX_global <- crossprod(X)
  b_fe <- solve(XX_global, crossprod(X, Y))
  resid_fe <- Y - X %*% b_fe

  # Storage for unit-specific quantities
  E <- matrix(NA_real_, nrow = T_g, ncol = N_g)  # residuals
  sigma2 <- numeric(N_g)
  beta_i <- matrix(NA_real_, nrow = N_g, ncol = K)
  beta_wfe_up <- rep(0, K)
  beta_wfe_low <- matrix(0, K, K)

  for (i in seq_len(N_g)) {
    idx <- index[[i]]
    Yi <- Y[idx]
    Xi <- X[idx, , drop = FALSE]

    tmp_xx <- crossprod(Xi)
    tmp_xx1 <- .safe_solve(tmp_xx)
    tmp_xy <- crossprod(Xi, Yi)

    tmp_xx_list[[i]] <- tmp_xx
    tmp_xx1_list[[i]] <- tmp_xx1

    # Individual OLS estimates
    beta_i_vec <- tmp_xx1 %*% tmp_xy
    beta_i[i, ] <- as.vector(beta_i_vec)

    # Residuals
    resid_i <- Yi - Xi %*% beta_i_vec
    E[, i] <- resid_i

    # Variance estimate using FE residuals
    resid_fe_i <- resid_fe[idx]
    sigma2[i] <- sum(resid_fe_i^2) / (T_g - Kpartial)

    # Weighted FE aggregation
    beta_wfe_up <- beta_wfe_up + as.vector(tmp_xy) / sigma2[i]
    beta_wfe_low <- beta_wfe_low + tmp_xx / sigma2[i]
  }

  # Weighted FE estimate
  beta_wfe <- solve(beta_wfe_low, beta_wfe_up)

  # Compute S_tilde and Delta statistics
  S_tilde <- 0
  for (i in seq_len(N_g)) {
    beta_diff <- beta_i[i, ] - beta_wfe
    tmp_xx_scaled <- tmp_xx_list[[i]] / sigma2[i]
    S_tilde <- S_tilde + as.numeric(t(beta_diff) %*% tmp_xx_scaled %*% beta_diff)
  }

  delta_orig <- sqrt(N_g) * (S_tilde / N_g - K) / sqrt(2 * K)

  # Adjusted variance
  var_st <- 2 * K * (T_g - K - Kpartial - 1) / (T_g - Kpartial + 1)
  if (var_st <= 0) var_st <- 2 * K  # fallback
  delta_adj_orig <- sqrt(N_g) * (S_tilde / N_g - K) / sqrt(var_st)

  # Bootstrap
  delta_stars <- numeric(reps)
  delta_adj_stars <- numeric(reps)

  for (b in seq_len(reps)) {
    # Block bootstrap resampling of residuals
    E_star <- matrix(NA_real_, nrow = T_g, ncol = N_g)
    t_start <- 1L

    while (t_start <= T_g) {
      rand_idx <- sample.int(T_g - blocklength + 1L, 1L)
      len <- min(blocklength, T_g - t_start + 1L)
      E_star[t_start:(t_start + len - 1L), ] <-
        E[rand_idx:(rand_idx + len - 1L), ]
      t_start <- t_start + len
    }

    # Construct bootstrap Y* under H0 (homogeneous slopes)
    Y_star <- numeric(length(Y))

    for (i in seq_len(N_g)) {
      idx <- index[[i]]
      Xi <- X[idx, , drop = FALSE]
      Ei_star <- E_star[, i]

      # Partial out Z from bootstrap residuals if needed
      if (Kpartial > 0L && !is.null(Z)) {
        Zi <- Z[idx, , drop = FALSE]
        complete <- stats::complete.cases(Zi)
        if (!all(complete)) Zi[!complete, ] <- 0
        tmp_zz1 <- tmp_zz1_list[[i]]
        Ei_star <- Ei_star - Zi %*% tmp_zz1 %*% crossprod(Zi, Ei_star)
      }

      Y_star[idx] <- Xi %*% beta_wfe + Ei_star
    }

    # Re-estimate on bootstrap sample
    b_fe_star <- solve(XX_global, crossprod(X, Y_star))
    resid_fe_star <- Y_star - X %*% b_fe_star

    sigma2_star <- numeric(N_g)
    beta_star <- matrix(NA_real_, nrow = N_g, ncol = K)
    beta_wfe_up_star <- rep(0, K)
    beta_wfe_low_star <- matrix(0, K, K)

    for (i in seq_len(N_g)) {
      idx <- index[[i]]
      Y_star_i <- Y_star[idx]
      Xi <- X[idx, , drop = FALSE]

      tmp_xx <- tmp_xx_list[[i]]
      tmp_xx1 <- tmp_xx1_list[[i]]
      tmp_xy_star <- crossprod(Xi, Y_star_i)

      beta_i_star <- tmp_xx1 %*% tmp_xy_star
      beta_star[i, ] <- as.vector(beta_i_star)

      resid_fe_star_i <- resid_fe_star[idx]
      sigma2_star[i] <- sum(resid_fe_star_i^2) / (T_g - Kpartial)

      beta_wfe_up_star <- beta_wfe_up_star + as.vector(tmp_xy_star) / sigma2_star[i]
      beta_wfe_low_star <- beta_wfe_low_star + tmp_xx / sigma2_star[i]
    }

    beta_wfe_star <- solve(beta_wfe_low_star, beta_wfe_up_star)

    S_tilde_star <- 0
    for (i in seq_len(N_g)) {
      beta_diff <- beta_star[i, ] - beta_wfe_star
      tmp_xx_scaled <- tmp_xx_list[[i]] / sigma2_star[i]
      S_tilde_star <- S_tilde_star +
        as.numeric(t(beta_diff) %*% tmp_xx_scaled %*% beta_diff)
    }

    delta_stars[b] <- sqrt(N_g) * (S_tilde_star / N_g - K) / sqrt(2 * K)
    delta_adj_stars[b] <- sqrt(N_g) * (S_tilde_star / N_g - K) / sqrt(var_st)
  }

  # Bootstrap p-values (proportion of bootstrap stats exceeding observed)
  pval <- mean(delta_stars > delta_orig)
  pval_adj <- mean(delta_adj_stars > delta_adj_orig)

  list(
    delta = delta_orig,
    delta_adj = delta_adj_orig,
    pval = pval,
    pval_adj = pval_adj,
    beta_i = beta_i,
    beta_fe = beta_wfe,
    delta_stars = delta_stars,
    delta_adj_stars = delta_adj_stars
  )
}


#' Safe matrix inversion with regularization fallback
#'
#' @param x A square matrix to invert.
#' @param tol Tolerance for singularity detection.
#' @return The inverse of x, or a regularized inverse if singular.
#' @keywords internal
#' @noRd
.safe_solve <- function(x, tol = .Machine$double.eps^0.5) {
  tryCatch(
    solve(x),
    error = function(e) {
      # Use regularized inverse (ridge-like approach)
      n <- nrow(x)
      diag_add <- max(abs(diag(x))) * tol
      if (diag_add < tol) diag_add <- tol
      solve(x + diag_add * diag(n))
    }
  )
}
