#' Print method for xtbhst objects
#'
#' @param x An object of class \code{"xtbhst"}.
#' @param digits Number of digits to display (default: 4).
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns the input object.
#' @export
#' @method print xtbhst
print.xtbhst <- function(x, digits = 4L, ...) {
  cat("\n")
  cat("Bootstrap test for slope heterogeneity\n")
  cat("(Blomquist and Westerlund, 2016, Empirical Economics)\n")
  cat("H0: slope coefficients are homogeneous\n")
  cat(strrep("-", 52), "\n")
  cat(sprintf("  %-16s %12s %14s\n", "Statistic", "Value", "p-value"))
  cat(sprintf("  %-16s %12.*f %14.*f  (bootstrap)\n",
              "S", digits, x$S, digits, x$pval))
  cat(sprintf("  %-16s %12.*f %14.*f  (asymptotic)\n",
              "Delta", digits, x$delta, digits, x$pval_delta_asy))
  cat(sprintf("  %-16s %12.*f %14.*f  (asymptotic)\n",
              "Delta (adjusted)", digits, x$delta_adj, digits,
              x$pval_delta_adj_asy))
  cat(strrep("-", 52), "\n")

  cat("Bootstrap replications:", x$reps, "\n")
  cat("Block length:", x$blocklength, "\n")
  cat("Variance estimator:",
      if (identical(x$variance, "py")) "Pesaran and Yamagata (2008)"
      else "Blomquist and Westerlund (2016)", "\n")
  cat("Panel: N =", x$N, ", T =", x$T, ", K =", x$K, "\n")

  if (x$Kpartial > 0L) {
    cat("Variables partialled out:", x$Kpartial, "\n")
  }
  if (!is.null(x$n_dropped) && x$n_dropped > 0L) {
    cat("Rows dropped for missing values:", x$n_dropped, "\n")
  }

  invisible(x)
}


#' Summary method for xtbhst objects
#'
#' @param object An object of class \code{"xtbhst"}.
#' @param digits Number of digits to display (default: 4).
#' @param ... Additional arguments (ignored).
#' @return Invisibly returns a list with summary statistics.
#' @export
#' @method summary xtbhst
summary.xtbhst <- function(object, digits = 4L, ...) {
  cat("\n")
  cat("============================================\n")
  cat("Bootstrap Slope Heterogeneity Test (xtbhst)\n")
  cat("============================================\n\n")

  cat("Call:\n")
  print(object$call)
  cat("\n")

  cat("Panel structure:\n")
  cat("  Cross-sectional units (N):", object$N, "\n")
  cat("  Time periods (T):", object$T, "\n")
  cat("  Regressors (K):", object$K, "\n")
  if (object$Kpartial > 0L) {
    cat("  Partialled variables:", object$Kpartial, "\n")
  }
  if (!is.null(object$n_dropped) && object$n_dropped > 0L) {
    cat("  Rows dropped for missing values:", object$n_dropped, "\n")
  }
  cat("\n")

  cat("Test Results:\n")
  cat(strrep("-", 62), "\n")
  cat(sprintf("  %-20s %12s %12s  %s\n", "Statistic", "Value", "p-value",
              "p-value type"))
  cat(strrep("-", 62), "\n")
  cat(sprintf("  %-20s %12.*f %12.*f  %s\n",
              "S", digits, object$S, digits, object$pval, "bootstrap"))
  cat(sprintf("  %-20s %12.*f %12.*f  %s\n",
              "Delta", digits, object$delta, digits, object$pval_delta_asy,
              "asymptotic"))
  cat(sprintf("  %-20s %12.*f %12.*f  %s\n",
              "Delta (adjusted)", digits, object$delta_adj, digits,
              object$pval_delta_adj_asy, "asymptotic"))
  cat(strrep("-", 62), "\n\n")

  cat("Bootstrap settings:\n")
  cat("  Replications:", object$reps, "\n")
  cat("  Block length:", object$blocklength, "\n")
  cat("  Variance estimator:",
      if (identical(object$variance, "py")) "Pesaran and Yamagata (2008)"
      else "Blomquist and Westerlund (2016)", "\n")
  cat("\n")

  # Interpretation, based on the bootstrap p-value of S
  alpha <- 0.05
  if (object$pval < alpha) {
    cat("Conclusion (bootstrap test at 5% level): Reject H0 - evidence of ",
        "slope heterogeneity.\n", sep = "")
  } else {
    cat("Conclusion (bootstrap test at 5% level): Fail to reject H0 - ",
        "slopes appear homogeneous.\n", sep = "")
  }
  cat("\n")

  # Summary of individual slopes
  cat("Individual slope coefficient summary:\n")
  beta_summary <- apply(object$beta_i, 2L, function(col) {
    c(Mean = mean(col), SD = stats::sd(col),
      Min = min(col), Max = max(col))
  })

  if (is.null(colnames(object$beta_i))) {
    cnames <- paste0("X", seq_len(ncol(object$beta_i)))
  } else {
    cnames <- colnames(object$beta_i)
  }
  colnames(beta_summary) <- cnames

  print(round(beta_summary, digits))
  cat("\n")

  cat("Weighted FE (pooled) estimates:\n")
  fe_df <- data.frame(Estimate = round(as.vector(object$beta_fe), digits))
  rownames(fe_df) <- cnames
  print(fe_df)

  invisible(list(
    S = object$S,
    pval = object$pval,
    delta = object$delta,
    delta_adj = object$delta_adj,
    pval_delta_asy = object$pval_delta_asy,
    pval_delta_adj_asy = object$pval_delta_adj_asy,
    beta_summary = beta_summary,
    beta_fe = object$beta_fe
  ))
}


#' Plot method for xtbhst objects
#'
#' Produces diagnostic plots for the bootstrap slope heterogeneity test.
#'
#' @param x An object of class \code{"xtbhst"}.
#' @param which Integer vector specifying which plots to produce:
#'   1 = Bootstrap distribution of the statistic S,
#'   2 = Bootstrap distribution of S on the Delta scale (the fixed
#'   rescaling \code{sqrt(N) (S/N - K)/sqrt(2K)} of S and S*),
#'   3+ = Individual coefficient distributions (plot 3 is the first
#'   regressor, plot 4 the second, and so on).
#'   Default is \code{c(1, 2)}.
#' @param ask Logical. If \code{TRUE}, prompt before each plot (default: \code{TRUE}
#'   if multiple plots and interactive session).
#' @param ... Additional arguments passed to plotting functions.
#' @return Invisibly returns \code{NULL}.
#' @export
#' @method plot xtbhst
plot.xtbhst <- function(x, which = c(1L, 2L), ask = NULL, ...) {
  if (is.null(ask)) {
    ask <- length(which) > 1L && grDevices::dev.interactive()
  }

  if (ask) {
    oask <- grDevices::devAskNewPage(TRUE)
    on.exit(grDevices::devAskNewPage(oask))
  }

  if (1L %in% which) {
    .plot_bootstrap_dist(
      x$S_stars, x$S,
      main = "Bootstrap Distribution: S",
      xlab = "S",
      observed_label = sprintf("Observed = %.3f", x$S),
      ...
    )
  }

  if (2L %in% which) {
    obs <- sqrt(x$N) * (x$S / x$N - x$K) / sqrt(2 * x$K)
    .plot_bootstrap_dist(
      x$delta_stars, obs,
      main = "Bootstrap Distribution: S on the Delta scale",
      xlab = "sqrt(N) (S/N - K) / sqrt(2K)",
      observed_label = sprintf("Observed = %.3f", obs),
      ...
    )
  }

  # Individual coefficient plots
  K <- x$K
  coef_plots <- which[which > 2L]
  for (k in coef_plots) {
    coef_idx <- k - 2L
    if (coef_idx <= K) {
      .plot_coef_dist(
        x$beta_i[, coef_idx], x$beta_fe[coef_idx],
        main = sprintf("Individual Slopes: Variable %d", coef_idx),
        xlab = "Estimate",
        ...
      )
    }
  }

  invisible(NULL)
}


#' Plot bootstrap distribution
#' @keywords internal
#' @noRd
.plot_bootstrap_dist <- function(bs_values, observed, main, xlab,
                                  observed_label, ...) {
  hist_out <- graphics::hist(bs_values, plot = FALSE)

  graphics::hist(
    bs_values,
    freq = FALSE,
    col = grDevices::rgb(0.53, 0.81, 0.98, 0.6),  # light blue
    border = "white",
    main = main,
    xlab = xlab,
    ylab = "Density",
    ...
  )

  # Add kernel density
  dens <- stats::density(bs_values)
  graphics::lines(dens, col = "navy", lwd = 2)

  # Add observed value line
  graphics::abline(v = observed, col = "darkred", lwd = 2, lty = 2)

  # Add legend
  graphics::legend(
    "topright",
    legend = c("Bootstrap density", observed_label),
    col = c("navy", "darkred"),
    lwd = 2,
    lty = c(1, 2),
    bty = "n",
    cex = 0.8
  )
}


#' Plot coefficient distribution
#' @keywords internal
#' @noRd
.plot_coef_dist <- function(coef_values, fe_value, main, xlab, ...) {
  graphics::hist(
    coef_values,
    freq = FALSE,
    col = grDevices::rgb(0, 0.5, 0.5, 0.4),  # teal
    border = "white",
    main = main,
    xlab = xlab,
    ylab = "Density",
    ...
  )

  dens <- stats::density(coef_values)
  graphics::lines(dens, col = "teal", lwd = 2)

  graphics::abline(v = fe_value, col = "darkred", lwd = 2, lty = 2)

  graphics::legend(
    "topright",
    legend = c("Individual slopes", sprintf("Pooled = %.3f", fe_value)),
    col = c("teal", "darkred"),
    lwd = 2,
    lty = c(1, 2),
    bty = "n",
    cex = 0.8
  )
}
