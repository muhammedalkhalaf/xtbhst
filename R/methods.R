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
  cat("(Blomquist & Westerlund, 2015. Empirical Economics)\n")
  cat("H0: slope coefficients are homogeneous\n")
  cat(strrep("-", 45), "\n")

  # Results table
  cat(sprintf("  %12s   %12s\n", "Delta", "BS p-value"))
  cat(sprintf("  %12.*f   %12.*f\n", digits, x$delta, digits, x$pval))
  cat(sprintf("adj. %9.*f   %12.*f\n", digits, x$delta_adj, digits, x$pval_adj))
  cat(strrep("-", 45), "\n")

  cat("Bootstrap replications:", x$reps, "\n")
  cat("Block length:", x$blocklength, "\n")
  cat("Panel: N =", x$N, ", T =", x$T, ", K =", x$K, "\n")

  if (x$Kpartial > 0L) {
    cat("Variables partialled out:", x$Kpartial, "\n")
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
  cat("\n")

  cat("Test Results:\n")
  cat(strrep("-", 50), "\n")
  cat(sprintf("  %-20s %12s %12s\n", "Statistic", "Value", "BS p-value"))
  cat(strrep("-", 50), "\n")
  cat(sprintf("  %-20s %12.*f %12.*f\n",
              "Delta", digits, object$delta, digits, object$pval))
  cat(sprintf("  %-20s %12.*f %12.*f\n",
              "Delta (adjusted)", digits, object$delta_adj, digits, object$pval_adj))
  cat(strrep("-", 50), "\n\n")

  cat("Bootstrap settings:\n")
  cat("  Replications:", object$reps, "\n")
  cat("  Block length:", object$blocklength, "\n")
  cat("\n")

  # Interpretation
  alpha <- 0.05
  if (object$pval < alpha || object$pval_adj < alpha) {
    cat("Conclusion (at 5% level): Reject H0 - evidence of slope heterogeneity.\n")
  } else {
    cat("Conclusion (at 5% level): Fail to reject H0 - slopes appear homogeneous.\n")
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
  fe_df <- data.frame(Estimate = round(object$beta_fe, digits))
  rownames(fe_df) <- cnames
  print(fe_df)

  invisible(list(
    delta = object$delta,
    delta_adj = object$delta_adj,
    pval = object$pval,
    pval_adj = object$pval_adj,
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
#'   1 = Bootstrap distribution of Delta,
#'   2 = Bootstrap distribution of adjusted Delta,
#'   3+ = Individual coefficient distributions.
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
      x$delta_stars, x$delta,
      main = "Bootstrap Distribution: Delta",
      xlab = "Delta",
      observed_label = sprintf("Observed = %.3f", x$delta),
      ...
    )
  }

  if (2L %in% which) {
    .plot_bootstrap_dist(
      x$delta_adj_stars, x$delta_adj,
      main = "Bootstrap Distribution: Adjusted Delta",
      xlab = "Adjusted Delta",
      observed_label = sprintf("Observed = %.3f", x$delta_adj),
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
