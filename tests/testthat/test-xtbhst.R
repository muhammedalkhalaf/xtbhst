make_panel <- function(N, T_periods, seed, K = 1L) {
  set.seed(seed)
  d <- data.frame(
    id = rep(seq_len(N), each = T_periods),
    time = rep(seq_len(T_periods), N)
  )
  d$x1 <- rnorm(N * T_periods)
  if (K >= 2L) d$x2 <- rnorm(N * T_periods)
  theta <- rnorm(N)
  d$y <- theta[d$id] + d$x1 + (if (K >= 2L) 0.5 * d$x2 else 0) +
    rnorm(N * T_periods) * rep(runif(N, 0.5, 2), each = T_periods)
  d
}

# Hand computation of S, Delta and Delta_adj from lm() fits
hand_stats <- function(d, variance, xvars = c("x1", "x2")) {
  units <- unique(d$id)
  N <- length(units)
  Tt <- length(unique(d$time))
  K <- length(xvars)
  f <- stats::reformulate(xvars, response = "y")
  dm <- function(v) v - stats::ave(v, d$id)
  X <- sapply(xvars, function(v) dm(d[[v]]))
  Y <- dm(d$y)
  fe <- stats::lm(update(f, . ~ . + factor(id)), data = d)
  bi <- matrix(unlist(lapply(units, function(i) stats::coef(stats::lm(f, d[d$id == i, ]))[-1])),
               ncol = K, byrow = TRUE)
  if (variance == "bw") {
    s2 <- sapply(units, function(i) sum(stats::resid(stats::lm(f, d[d$id == i, ]))^2) / Tt)
  } else {
    s2 <- sapply(units, function(i) {
      e <- stats::resid(fe)[d$id == i]
      sum(e^2) / (Tt - K - 1)
    })
  }
  A <- matrix(0, K, K)
  bb <- numeric(K)
  for (j in seq_along(units)) {
    Xi <- X[d$id == units[j], , drop = FALSE]
    A <- A + crossprod(Xi) / s2[j]
    bb <- bb + crossprod(Xi, Y[d$id == units[j]]) / s2[j]
  }
  bw <- as.vector(solve(A, bb))
  S <- 0
  for (j in seq_along(units)) {
    Xi <- X[d$id == units[j], , drop = FALSE]
    dd <- bi[j, ] - bw
    S <- S + as.numeric(t(dd) %*% crossprod(Xi) %*% dd) / s2[j]
  }
  list(S = S, beta_i = bi, beta_wfe = bw, sigma2 = s2,
       delta = sqrt(N) * (S / N - K) / sqrt(2 * K),
       delta_adj = sqrt(N) * (S / N - K) / sqrt(2 * K * (Tt - K - 1) / (Tt + 1)))
}

test_that("xtbhst runs on balanced panel and returns the documented elements", {
  d <- make_panel(10, 20, seed = 123)
  result <- xtbhst(y ~ x1, data = d, id = "id", time = "time",
                   reps = 99, seed = 42)

  expect_s3_class(result, "xtbhst")
  expect_true(is.numeric(result$S))
  expect_true(is.numeric(result$delta))
  expect_true(is.numeric(result$delta_adj))
  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_true(result$pval_delta_asy >= 0 && result$pval_delta_asy <= 1)
  expect_equal(result$pval, mean(result$S_stars >= result$S))
  expect_length(result$S_stars, 99)
  expect_equal(result$N, 10)
  expect_equal(result$T, 20)
  expect_equal(result$K, 1)
  expect_equal(result$reps, 99)
  expect_equal(result$variance, "bw")
  expect_equal(result$n_dropped, 0)
  expect_null(result$pval_adj)
})

test_that("S, Delta and Delta_adj match hand computation, variance = 'bw'", {
  d <- make_panel(6, 20, seed = 1, K = 2L)
  r <- xtbhst(y ~ x1 + x2, data = d, id = "id", time = "time",
              reps = 19, seed = 3, variance = "bw")
  h <- hand_stats(d, "bw")
  expect_equal(r$S, h$S, tolerance = 1e-10)
  expect_equal(unname(r$beta_i), unname(h$beta_i), tolerance = 1e-10)
  expect_equal(unname(r$beta_fe), h$beta_wfe, tolerance = 1e-10)
  expect_equal(r$sigma2, h$sigma2, tolerance = 1e-10)
  # Delta and Delta_adj are always the Pesaran and Yamagata statistics
  h_py <- hand_stats(d, "py")
  expect_equal(r$S_py, h_py$S, tolerance = 1e-10)
  expect_equal(r$delta, h_py$delta, tolerance = 1e-10)
  expect_equal(r$delta_adj, h_py$delta_adj, tolerance = 1e-10)
  expect_false(isTRUE(all.equal(r$delta, sqrt(6) * (h$S / 6 - 2) / sqrt(4))))
  expect_equal(r$pval_delta_asy, stats::pnorm(h_py$delta, lower.tail = FALSE))
  expect_equal(r$pval_delta_adj_asy, stats::pnorm(h_py$delta_adj, lower.tail = FALSE))
  # Delta_adj uses the Pesaran and Yamagata (2008) variance 2K(T-K-1)/(T+1)
  expect_equal(r$delta_adj / r$delta,
               sqrt(2 * 2) / sqrt(2 * 2 * (20 - 2 - 1) / (20 + 1)),
               tolerance = 1e-12)
  # bootstrap Delta* is the rescaled S*
  expect_equal(r$delta_stars, sqrt(6) * (r$S_stars / 6 - 2) / sqrt(4))
})

test_that("S, Delta and Delta_adj match hand computation, variance = 'py'", {
  d <- make_panel(6, 20, seed = 1, K = 2L)
  r <- xtbhst(y ~ x1 + x2, data = d, id = "id", time = "time",
              reps = 19, seed = 3, variance = "py")
  h <- hand_stats(d, "py")
  expect_equal(r$variance, "py")
  expect_equal(r$S, h$S, tolerance = 1e-10)
  expect_equal(r$S_py, r$S)
  expect_equal(r$delta, h$delta, tolerance = 1e-10)
  expect_equal(r$delta_adj, h$delta_adj, tolerance = 1e-10)
  expect_equal(unname(r$beta_fe), h$beta_wfe, tolerance = 1e-10)
  expect_equal(r$sigma2, h$sigma2, tolerance = 1e-10)
  # The two variance options give different S
  r_bw <- xtbhst(y ~ x1 + x2, data = d, id = "id", time = "time",
                 reps = 19, seed = 3, variance = "bw")
  expect_false(isTRUE(all.equal(r$S, r_bw$S)))
})

test_that("bootstrap statistics use the chosen variance estimator", {
  # With reps = 1 and a fixed seed, recompute S* by hand from the
  # resampled residuals for both variance options.
  d <- make_panel(4, 12, seed = 7, K = 1L)
  for (v in c("bw", "py")) {
    r <- xtbhst(y ~ x1, data = d, id = "id", time = "time",
                reps = 1, seed = 11, blocklength = 3, variance = v)
    N <- 4; Tt <- 12
    # unit LS residuals and theta_i, then the same block draw
    fits <- lapply(1:N, function(i) stats::lm(y ~ x1, d[d$id == i, ]))
    E <- sapply(fits, stats::resid)
    set.seed(11)
    E_star <- matrix(NA_real_, Tt, N)
    t_start <- 1
    while (t_start <= Tt) {
      ri <- sample.int(Tt - 3 + 1, 1)
      len <- min(3, Tt - t_start + 1)
      E_star[t_start:(t_start + len - 1), ] <- E[ri:(ri + len - 1), ]
      t_start <- t_start + len
    }
    d_star <- d
    for (i in 1:N) {
      idx <- d$id == i
      theta_i <- mean(d$y[idx]) - r$beta_fe * mean(d$x1[idx])
      d_star$y[idx] <- theta_i + r$beta_fe * d$x1[idx] + E_star[, i]
    }
    h <- hand_stats(d_star, v, xvars = "x1")
    expect_equal(r$S_stars, h$S, tolerance = 1e-10)
  }
})

test_that("Delta_adj is NA with a warning when T - K - 1 <= 0", {
  set.seed(5)
  N <- 6; Tt <- 3
  d <- data.frame(id = rep(1:N, each = Tt), time = rep(1:Tt, N),
                  x1 = rnorm(N * Tt))
  d$y <- rnorm(N * Tt)
  # T - K - 1 = 3 - 1 - 1 = 1 > 0 with a constant, so use constant = FALSE
  # and T = 2 so that T - K - 1 = 0 while T > K + Kpartial
  d2 <- d[d$time <= 2, ]
  expect_warning(
    r <- xtbhst(y ~ x1, data = d2, id = "id", time = "time", reps = 5,
                constant = FALSE, blocklength = 1, seed = 1),
    "adjusted Delta"
  )
  expect_true(is.na(r$delta_adj))
  expect_true(is.na(r$pval_delta_adj_asy))
  expect_false(is.na(r$delta))
})

test_that("xtbhst detects heterogeneous slopes", {
  set.seed(456)
  N <- 15
  T_periods <- 30
  d <- data.frame(id = rep(1:N, each = T_periods),
                  time = rep(1:T_periods, N))
  d$x <- rnorm(N * T_periods)
  slopes <- rnorm(N, mean = 0.5, sd = 1)
  d$y <- 1 + slopes[d$id] * d$x + rnorm(N * T_periods, sd = 0.5)

  result <- xtbhst(y ~ x, data = d, id = "id", time = "time",
                   reps = 99, seed = 123)
  expect_true(result$delta > 0)
  expect_true(result$pval < 0.05)
})

test_that("xtbhst fails on unbalanced panel and on duplicated (id, time)", {
  d <- make_panel(10, 20, seed = 789)
  expect_error(
    xtbhst(y ~ x1, data = d[-c(1, 2, 3), ], id = "id", time = "time", reps = 9),
    "strongly balanced"
  )
  d_dup <- rbind(d, d[1, ])
  expect_error(
    xtbhst(y ~ x1, data = d_dup, id = "id", time = "time", reps = 9),
    "strongly balanced"
  )
})

test_that("missing values are dropped before the balance check", {
  d <- make_panel(6, 20, seed = 1, K = 2L)
  # NA in one row only: the panel becomes unbalanced, clear error
  d_na <- d
  d_na$y[3] <- NA
  expect_error(
    xtbhst(y ~ x1 + x2, data = d_na, id = "id", time = "time", reps = 9),
    "strongly balanced"
  )
  # NA in period 1 for every unit: the remaining panel is balanced with T = 19
  d_na2 <- d
  d_na2$x1[d_na2$time == 1] <- NA
  r <- xtbhst(y ~ x1 + x2, data = d_na2, id = "id", time = "time",
              reps = 9, seed = 1)
  expect_equal(r$T, 19)
  expect_equal(r$n_dropped, 6)
  h <- hand_stats(d[d$time > 1, ], "bw")
  expect_equal(r$S, h$S, tolerance = 1e-10)
  expect_equal(unname(r$beta_i), unname(h$beta_i), tolerance = 1e-10)
  # NA in a variable not used by the model is ignored
  d_na3 <- d
  d_na3$x2[5] <- NA
  r3 <- xtbhst(y ~ x1, data = d_na3, id = "id", time = "time", reps = 9, seed = 1)
  expect_equal(r3$n_dropped, 0)
  expect_equal(r3$T, 20)
})

test_that("xtbhst handles multiple regressors", {
  d <- make_panel(10, 25, seed = 111, K = 2L)
  result <- xtbhst(y ~ x1 + x2, data = d, id = "id", time = "time",
                   reps = 19, seed = 42)
  expect_equal(result$K, 2)
  expect_equal(dim(result$beta_i), c(10, 2))
  expect_equal(colnames(result$beta_i), c("x1", "x2"))
  expect_equal(names(result$beta_fe), c("x1", "x2"))
})

test_that("partial option matches lm() with the control variable", {
  d <- make_panel(8, 25, seed = 222, K = 1L)
  d$z <- rnorm(nrow(d))
  d$y <- d$y + 0.2 * d$z
  result <- xtbhst(y ~ x1, data = d, id = "id", time = "time",
                   partial = ~ z, reps = 9, seed = 42)
  expect_equal(result$K, 1)
  expect_equal(result$Kpartial, 2)
  bi <- sapply(1:8, function(i) stats::coef(stats::lm(y ~ x1 + z, d[d$id == i, ]))["x1"])
  expect_equal(unname(result$beta_i[, 1]), unname(bi), tolerance = 1e-10)
  s2 <- sapply(1:8, function(i) sum(stats::resid(stats::lm(y ~ x1 + z, d[d$id == i, ]))^2) / 25)
  expect_equal(result$sigma2, s2, tolerance = 1e-10)
})

test_that("csa_lags drops the first periods and matches lm() on the common sample", {
  set.seed(2)
  N <- 5; Tt <- 15
  d <- data.frame(id = rep(1:N, each = Tt), time = rep(1:Tt, N),
                  x = rnorm(N * Tt))
  d$y <- rnorm(N)[d$id] + d$x + rnorm(N * Tt)
  r <- xtbhst(y ~ x, data = d, id = "id", time = "time", reps = 9,
              csa = ~ y, csa_lags = 1, seed = 1)
  expect_equal(r$T, Tt - 1)
  expect_equal(r$Kpartial, 3)  # constant, csa_y, csa_y_L1
  d$c0 <- stats::ave(d$y, d$time)
  d$c1 <- unlist(tapply(d$c0, d$id, function(v) c(NA, v[-length(v)])))
  bi <- sapply(1:N, function(i) stats::coef(stats::lm(y ~ x + c0 + c1, d[d$id == i, ]))["x"])
  expect_equal(unname(r$beta_i[, 1]), unname(bi), tolerance = 1e-10)
  s2 <- sapply(1:N, function(i) sum(stats::resid(stats::lm(y ~ x + c0 + c1, d[d$id == i, ]))^2) / (Tt - 1))
  expect_equal(r$sigma2, s2, tolerance = 1e-10)
  # csa of the regressor only: the same sample, but Z fixed in the bootstrap
  r2 <- xtbhst(y ~ x, data = d, id = "id", time = "time", reps = 9,
               csa = ~ x, csa_lags = 2, seed = 1)
  expect_equal(r2$T, Tt - 2)
  expect_equal(r2$Kpartial, 4)
  expect_error(
    xtbhst(y ~ x, data = d, id = "id", time = "time", reps = 9, csa_lags = 1),
    "requires 'csa'"
  )
  expect_error(
    xtbhst(log(abs(y)) ~ x, data = d, id = "id", time = "time", reps = 9,
           csa = ~ y, csa_lags = 1),
    "plain column"
  )
})

test_that("print and summary methods work", {
  d <- make_panel(10, 20, seed = 333)
  result <- xtbhst(y ~ x1, data = d, id = "id", time = "time",
                   reps = 19, seed = 42)
  output <- capture.output(print(result))
  expect_true(any(grepl("Bootstrap test", output)))
  expect_true(any(grepl("bootstrap\\)", output)))
  expect_true(any(grepl("asymptotic\\)", output)))
  expect_true(any(grepl("2016", output)))

  summary_output <- capture.output(s <- summary(result))
  expect_true(any(grepl("Panel structure", summary_output)))
  expect_true(any(grepl("Conclusion", summary_output)))
  expect_equal(s$S, result$S)
})

test_that("blocklength defaults to round(2 T^(1/3)) and is validated", {
  d <- make_panel(6, 25, seed = 444)
  result <- xtbhst(y ~ x1, data = d, id = "id", time = "time",
                   reps = 9, seed = 42)
  expect_equal(result$blocklength, 6L)  # round(2 * 25^(1/3)) = 6, not floor = 5
  expect_equal(round(2 * c(25, 50, 100)^(1 / 3)), c(6, 7, 9))

  result2 <- xtbhst(y ~ x1, data = d, id = "id", time = "time",
                    reps = 9, blocklength = 10, seed = 42)
  expect_equal(result2$blocklength, 10L)
  expect_error(
    xtbhst(y ~ x1, data = d, id = "id", time = "time", reps = 9, blocklength = 26),
    "must not exceed"
  )
  expect_error(
    xtbhst(y ~ x1, data = d, id = "id", time = "time", reps = 9, blocklength = 0),
    "positive integer"
  )
})

test_that("singular unit regressions give an informative error", {
  d <- make_panel(6, 20, seed = 1, K = 2L)
  d$x2[d$id == 3] <- 1
  expect_error(
    xtbhst(y ~ x1 + x2, data = d, id = "id", time = "time", reps = 9),
    "singular for unit '3'"
  )
})

test_that("seed produces reproducible results", {
  d <- make_panel(10, 20, seed = 555)
  result1 <- xtbhst(y ~ x1, data = d, id = "id", time = "time",
                    reps = 29, seed = 12345)
  result2 <- xtbhst(y ~ x1, data = d, id = "id", time = "time",
                    reps = 29, seed = 12345)
  expect_equal(result1$S, result2$S)
  expect_equal(result1$pval, result2$pval)
  expect_equal(result1$S_stars, result2$S_stars)
})
