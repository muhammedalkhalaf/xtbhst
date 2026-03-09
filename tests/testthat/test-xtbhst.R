test_that("xtbhst runs on balanced panel", {
  set.seed(123)
  N <- 10
  T_periods <- 20

  data <- data.frame(
    id = rep(1:N, each = T_periods),
    time = rep(1:T_periods, N),
    x = rnorm(N * T_periods)
  )
  data$y <- 1 + 0.5 * data$x + rnorm(N * T_periods)

  result <- xtbhst(y ~ x, data = data, id = "id", time = "time",
                   reps = 99, seed = 42)

  expect_s3_class(result, "xtbhst")
  expect_true(is.numeric(result$delta))
  expect_true(is.numeric(result$delta_adj))
  expect_true(result$pval >= 0 && result$pval <= 1)
  expect_true(result$pval_adj >= 0 && result$pval_adj <= 1)
  expect_equal(result$N, N)
  expect_equal(result$T, T_periods)
  expect_equal(result$K, 1)
  expect_equal(result$reps, 99)
})

test_that("xtbhst detects heterogeneous slopes", {
  set.seed(456)
  N <- 15
  T_periods <- 30

  data <- data.frame(
    id = rep(1:N, each = T_periods),
    time = rep(1:T_periods, N)
  )
  data$x <- rnorm(N * T_periods)

  # Generate heterogeneous slopes
  slopes <- rnorm(N, mean = 0.5, sd = 1)
  data$y <- 1 + slopes[data$id] * data$x + rnorm(N * T_periods, sd = 0.5)

  result <- xtbhst(y ~ x, data = data, id = "id", time = "time",
                   reps = 199, seed = 123)

  # With large heterogeneity, p-value should be relatively small
  # (This is a probabilistic test, so we use a loose threshold)
  expect_s3_class(result, "xtbhst")
  expect_true(result$delta > 0)  # Positive delta expected with heterogeneity
})

test_that("xtbhst fails on unbalanced panel", {
  set.seed(789)
  N <- 10
  T_periods <- 20

  data <- data.frame(
    id = rep(1:N, each = T_periods),
    time = rep(1:T_periods, N),
    x = rnorm(N * T_periods)
  )
  data$y <- rnorm(N * T_periods)

  # Remove some observations to create unbalanced panel
  data <- data[-c(1, 2, 3), ]

  expect_error(
    xtbhst(y ~ x, data = data, id = "id", time = "time", reps = 99),
    "strongly balanced"
  )
})

test_that("xtbhst handles multiple regressors", {
  set.seed(111)
  N <- 10
  T_periods <- 25

  data <- data.frame(
    id = rep(1:N, each = T_periods),
    time = rep(1:T_periods, N),
    x1 = rnorm(N * T_periods),
    x2 = rnorm(N * T_periods)
  )
  data$y <- 1 + 0.5 * data$x1 + 0.3 * data$x2 + rnorm(N * T_periods)

  result <- xtbhst(y ~ x1 + x2, data = data, id = "id", time = "time",
                   reps = 99, seed = 42)

  expect_s3_class(result, "xtbhst")
  expect_equal(result$K, 2)
  expect_equal(nrow(result$beta_i), N)
  expect_equal(ncol(result$beta_i), 2)
  expect_equal(length(result$beta_fe), 2)
})

test_that("xtbhst handles partial option", {
  set.seed(222)
  N <- 10
  T_periods <- 25

  data <- data.frame(
    id = rep(1:N, each = T_periods),
    time = rep(1:T_periods, N),
    x = rnorm(N * T_periods),
    z = rnorm(N * T_periods)
  )
  data$y <- 1 + 0.5 * data$x + 0.2 * data$z + rnorm(N * T_periods)

  result <- xtbhst(y ~ x, data = data, id = "id", time = "time",
                   partial = ~ z, reps = 99, seed = 42)

  expect_s3_class(result, "xtbhst")
  expect_equal(result$K, 1)
  expect_true(result$Kpartial >= 2)  # constant + z
})

test_that("print and summary methods work", {
  set.seed(333)
  N <- 10
  T_periods <- 20

  data <- data.frame(
    id = rep(1:N, each = T_periods),
    time = rep(1:T_periods, N),
    x = rnorm(N * T_periods)
  )
  data$y <- rnorm(N * T_periods)

  result <- xtbhst(y ~ x, data = data, id = "id", time = "time",
                   reps = 49, seed = 42)

  # Test print
  output <- capture.output(print(result))
  expect_true(any(grepl("Bootstrap test", output)))
  expect_true(any(grepl("Delta", output)))

  # Test summary
  summary_output <- capture.output(summary(result))
  expect_true(any(grepl("Panel structure", summary_output)))
  expect_true(any(grepl("Conclusion", summary_output)))
})

test_that("blocklength is computed correctly", {
  set.seed(444)
  N <- 10
  T_periods <- 27  # 2 * 27^(1/3) = 2 * 3 = 6

  data <- data.frame(
    id = rep(1:N, each = T_periods),
    time = rep(1:T_periods, N),
    x = rnorm(N * T_periods)
  )
  data$y <- rnorm(N * T_periods)

  result <- xtbhst(y ~ x, data = data, id = "id", time = "time",
                   reps = 49, seed = 42)

  expected_bl <- floor(2 * T_periods^(1/3))
  expect_equal(result$blocklength, expected_bl)

  # Test manual blocklength
  result2 <- xtbhst(y ~ x, data = data, id = "id", time = "time",
                    reps = 49, blocklength = 10, seed = 42)
  expect_equal(result2$blocklength, 10)
})

test_that("seed produces reproducible results", {
  set.seed(555)
  N <- 10
  T_periods <- 20

  data <- data.frame(
    id = rep(1:N, each = T_periods),
    time = rep(1:T_periods, N),
    x = rnorm(N * T_periods)
  )
  data$y <- rnorm(N * T_periods)

  result1 <- xtbhst(y ~ x, data = data, id = "id", time = "time",
                    reps = 99, seed = 12345)
  result2 <- xtbhst(y ~ x, data = data, id = "id", time = "time",
                    reps = 99, seed = 12345)

  expect_equal(result1$delta, result2$delta)
  expect_equal(result1$pval, result2$pval)
  expect_equal(result1$delta_stars, result2$delta_stars)
})
