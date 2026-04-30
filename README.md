# xtbhst: Bootstrap Slope Heterogeneity Test for Panel Data

<!-- badges: start -->
[![CRAN status](https://www.r-pkg.org/badges/version/xtbhst)](https://CRAN.R-project.org/package=xtbhst)
<!-- badges: end -->

## Overview

`xtbhst` implements the bootstrap slope heterogeneity test for panel data based on Blomquist and Westerlund (2015). The test examines whether slope coefficients are homogeneous across cross-sectional units.

**Reference:**
Blomquist, J., & Westerlund, J. (2015). Panel bootstrap tests of slope homogeneity. *Empirical Economics*, 48(1), 1191-1204. [doi:10.1007/s00181-015-0978-z](https://doi.org/10.1007/s00181-015-0978-z)

## Installation

```r
# Install from CRAN (when available)
install.packages("xtbhst")

# Or install development version from GitHub
# install.packages("devtools")
```

## Usage

### Basic Example

```r
library(xtbhst)

# Generate panel data with homogeneous slopes
set.seed(123)
N <- 20   # cross-sectional units
T <- 30   # time periods

data <- data.frame(
  id = rep(1:N, each = T),
  time = rep(1:T, N),
  x = rnorm(N * T)
)
data$y <- 1 + 0.5 * data$x + rnorm(N * T)

# Test for slope heterogeneity
result <- xtbhst(y ~ x, data = data, id = "id", time = "time",
                 reps = 999, seed = 42)
print(result)
```

Output:
```
Bootstrap test for slope heterogeneity
(Blomquist & Westerlund, 2015. Empirical Economics)
H0: slope coefficients are homogeneous
---------------------------------------------
         Delta      BS p-value
       -1.2345        0.8900
adj.   -1.3456        0.8700
---------------------------------------------
Bootstrap replications: 999 
Block length: 6 
Panel: N = 20 , T = 30 , K = 1 
```

### With Control Variables

```r
# Partial out control variables
data$z <- rnorm(N * T)
result <- xtbhst(y ~ x, data = data, id = "id", time = "time",
                 partial = ~ z, reps = 999)
```

### Handling Cross-Sectional Dependence

```r
# Include cross-sectional averages (Pesaran, 2006)
result <- xtbhst(y ~ x, data = data, id = "id", time = "time",
                 csa = ~ x, csa_lags = 2, reps = 999)
```

### Diagnostic Plots

```r
# Plot bootstrap distributions
plot(result)

# Plot all diagnostics including individual slopes
plot(result, which = 1:4)
```

## Interpretation

- **Null hypothesis (H0):** Slope coefficients are homogeneous across all cross-sectional units
- **Alternative (H1):** Slopes differ across units

If the p-value is small (e.g., < 0.05), reject H0 and conclude there is evidence of slope heterogeneity. This suggests that pooled OLS or standard fixed effects estimators may be inappropriate, and heterogeneous coefficient models (e.g., mean group estimator) should be considered.

## Requirements

- Strongly balanced panel (all units observed for all time periods)
- At least one regressor beyond the constant

## Citation

If you use this package, please cite:

```
Blomquist, J., & Westerlund, J. (2015). Panel bootstrap tests of slope 
homogeneity. Empirical Economics, 48(1), 1191-1204.
https://doi.org/10.1007/s00181-015-0978-z
```

## How to cite

If you use **xtbhst** in your work, please cite the package as:

> Alkhalaf, M.A. (2026). *xtbhst: Bootstrap Slope Heterogeneity Test for
> Panel Data*. R package version 1.0.2.
> <https://CRAN.R-project.org/package=xtbhst>.
> doi:10.32614/CRAN.package.xtbhst

In R, the suggested citation is also available via `citation("xtbhst")`.

The R port is based on the Stata `xthst` command; please also cite the
original Stata implementation when relevant.

## License

GPL (>= 3)

## Author

R port based on the Stata implementation `xtbhst`, which was adapted from the Stata `xthst` command.
