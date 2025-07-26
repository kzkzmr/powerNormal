# powerNormal

The **powerNormal** package provides functions to work with the power normal distribution (Goto & Inoue, 1980), including:

- Probability density, distribution function, quantiles, and random number generation
- Reparametrized versions using the median and interquantile range (Maruo et al., 2017)
- Multivariate extensions
- Parameter estimation based on maximum likelihood with normal or truncated normal assumptions (Maruo et al., 2011)

## Installation

You can install the development version from GitHub:

```r
# install.packages("devtools")
devtools::install_github("kzkzmr/powerNormal")
```

## Example

```r
library(powerNormal)

# Generate random values from a reparametrized PND
x <- rpnd_rp(100, lambda = 0.5, xi = 100, tau = 1.2)

# Estimate parameters
est <- pnd_est(x)

# Use multivariate PND
mu <- c(0, 0)
sigma <- c(1, 1)
lambda <- c(0.5, -0.2)
corr <- diag(2)
x_mv <- rmvpnd(100, lambda, mu, sigma, corr)
```

## References

- Goto M, Inoue T. Some properties of the power normal distribution. *Japanese Journal of Biometrics*, 1, 28–54, 1980. https://doi.org/10.5691/jjb.1.28
- Maruo K, Shirahata S, Goto M. Underlying assumptions of the power-normal distribution. *Behaviormetrika*, 38(1), 85–95, 2011. https://doi.org/10.2333/bhmk.38.85
- Maruo K, Yamabe T, Yamaguchi Y. Statistical simulation based on right-skewed distributions. *Computational Statistics*, 32(3), 889–907, 2017. https://doi.org/10.1007/s00180-016-0664-4

## Author

Kazushi Maruo
<https://github.com/kzkzmr>
---
