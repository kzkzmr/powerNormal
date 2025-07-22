#' The Multivariate Power Normal Distribution
#'
#' Density and random generation for the multivariate (k-variate)
#' power normal distribution (MVPND).
#'
#' @param x Numeric vector (length k) or matrix (n × k) of quantiles.
#' @param n Numeric scalar: number of observations.
#' @param lambda Numeric vector of shape parameters (length k).
#' @param mu Numeric vector of location parameters (length k).
#' @param sigma Numeric vector of scale parameters (length k).
#' @param corr Correlation parameter matrix (k × k).
#' @param log Logical; if \code{TRUE}, log-density is returned.
#' @param ... Additional arguments passed to \code{rmvnorm} (used in
#'   \code{rmvpnd}).
#'
#' @details
#' Each marginal distribution in the MVPND follows the univariate power normal
#' distribution, as defined in the \code{pnd_functions}.
#' The dependence structure across variables is parameterized via a
#' correlation parameter matrix \code{corr}, analogous to the multivariate
#' (truncated) normal distribution.
#'
#' @return
#' \code{dmvpnd()} returns the joint density of the MVPND.
#' \code{rmvpnd()} generates random variates from the distribution.
#'
#' @seealso
#' \code{\link[powerNormal]{dpnd}},
#' \code{\link[powerNormal]{ppnd}},
#' \code{\link[powerNormal]{qpnd}},
#' \code{\link[powerNormal]{rpnd}}
#'
#' @examples
#' dmvpnd(x = c(1, 1, 1), lambda = c(0.4, 0, -0.2), mu = c(1, 1, 1),
#'        sigma = c(1, 1, 1), corr = diag(3))
#' dmvpnd(x = matrix(c(1, 1, 1, 1, 1, 1), 2, 3), lambda = c(0.4, 0, -0.2),
#'        mu = c(1, 1, 1), sigma = c(1, 1, 1), corr = diag(3))
#' rmvpnd(n = 10, lambda = c(0.4, 0, -0.2), mu = c(1, 1, 1),
#'        sigma = c(1, 1, 1), corr = diag(3))
#'
#' @importFrom bcmixed bct
#' @importFrom mvtnorm dmvnorm pmvnorm rmvnorm
#'
#' @name mvpnd_functions
#' @rdname mvpnd_functions
#' @export
dmvpnd <- function(x, lambda, mu, sigma, corr, log = FALSE){
  k <- length(lambda)
  if (class(x)[1] == "numeric") {
    if (length(x) != k) {
      stop("if y is numeric vector, dimension of y must equal to that of lambda.")
    } else {
      x <- t(x)
    }
  }
  sK <- numeric(k) + Inf
  for (i in 1:k) {
    if (lambda[i] != 0) {
      sK[i] <- sign(lambda[i]) * (1 + lambda[i] * mu[i]) / (lambda[i] * sigma[i])
    }
  }
  AK <- pmvnorm(upper = sK, corr = corr)

  z <- x
  zm1 <- x
  for (i in 1:k){
    z[, i] <- bct(x[, i], lambda[i])
    zm1[, i] <- x[, i] ^ (lambda[i] - 1)
  }
  Sigma <- sigma %*% t(sigma) * corr
  if (log) {
    fn <- apply(log(zm1), 1, sum) - log (AK) + dmvnorm(z, mu, Sigma, log = TRUE)
  } else {
    fn <- apply(zm1, 1, prod) / AK * dmvnorm(z, mu, Sigma)
  }
  return(as.numeric(fn))
}

#' @rdname mvpnd_functions
#' @export
rmvpnd <- function(n, lambda, mu, sigma, corr, ...){
  k <- length(lambda)
  sK <- numeric(k) + Inf
  for (i in 1:k) {
    if (lambda[i] != 0) {
      sK[i] <- sign(lambda[i]) * (1 + lambda[i] * mu[i]) / (lambda[i] * sigma[i])
    }
  }
  AK <- pmvnorm(upper = sK, corr = corr)
  n2 <- round(n / AK * 2)

  z0 <- rmvnorm(n2, sigma = corr, )
  z1 <- z0

  cut0 <- matrix(TRUE, n2, k)
  for (i in 1:k){
    z1[, i] <- z0[, i] * sigma[i] + mu[i]
    if (lambda[i] < 0) {
      cut0[, i] <- (z1[, i] < -1 / lambda[i])
    }
    if (lambda[i] > 0) {
      cut0[, i] <- (z1[, i] > -1 / lambda[i])
    }
  }
  z <- z1[apply(cut0, 1, prod) == 1, ][1:n, ]
  y <- z
  for (i in 1:k){
    if (lambda[i] != 0){
      y[, i] <- (lambda[i] * z[, i] + 1) ^ (1 / lambda[i])
    } else {
      y[, i] <- exp(z[, i])
    }
  }
  return(y)
}

