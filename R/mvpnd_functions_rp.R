#' The Multivariate Power Normal Distribution with Reparametrization
#'
#' Density and random generation for the multivariate (k variate) power normal
#' distribution (MVPND) with the reparametrization proposed by
#' Maruo et al. (2017).
#'
#' @param x Numeric vector (length: \code{k}) or matrix (\code{n} × \code{k})
#'   of quantiles.
#' @param n Numeric scalar for number of observations.
#' @param lambda Numeric vector of shape parameters.
#' @param xi Numeric vector of reparametrized location parameters (medians).
#' @param tau Numeric vector of reparametrized scale parameters
#'   (interquantile range divided by median).
#' @param corr Correlation parameter matrix.
#' @param log Logical; if \code{TRUE}, log-density is returned.
#' @param ... Additional arguments passed to \code{rmvnorm} (used in
#'   \code{rmvpnd_rp}).
#'
#' @details
#' These functions are wrapper functions for the MPND, using the
#' reparametrization method proposed by Maruo et al. (2017).
#' Specifically, they convert the reparametrized parameters \eqn{\xi_{0.5}}
#' (median) and \eqn{\tau = (\xi_{0.75} - \xi_{0.25}) / \xi_{0.5}} into
#' the original parameters \eqn{\mu} and \eqn{\sigma} for each margin via
#' \code{pnd_reparm}, and then call the underlying \code{dmvpnd} and
#' \code{rmvpnd} functions.
#'
#' Each marginal distribution follows a power normal distribution (PND) with
#' its own shape, location, and scale parameters. The dependence structure
#' between variables is modeled via the correlation matrix \code{corr}, as in
#' the multivariate normal distribution.
#'
#' @return
#' \code{dmvpnd_rp()} returns the probability density function (PDF),
#' and \code{rmvpnd_rp()} generates random samples from the multivariate
#' power normal distribution.
#'
#' @references
#' Maruo, K., Yamabe, T., & Yamaguchi, Y. (2017). Statistical simulation based on
#' right skewed distributions. \emph{Computational Statistics}, 32(3), 889–907.
#' \url{https://doi.org/10.1007/s00180-016-0664-4}
#'
#' @seealso \code{\link[powerNormal]{pnd_reparm}},
#'   \code{\link[powerNormal]{dmvpnd}}, \code{\link[powerNormal]{rmvpnd}}
#'
#' @examples
#' dmvpnd_rp(
#'   x = c(100, 100, 100),
#'   lambda = c(0.4, 0, -0.2),
#'   xi = c(100, 90, 80),
#'   tau = c(1, 0.9, 0.8),
#'   corr = diag(3)
#' )
#'
#' dmvpnd_rp(
#'   x = matrix(c(100, 100, 100, 98, 95, 90), 2, 3),
#'   lambda = c(0.4, 0, -0.2),
#'   xi = c(100, 90, 80),
#'   tau = c(1, 0.9, 0.8),
#'   corr = diag(3)
#' )
#'
#' rmvpnd_rp(
#'   n = 10,
#'   lambda = c(0.4, 0, -0.2),
#'   xi = c(100, 90, 80),
#'   tau = c(1, 0.9, 0.8),
#'   corr = diag(3)
#' )
#'
#' @name mvpnd_functions_rp
#' @rdname mvpnd_functions_rp
#' @export
dmvpnd_rp <- function(x, lambda, xi, tau, corr, log = FALSE){
  p <- length(lambda)
  if (class(x)[1] == "numeric") {
    if (length(x) != p) {
      stop("if y is numeric vector, dimension of y must equal to that of lambda.")
    } else {
      x <- t(x)
    }
  }
  mu <- numeric(p)
  sigma <- numeric(p)
  for (i in 1:p) {
    ms <- pnd_reparm(lambda[i], xi[i], tau[i])
    mu[i] <- ms$mu
    sigma[i] <- ms$sigma
  }
  fn <- dmvpnd(x, lambda, mu, sigma, corr, log)
  return(fn)
}

#' @rdname mvpnd_functions_rp
#' @export
rmvpnd_rp <- function(n, lambda, xi, tau, corr, ...){
  p <- length(lambda)
  mu <- numeric(p)
  sigma <- numeric(p)
  for (i in 1:p) {
    ms <- pnd_reparm(lambda[i], xi[i], tau[i])
    mu[i] <- ms$mu
    sigma[i] <- ms$sigma
  }
  x <- rmvpnd(n, lambda, mu, sigma, corr, ...)
  return(x)
}

