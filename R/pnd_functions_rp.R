#' The Power Normal Distribution with Reparametrization
#'
#' Probability density function (PDF), cumulative distribution function (CDF),
#' quantile function, and random number generation for the power normal
#' distribution (PND) with the reparametrization proposed by Maruo et al.
#' (2017).
#'
#' @param x Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param lambda Numeric scalar: shape parameter.
#' @param xi Numeric scalar: reparametrized location parameter (i.e.,
#'   the median).
#' @param tau Numeric scalar: reparametrized scale parameter
#'   (i.e., interquantile range divided by median).
#' @param log Logical; if \code{TRUE}, log-density is returned.
#'
#' @details
#' These functions are wrapper functions that calculate the original
#' parameters \eqn{\mu} and \eqn{\sigma} from the reparametrized parameters
#' \eqn{\xi_{0.5}} and \eqn{\tau}, using the \code{pnd_reparm} function.
#' The resulting values are passed to the underlying functions:
#' \code{dpnd}, \code{ppnd}, \code{qpnd}, and \code{rpnd}.
#'
#' @return
#' \code{dpnd_rp()} returns the PDF,
#' \code{ppnd_rp()} returns the CDF,
#' \code{qpnd_rp()} returns the quantile function, and
#' \code{rpnd_rp()} generates random values from the distribution.
#'
#' @references
#' Goto, M., & Inoue, T. (1980). Some properties of the power normal distribution.
#' \emph{Japanese Journal of Biometrics}, 1, 28–54.
#' \url{https://doi.org/10.5691/jjb.1.28}
#'
#' Maruo, K., Yamabe, T., & Yamaguchi, Y. (2017). Statistical simulation based on
#' right-skewed distributions. \emph{Computational Statistics}, 32(3), 889–907.
#' \url{https://doi.org/10.1007/s00180-016-0664-4}
#'
#' @seealso \code{\link[powerNormal]{pnd_reparm}}
#'
#' @examples
#' dpnd_rp(x = 50, lambda = 0.2, xi = 100, tau = 0.8)
#' ppnd_rp(x = 50, lambda = 0.2, xi = 100, tau = 0.8)
#' qpnd_rp(p = 0.4, lambda = 0.2, xi = 100, tau = 0.8)
#' rpnd_rp(n = 10, lambda = 0.2, xi = 100, tau = 0.8)
#'
#' @name pnd_functions_rp
#' @rdname pnd_functions_rp
#' @export
dpnd_rp <- function(x, lambda, xi, tau, log = FALSE) {
  if (!is.numeric(lambda) || length(lambda) != 1 ||
      !is.numeric(xi) || length(xi) != 1 ||
      !is.numeric(tau) || length(tau) != 1) {
    stop("Arguments 'lambda', 'xi', and 'tau' must be scalar numeric values.")
  }
  if (!is.numeric(x) || any(x < 0)) {
    stop("Argument 'x' must be a non-negative numeric vector.")
  }
  ms <- pnd_reparm(lambda, xi, tau)
  fx <- dpnd(x, lambda, ms$mu, ms$sigma, log)
  return(fx)
}

#' @rdname pnd_functions_rp
#' @export
ppnd_rp <- function(x, lambda, xi, tau) {
  if (!is.numeric(lambda) || length(lambda) != 1 ||
      !is.numeric(xi) || length(xi) != 1 ||
      !is.numeric(tau) || length(tau) != 1) {
    stop("Arguments 'lambda', 'xi', and 'tau' must be scalar numeric values.")
  }
  if (!is.numeric(x) || any(x < 0)) {
    stop("Argument 'x' must be a non-negative numeric vector.")
  }
  ms <- pnd_reparm(lambda, xi, tau)
  Fx <- ppnd(x, lambda, ms$mu, ms$sigma)
  return(Fx)
}


#' @rdname pnd_functions_rp
#' @export
qpnd_rp <- function(p, lambda, xi, tau) {
  if (!is.numeric(lambda) || length(lambda) != 1 ||
      !is.numeric(xi) || length(xi) != 1 ||
      !is.numeric(tau) || length(tau) != 1) {
    stop("Arguments 'lambda', 'xi', and 'tau' must be scalar numeric values.")
  }
  if (!is.numeric(p)) {
    stop("Argument 'p' must be a numeric vector.")
  } else if (any(p < 0 | p > 1)) {
    stop("All elements of 'p' must be between 0 and 1.")
  }
  ms <- pnd_reparm(lambda, xi, tau)
  q <- qpnd(p, lambda, ms$mu, ms$sigma)
  return(q)
}


#' @rdname pnd_functions_rp
#' @export
rpnd_rp <- function(n, lambda, xi, tau) {
  if (!is.numeric(lambda) || length(lambda) != 1 ||
      !is.numeric(xi) || length(xi) != 1 ||
      !is.numeric(tau) || length(tau) != 1) {
    stop("Arguments 'lambda', 'xi', and 'tau' must be scalar numeric values.")
  }
  if (!is.numeric(n)) {
    stop("Argument 'n' must be numeric.")
  }
  ms <- pnd_reparm(lambda, xi, tau)
  x <- rpnd(n, lambda, ms$mu, ms$sigma)
  return(x)
}
