#' Reparametrization for the Power Normal Distribution
#'
#' Computes the original location and scale parameters from the reparametrized
#' location (median) and scale (interquantile range / median) parameters,
#' as proposed by Maruo et al. (2017).
#'
#' @param lambda Numeric scalar: shape parameter.
#' @param xi Numeric scalar: reparametrized location parameter (i.e., median).
#' @param tau Numeric scalar: reparametrized scale parameter (i.e.,
#'   interquantile range divided by median).
#'
#' @details
#' In the power normal distribution, the location and scale parameters
#' \eqn{\mu} and \eqn{\sigma} correspond to the mean and standard deviation
#' parameters on the transformed scale. Since both depend strongly on the
#' shape parameter \eqn{\lambda}, they are not easy to interpret.
#'
#' This makes it difficult to directly specify meaningful parameter values in
#' simulations. To address this, Maruo et al. (2017) proposed a
#' reparametrization based on:
#'
#' \itemize{
#'   \item the median \eqn{\xi_{0.5}}, and
#'   \item the scale-related quantity
#'         \eqn{\tau = (\xi_{0.75} - \xi_{0.25}) / \xi_{0.5}},
#' }
#'
#' where \eqn{\xi_p} denotes the \eqn{100p} percentile of the distribution.
#'
#' The \code{pnd_reparm} function numerically calculates the original
#' parameters \eqn{\mu} and \eqn{\sigma} corresponding to the given
#' reparametrized values \eqn{\xi_{0.5}} and \eqn{\tau}.
#' See Maruo et al. (2017) for the algorithmic details.
#'
#' @return A data frame with two elements:
#'   \code{mu} and \code{sigma}, corresponding to the original parameters.
#'
#' @references
#' Goto, M., & Inoue, T. (1980). Some properties of the power normal
#' distribution. \emph{Japanese Journal of Biometrics}, 1, 28–54.
#' \url{https://doi.org/10.5691/jjb.1.28}
#'
#' Maruo, K., Yamabe, T., & Yamaguchi, Y. (2017). Statistical simulation based on
#' right-skewed distributions. \emph{Computational Statistics}, 32(3), 889–907.
#' \url{https://doi.org/10.1007/s00180-016-0664-4}
#'
#' @seealso \code{\link[powerNormal]{dpnd_rp}},
#'   \code{\link[powerNormal]{ppnd_rp}},
#'   \code{\link[powerNormal]{qpnd_rp}},
#'   \code{\link[powerNormal]{rpnd_rp}}
#'
#' @examples
#' pnd_reparm(lambda = 0.5, xi = 100, tau = 1)
#'
#' @name pnd_reparm
#' @rdname pnd_reparm
#' @export
pnd_reparm <- function(lambda, xi, tau){
  if (!is.numeric(lambda) | length(lambda) != 1 | !is.numeric(xi) |
      length(xi) != 1 | !is.numeric(tau) | length(tau) != 1){
    stop("Arguments 'lambda', 'xi', and 'tau' must be scalar numeric values.")
  }
  zps <- function(lambda, K, p){
    AK <- pnorm(sign(lambda) * K)
    if (lambda < 0) ps <- AK * p
    if (lambda > 0) ps <- 1 - AK * (1 - p)
    zps <- qnorm(ps)
    return(zps)
  }
  xik_ms <- function(lambda, xi, K){
    z5 <- zps(lambda, K, 0.5)
    mu <- (K * (xi ^ lambda - 1) - z5)/(lambda * (K + z5))
    sigma <- (1 + lambda * mu) / (lambda * K)
    return(list(mu, sigma))
  }

  if (abs(lambda) > 0.01){
    K1 <- sign(lambda) * seq(0.1, 100, 0.1)
    kp <- length(K1)
    ms1 <- xik_ms(lambda, xi, K1)
    mu1 <- ms1[[1]]
    sigma1 <- ms1[[2]]
    d1 <- qpnd(0.75, lambda, mu1, sigma1) -
      qpnd(0.25, lambda , mu1, sigma1) - xi * tau
    ld1 <- which.min(abs(d1))
    if (sign(d1[ld1 - 1]) == sign(d1[ld1])) {
      K2 <- seq(K1[ld1], K1[ld1 + 1], sign(lambda) * 1e-4)
    } else {
      K2 <- seq(K1[ld1 - 1], K1[ld1], sign(lambda) * 1e-4)
    }
    kp <- length(K2)
    ms2 <- xik_ms(lambda, xi, K2)
    mu2 <- ms2[[1]]
    sigma2 <- ms2[[2]]
    d2 <- qpnd(0.75, lambda, mu2, sigma2)-
      qpnd(0.25, lambda, mu2, sigma2) - xi * tau
    ld2 <- which.min(abs(d2))
    mu <- mu2[ld2]
    sigma <- sigma2[ld2]
  }
  else{
    mu <- log(xi)
    sigma <- log((tau + sqrt(tau ^ 2 + 4)) / 2) / qnorm(0.75)
  }
  return(data.frame(mu = mu, sigma = sigma))
}
