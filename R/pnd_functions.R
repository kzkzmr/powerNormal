#' The Power Normal Distribution
#'
#' Probability density function (PDF), cumulative distribution function (CDF),
#' quantile function, and random number generation for the power normal
#' distribution (PND).
#'
#' @param x Vector of quantiles.
#' @param p Vector of probabilities.
#' @param n Number of observations.
#' @param lambda Numeric scalar; the shape parameter.
#' @param mu Numeric scalar; the location parameter.
#' @param sigma Numeric scalar; the scale parameter.
#' @param log Logical; if \code{TRUE}, returns the log-density.
#'
#' @details
#' The power normal distribution (PND; Goto and Inoue, 1980) describes the
#' distribution of a positive random variable \eqn{X} such that the Box–Cox
#' transformation \eqn{Z(\lambda)} follows a (possibly truncated) normal
#' distribution (Box & Cox, 1964):
#'
#' \deqn{
#' Z(\lambda) =
#' \begin{cases}
#' (X^\lambda - 1)/\lambda, & \lambda \neq 0, \\
#' \log(X), & \lambda = 0.
#' \end{cases}
#' }
#'
#' The PDF of the PND is:
#' \deqn{
#' f_{\mathrm{PN}}(x; \lambda, \mu, \sigma) =
#' \frac{x^{\lambda - 1}}{A(K)} f_{\mathrm{N}}(z(\lambda); \mu, \sigma),
#' }
#' where \eqn{f_{\mathrm{N}}} is the normal PDF, \eqn{z(\lambda)} is the
#' Box–Cox transformed value, and \eqn{A(K)} is the normalization constant:
#'
#' \deqn{
#' A(K) =
#' \begin{cases}
#' \Phi(\mathrm{sgn}(\lambda)K), & \lambda \neq 0, \\
#' 1, & \lambda = 0,
#' \end{cases}
#' }
#' with \eqn{K = (1 + \lambda \mu)/(\lambda \sigma)} and \eqn{\Phi} the
#' standard normal CDF.
#'
#' The CDF is:
#' \deqn{
#' F_{\mathrm{PN}}(x; \lambda, \mu, \sigma) =
#' \begin{cases}
#' \frac{1}{A(K)}\left\{ \Phi\left( \frac{z(\lambda) - \mu}{\sigma} \right) -
#' \Phi(-K) \right\}, & \lambda > 0, \\
#' \frac{1}{A(K)} \Phi\left( \frac{z(\lambda) - \mu}{\sigma} \right), & \lambda
#' \le 0.
#' \end{cases}
#' }
#'
#' The quantile function (i.e., the \eqn{100p} percentile \eqn{\xi_p}) is:
#' \deqn{
#' \xi_p(\lambda, \mu, \sigma) =
#' \begin{cases}
#' \left\{ \lambda (\mu + \sigma z_{p^*}) + 1 \right\}^{1/\lambda}, & \lambda
#' \ne 0, \\
#' \exp(\mu + \sigma z_p), & \lambda = 0,
#' \end{cases}
#' }
#' where \eqn{z_p} and \eqn{z_{p^*}} are quantiles of the standard normal
#' distribution, and
#'
#' \deqn{
#' p^* =
#' \begin{cases}
#' 1 - A(K)(1 - p), & \lambda > 0, \\
#' A(K)p, & \lambda < 0.
#' \end{cases}
#' }
#'
#' In practice, parameters are often estimated assuming \eqn{A(K) = 1}
#' (i.e., the truncation is negligible). This simplification is generally
#' acceptable unless extreme quantiles (e.g., 97.5%) are of interest.
#' Estimation under this assumption can be performed using \code{bct.v()}
#' from the \pkg{bcmixed} package.
#'
#' @return
#' \code{dpnd()} returns the PDF,
#' \code{ppnd()} returns the CDF,
#' \code{qpnd()} returns the quantile function, and
#' \code{rpnd()} generates random values from the distribution.
#'
#' @references
#' Box, G.E.P., & Cox, D.R. (1964). An analysis of transformations.
#' \emph{Journal of the Royal Statistical Society: Series B}, 26, 211–246.
#' \url{https://doi.org/10.1111/j.2517-6161.1964.tb00553.x}
#'
#' Goto, M., & Inoue, T. (1980). Some properties of the power normal
#' distribution. \emph{Japanese Journal of Biometrics}, 1, 28–54.
#' \url{https://doi.org/10.5691/jjb.1.28}
#'
#' Maruo, K., Shirahata, S., & Goto, M. (2011). Underlying assumptions of the
#' power-normal distribution. \emph{Behaviormetrika}, 38(1), 85–95.
#' \url{https://doi.org/10.2333/bhmk.38.85}
#'
#' Maruo, K., Ishii, R., Yamaguchi, Y., & Gosho, M. (2021). bcmixed: A package
#' for median inference on longitudinal data with the Box–Cox transformation.
#' \emph{The R Journal}, 13(2), 253–265.
#' \url{https://doi.org/10.32614/RJ-2021-083}
#'
#' @seealso
#' \code{\link[powerNormal]{dpnd}},
#' \code{\link[powerNormal]{ppnd}},
#' \code{\link[powerNormal]{qpnd}},
#' \code{\link[powerNormal]{rpnd}},
#' \code{\link[bcmixed]{bct.v}}
#'
#' @examples
#' dpnd(x = 5, lambda = 0.5, mu = 10, sigma = 2)
#' ppnd(x = 5, lambda = 0.5, mu = 10, sigma = 2)
#' qpnd(p = 0.6, lambda = -0.5, mu = 1, sigma = 0.1)
#' rpnd(n = 20, lambda = -0.5, mu = 1, sigma = 0.1)
#'
#' @importFrom bcmixed bct
#' @importFrom stats dnorm pnorm qnorm rnorm
#'
#' @name pnd_functions
#' @rdname pnd_functions
#' @export
dpnd <- function(x, lambda, mu, sigma, log= FALSE) {
  if (!is.numeric(lambda) | length(lambda) != 1 | !is.numeric(mu) |
      length(mu) != 1 | !is.numeric(sigma) | length(sigma) != 1){
    stop("Arguments 'lambda', 'mu', and 'sigma' must be scalar numeric values.")
  }
  if (!is.numeric(x) | sum(x < 0) != 0) {
    stop("Argument 'x' must be a positive numeric vector.")
  }
  z <- bct(x, lambda)
  if (lambda == 0) {
    AK <- 1
  } else {
    K <- (1 + lambda * mu) / (lambda * sigma)
    AK <- pnorm(sign(lambda) * K)
  }
  if (log) {
    fx <- (lambda - 1) * log(x) - log(AK) + dnorm(z, mu, sigma, log = TRUE)
  } else {
    fx <- x ^ (lambda - 1) / AK * dnorm(z, mu, sigma)
  }
  return(fx)
}

#' @rdname pnd_functions
#' @export
ppnd <- function(x, lambda, mu, sigma) {
  if (!is.numeric(lambda) | length(lambda) != 1 | !is.numeric(mu) |
      length(mu) != 1 | !is.numeric(sigma) | length(sigma) != 1){
    stop("Arguments 'lambda', 'mu', and 'sigma' must be scalar numeric values.")
  }
  if (!is.numeric(x) | sum(x < 0) != 0) {
    stop("Argument 'x' must be a positive numeric vector.")
  }
  z <- bct(x, lambda)
  if (lambda == 0) {
    AK <- 1
  } else {
    K <- (1 + lambda * mu) / (lambda * sigma)
    AK <- pnorm(sign(lambda) * K)
  }
  if (lambda > 0) {
    Fx <- 1 / AK * (pnorm(z, mu, sigma) - pnorm(-K))
  } else {
    Fx <- 1 / AK * pnorm(z, mu, sigma)
  }
  return(Fx)
}

#' @rdname pnd_functions
#' @export
qpnd <- function(p, lambda, mu, sigma){
  if (!is.numeric(lambda) | length(lambda) != 1 | !is.numeric(mu) |
      !is.numeric(sigma)){
    stop("Arguments 'lambda', 'mu', and 'sigma' must be scalar numeric values.")
  }
  if (!is.numeric(p)) {
    stop("Argument 'p' must be a numeric vector")
  } else if (any(p < 0 | p > 1)) {
    stop("All elements of 'p' must be between 0 and 1.")
  }
  zps <- function(lambda, K, p){
    AK <- pnorm(sign(lambda) * K)
    if (lambda < 0) ps <- AK * p
    if (lambda > 0) ps <- 1 - AK * (1 - p)
    zps <- qnorm(ps)
    return(zps)
  }
  if(lambda != 0){
    K <- (1 + lambda * mu) / (lambda * sigma)
    zps <- zps(lambda, K, p)
    q <- (lambda * (mu + sigma * zps) + 1) ^ (1 / lambda)
  }
  else {
    q <- exp(mu + sigma * qnorm(p))
  }
  return(q)
}


#' @rdname pnd_functions
#' @export
rpnd <- function(n, lambda, mu, sigma){
  if (!is.numeric(lambda) | length(lambda) != 1 | !is.numeric(mu) |
      length(mu) != 1 | !is.numeric(sigma) | length(sigma) != 1){
    stop("Arguments 'lambda', 'mu', and 'sigma' must be scalar numeric values.")
  }
  if (!is.numeric(n)) {
    stop("Argument 'n' must be numeric.")
  }
  n2 <- n
  if (length(n) > 1) {
    n2 <- length(n)
  }
  if (lambda == 0) {
    AK <- 1
  } else {
    K <- (1 + lambda * mu) / (lambda * sigma)
    AK <- pnorm(sign(lambda) * K)
  }

  n3 <- round(n2 / AK * 2)
  z0 <- rnorm(n3, mu, sigma)
  cut0 <- numeric(n3) == 0
  if (lambda < 0){
    cut0 <- z0 < -1 / lambda
  } else {
    cut0 <- z0 > -1 / lambda
  }
  z <- z0[cut0][1:n2]
  if (lambda == 0) {
    x <- exp(z)
  } else {
    x <- (lambda * z + 1) ^ (1 / lambda)
  }
  return(x)
}
