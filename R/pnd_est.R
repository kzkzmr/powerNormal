#' Parameter Estimation for the Power Normal Distribution
#'
#' Estimates the parameters \code{lambda}, \code{mu}, and \code{sigma}
#' of the power normal distribution based on a given numeric vector.
#'
#' @param x A positive numeric vector of observed values.
#' @param tn Logical. If \code{TRUE}, estimation is based on a
#' truncated normal likelihood on the transformed scale (see Details).
#' Default is \code{FALSE}.
#' @param lmdint A numeric vector of length 2 specifying the interval
#' within which to search for the optimal shape parameter \code{lambda}.
#' Default is \code{c(-3, 3)}.
#'
#' @details
#' This function estimates the parameters of the power normal
#' distribution, which is defined as the inverse Box–Cox transformation
#' of a normal distribution. The estimated parameters are:
#' \code{lambda} (shape), \code{mu} (location), and \code{sigma} (scale).
#'
#' When \code{tn = FALSE} (default), the transformed values are assumed
#' to follow a normal distribution, and \code{mu} and \code{sigma}
#' are estimated by the sample mean and standard deviation on the transformed
#' scale.
#'
#' When \code{tn = TRUE}, the transformed values are assumed to follow
#' a truncated normal distribution. In this case, maximum likelihood
#' estimation is performed using numerical optimization to obtain
#' \code{mu} and \code{sigma}, and a correction is applied via a
#' normalizing constant, \eqn{A(K)} (Maruo et al., 2011).
#'
#' In practice, parameters are often estimated assuming \eqn{A(K) = 1}
#' (i.e., the truncation is negligible). This simplification is generally
#' acceptable unless extreme quantiles (e.g., 97.5%) are of interest.
#' For example, the BCMMRM method (Maruo et al., 2017; 2021), which focuses
#' on the median, makes inferences based on the assumption that \eqn{A(K) = 1}.
#'
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{lambda}}{Estimated shape parameter.}
#'   \item{\code{mu}}{Estimated location parameter.}
#'   \item{\code{sigma}}{Estimated scale parameter.}
#' }
#'
#' @references
#'
#' Maruo, K., Shirahata, S., & Goto, M. Underlying assumptions of the
#' power-normal distribution. \emph{Behaviormetrika}, 38(1), 85–95, 2011.
#' \url{https://doi.org/10.2333/bhmk.38.85}
#'
#' Maruo, K., Yamaguchi, Y., Noma, H., Gosho, M. Interpretable inference on
#' the mixed effect model with the Box-Cox transformation. \emph{Statistics in
#' Medicine}, 36, 2420-2434, 2017. \url{https://doi.org/10.1002/sim.7279}
#'
#' Maruo, K., Ishii, R., Yamaguchi, Y., & Gosho, M. (2021). bcmixed: A package
#' for median inference on longitudinal data with the Box–Cox transformation.
#' \emph{The R Journal}, 13(2), 253–265.
#' \url{https://doi.org/10.32614/RJ-2021-083}
#'
#' @seealso \code{\link[powerNormal]{pnd_functions}},
#'  \code{\link[powerNormal]{pnd_reparm}}
#'
#' @examples
#' # A(K)=0.826: case where truncation might not be negligible
#' x <- rpnd_rp(100, lambda = 0.4, xi = 100, tau = 2)
#' pnd_est(x)
#' pnd_est(x, tn = TRUE)
#'
#' @importFrom stats optimise sd
#' @importFrom bcmixed bct bct.v
#' @importFrom MASS ginv
#'
#' @export
pnd_est <- function(x, tn = FALSE, lmdint = c(-3, 3)) {
  pnd_est_tn <- function(x, lmdint){
    if (sum(x <= 0) > 0) {
      stop("all elements of x must be positive.")
    }
    tnd_est <- function(z, lambda){
      n <- length(z)
      ln_tnd <- function(z, lambda, mu, sigma){
        n <- length(z)
        K <- (1 + lambda * mu) / (lambda * sigma)
        AK <- pnorm(sign(lambda) * K)
        K <- (1 + lambda * mu) / (lambda * sigma)
        ln <- -n * log(AK) -
          n / 2 * log(2 * pi) - n * log(sigma) -
          sum((z - mu) ^ 2) / (2 * sigma ^ 2)
        return(ln)
      }
      dln_tnd <- function(z, lambda, mu, sigma) {
        n <- length(z)
        K <- (1 + lambda * mu) / (lambda * sigma)
        AK <- pnorm(sign(lambda) * K)
        eta <- dnorm(K) / AK
        dln <- numeric(2)
        dln[1] <- 1 / sigma ^ 2 * sum(z - mu) - sign(lambda) * n * eta / sigma
        dln[2] <- -n / sigma + 1 / sigma ^ 3 * sum((z - mu) ^ 2) +
          sign(lambda) * n * eta * K / sigma
        return(dln)
      }
      ddln_tnd <- function(z, lambda, mu, sigma) {
        n <- length(z)
        K <- (1 + lambda * mu) / (lambda * sigma)
        AK <- pnorm(sign(lambda) * K)
        eta <- dnorm(K) / AK
        ddln1 <- matrix(0, 2, 2)
        ddln2 <- ddln1
        ddln3 <- ddln1
        ddln1[1, 1] <- -n / sigma ^ 2
        ddln1[1, 2] <- -2 / sigma ^ 3 * sum(z - mu)
        ddln1[2, 1] <- ddln1[1, 2]
        ddln1[2, 2] <- n / sigma ^ 2 - 3 / sigma ^ 4 * sum((z - mu) ^ 2)
        ddln2[1, 1] <- K
        ddln2[1, 2] <- 1 - K ^ 2
        ddln2[2, 1] <- ddln2[1, 2]
        ddln2[2, 2] <- K * (K ^ 2 ^ 2)
        ddln3[1, 1] <- 1
        ddln3[1, 2] <- -K
        ddln3[2, 1] <- -K
        ddln3[2, 2] <- K ^ 2
        ddln <- ddln1 + sign(lambda) * n * eta / sigma ^ 2 * ddln2 +
          n * eta ^ 2 / sigma ^ 2 * ddln3
        return(ddln)
      }
      if (lambda != 0) {
        m1 <- 1 / n * sum(sign(lambda) * (z + 1 / lambda))
        m2 <- 1 / n * sum((sign(lambda) * (z + 1 / lambda)) ^ 2)
        m3 <- 1 / n * sum((sign(lambda) * (z + 1 / lambda)) ^ 3)
        mu0 <- sign(lambda) * (2 * m1 * m2 - m3) / (2 * m1 ^ 2 - m2) - 1 / lambda
        sigma02 <- (m1 * m3 - m2 ^ 2) / (2 * m1 ^2 - m2)
        if (sigma02 > 0) {
          sigma0 <- sqrt(sigma02)
        } else {
          sigma0 <- sd(z)
        }
        dlt <- 1e-5 * sigma0
        stpflg <- 0
        ln0 <- ln_tnd(z, lambda, mu0, sigma0)
        if (is.infinite(ln0)) {
          mu <- NA
          sigma <- NA
          } else {
          dln0 <- dln_tnd(z, lambda, mu0, sigma0)
          ddln0 <- ddln_tnd(z, lambda, mu0, sigma0)
          if (det(ddln0) > 0) {
            ddln0 <- diag(2)
          }
          count <- 1
          while (stpflg == 0) {
            dd <- -ginv(ddln0) %*% dln0
            zeta <- 1
            stpflg2 <- 0
            eta <- 1 / 3
            while (stpflg2 == 0) {
              mu1 <- mu0 + zeta * dd[1]
              sigma1 <- sigma0 + zeta * dd[2]
              if (sigma1 < 0) {
                zeta <- zeta * 0.5
              } else {
                ln1 <- ln_tnd(z, lambda, mu1, sigma1)
                if (ln1 >= ln0 + zeta * eta * sum(dln0 * dd)) {
                  stpflg2 <- 1
                } else {
                  zeta <- zeta * 0.5
                }
              }
              if (zeta < 1e-7) {
                stpflg2 <- 1
              }
            }
            if (sqrt((mu1 - mu0) ^ 2 + (sigma1 - sigma0) ^ 2) < dlt |
                count > 100 | zeta < 1e-7) {
              stpflg <- 1
            } else {
              mu0 <- mu1
              sigma0 <- sigma1
              ln0 <- ln_tnd(z, lambda, mu0, sigma0)
              dln0 <- dln_tnd(z, lambda, mu0, sigma0)
              ddln0 <- ddln_tnd(z, lambda, mu0, sigma0)
              if (det(ddln0) > 0) {
                ddln0 <- diag(2)
              }
              count <- count + 1
            }
          }
          mu <- mu1
          sigma <- sigma1
          if (zeta < 1e-10 | count > 100) {
            mu <- NA
            sigma <- NA
          }
        }

      } else {
        mu <- mean(z)
        sigma <- sd(z) * sqrt((n - 1) / n)
      }
      return(c(mu, sigma))
    }
    ln_pnd_tn <- function(x, lambda) {
      n <- length(x)
      z <- bct(x, lambda)
      ms <- tnd_est(z, lambda)
      if (is.na(ms)[1]){
        ln <- -1e+10
      } else {
        mu <- ms[1]
        sigma <- ms[2]
        if (lambda == 0) {
          AK <- 1
        } else {
          K <- (1 + lambda * mu) / (lambda * sigma)
          AK <- pnorm(sign(lambda) * K)
        }
        ln <- -n / 2 * log(2 * pi) - n * log(sigma) -
          sum((z - mu) ^ 2 / (2 * sigma ^ 2)) + (lambda - 1 ) * sum(log(x)) -
          n * log(AK)
      }
      return(ln)
    }
    lambda <- as.numeric(optimise(ln_pnd_tn, interval = lmdint, x = x,
                       maximum = TRUE)$maximum)
    z <- bct(x, lambda)
    ms <- tnd_est(z, lambda)
    return(list(lambda = lambda, mu = ms[1], sigma = ms[2]))
  }
  if (sum(x <= 0) > 0) {
    stop("all elements of x must be positive.")
  }
  if (tn) {
    est <- pnd_est_tn(x, lmdint)
  } else {
    lambda <- bct.v(x, lmdint)$lambda
    z <- bct(x, lambda)
    mu <- mean(z)
    n <- length(x)
    sigma <- sd(z) * (n - 1) / n
    est <- list(lambda = lambda, mu = mu, sigma = sigma)
  }
  return(est)
}
