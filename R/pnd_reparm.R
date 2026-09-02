#' Reparameterization for the Power Normal Distribution
#'
#' Computes the conventional location and scale parameters from parameters
#' based on the median and relative interquartile range, as proposed by
#' Maruo et al. (2017).
#'
#' @param lambda Numeric scalar: shape parameter.
#' @param xi Numeric scalar: median on the original measurement scale.
#' @param tau Numeric scalar: relative interquartile range, defined as the
#'   interquartile range divided by the median.
#'
#' @details
#' In the power normal distribution, the location and scale parameters
#' \eqn{\mu} and \eqn{\sigma} correspond to the mean and standard deviation
#' on the transformed scale. Because both depend on the shape parameter
#' \eqn{\lambda}, they are not directly interpretable on the original scale.
#'
#' Maruo et al. (2017) proposed a reparameterization based on
#' \itemize{
#'   \item the median \eqn{\xi = Q(0.5)}, and
#'   \item the relative interquartile range
#'     \eqn{\tau = \{Q(0.75)-Q(0.25)\}/Q(0.5)},
#' }
#' where \eqn{Q(p)} denotes the \eqn{100p}th percentile of the normalized
#' power normal distribution.
#'
#' For \eqn{\lambda=0}, the distribution is log-normal and the corresponding
#' \eqn{\mu} and \eqn{\sigma} are obtained in closed form. For
#' \eqn{\lambda\ne0}, the function solves a one-dimensional root-finding
#' problem after introducing the numerically stable parameter
#' \deqn{c = \lambda K = (1+\lambda\mu)/\sigma.}
#' Unlike \eqn{K}, \eqn{c} remains finite as \eqn{\lambda} approaches zero.
#' Conditional normal probabilities are evaluated on the log-probability
#' scale, and logarithmic quantile ratios are calculated using
#' \code{log1p()} to reduce cancellation near \eqn{\lambda=0}.
#'
#' Not every combination of \code{lambda} and \code{tau} corresponds to a
#' finite power normal distribution. If the requested relative interquartile
#' range is outside the admissible range, or if a numerically stable solution
#' cannot be obtained, the function returns \code{NA} values with a warning.
#'
#' @return A list with components \code{mu} and \code{sigma},
#'   corresponding to the conventional location and scale parameters.
#'
#' @references
#' Goto, M., & Inoue, T. (1980). Some properties of the power normal
#' distribution. \emph{Japanese Journal of Biometrics}, 1, 28--54.
#' \url{https://doi.org/10.5691/jjb.1.28}
#'
#' Maruo, K., Yamabe, T., & Yamaguchi, Y. (2017). Statistical simulation based
#' on right-skewed distributions. \emph{Computational Statistics}, 32(3),
#' 889--907. \url{https://doi.org/10.1007/s00180-016-0664-4}
#'
#' @seealso \code{\link[powerNormal]{dpnd_rp}},
#'   \code{\link[powerNormal]{ppnd_rp}},
#'   \code{\link[powerNormal]{qpnd_rp}},
#'   \code{\link[powerNormal]{rpnd_rp}}
#'
#' @examples
#' pnd_reparm(lambda = 0.5, xi = 100, tau = 1)
#' pnd_reparm(lambda = 0.001, xi = 100, tau = 1)
#' pnd_reparm(lambda = 0, xi = 100, tau = 1)
#'
#' @name pnd_reparm
#' @rdname pnd_reparm
#' @export
pnd_reparm <- function(lambda, xi, tau) {

  ## ---- diagnostics ----
  debug_on <- isTRUE(getOption("powerNormal.debug"))

  warn_result <- function(message, debug = NULL) {
    if (debug_on && !is.null(debug)) {
      message <- paste0(message, " [debug: ", debug, "]")
    }
    warning(message, call. = FALSE)
    data.frame(mu = NA_real_, sigma = NA_real_)
  }

  ## ---- input validation ----
  args <- list(lambda = lambda, xi = xi, tau = tau)
  scalar_numeric <- vapply(
    args,
    function(x) is.numeric(x) && length(x) == 1L,
    logical(1)
  )

  if (!all(scalar_numeric)) {
    stop(
      "Arguments 'lambda', 'xi', and 'tau' must be scalar numeric values.",
      call. = FALSE
    )
  }
  if (!all(is.finite(c(lambda, xi, tau)))) {
    stop(
      "Arguments 'lambda', 'xi', and 'tau' must be finite.",
      call. = FALSE
    )
  }
  if (xi <= 0) {
    stop("Argument 'xi' must be positive.", call. = FALSE)
  }
  if (tau <= 0) {
    stop("Argument 'tau' must be positive.", call. = FALSE)
  }

  log_xi <- log(xi)
  z_075 <- stats::qnorm(0.75)

  ## Values indistinguishable from zero at double precision use the
  ## continuous log-normal limit.
  if (abs(lambda) <= .Machine$double.eps) {
    mu <- log_xi
    sigma <- asinh(tau / 2) / z_075
    return(data.frame(mu = mu, sigma = sigma))
  }

  probs <- c(0.25, 0.50, 0.75)
  log_tau_target <- log(tau)

  ## ---- conditional standard-normal quantiles ----
  ## With c = lambda*K, A(K) = Phi(c/abs(lambda)).  Log probabilities
  ## avoid underflow under pronounced truncation.
  z_star <- function(c_value, p = probs) {
    log_A <- stats::pnorm(c_value / abs(lambda), log.p = TRUE)

    if (lambda > 0) {
      ## Upper-tail probability is A(K)*(1-p).
      stats::qnorm(
        log_A + log1p(-p),
        lower.tail = FALSE,
        log.p = TRUE
      )
    } else {
      ## Lower-tail probability is A(K)*p.
      stats::qnorm(
        log_A + log(p),
        lower.tail = TRUE,
        log.p = TRUE
      )
    }
  }

  ## ---- log relative IQR for a given c ----
  log_tau_from_c <- function(c_value) {
    z <- z_star(c_value)
    denominator <- c_value + lambda * z[2L]

    if (!is.finite(denominator) || denominator <= 0) {
      return(NA_real_)
    }

    delta <- lambda * (z - z[2L]) / denominator
    if (any(!is.finite(delta)) || any(1 + delta <= 0)) {
      return(NA_real_)
    }

    ## log{Q(p)/xi}; log1p() is stable when lambda is near zero.
    log_q_ratio <- log1p(delta) / lambda
    log_ratio_difference <- log_q_ratio[1L] - log_q_ratio[3L]

    if (is.na(log_ratio_difference) || log_ratio_difference > 0) {
      return(NA_real_)
    }

    ## log(exp(a)-exp(b)), with a=log Q(.75)/xi and b=log Q(.25)/xi.
    log_q_ratio[3L] + log(-expm1(log_ratio_difference))
  }

  ## ---- admissible limiting relative IQR ----
  ## As c tends to -Inf, the relative IQR approaches a finite limit for
  ## fixed nonzero lambda.  Work on the log scale to avoid overflow.
  log_difference_exp <- function(log_x, log_y) {
    if (log_y > log_x) {
      stop("Internal error in log-difference calculation.", call. = FALSE)
    }
    if (is.infinite(log_x) && log_x > 0) {
      return(Inf)
    }
    log_x + log(-expm1(log_y - log_x))
  }

  r_limit <- log(4 / 3) / log(2)
  limit_terms <- c(log(2) / lambda, log(r_limit) / lambda)
  log_tau_limit <- log_difference_exp(
    max(limit_terms),
    min(limit_terms)
  )

  if (is.finite(log_tau_limit) && log_tau_target >= log_tau_limit) {
    tau_limit <- exp(log_tau_limit)
    return(warn_result(
      paste0(
        "No finite parameterization exists for the requested 'lambda' and ",
        "'tau'. For lambda = ", format(lambda, digits = 6),
        ", the limiting maximum of 'tau' is ",
        format(tau_limit, digits = 6), "."
      )
    ))
  }

  objective <- function(c_value) {
    log_tau_from_c(c_value) - log_tau_target
  }

  ## ---- construct a root bracket adaptively ----
  f_zero <- objective(0)
  if (is.na(f_zero)) {
    return(warn_result(
      "Failed to evaluate the reparameterization at the initial value.",
      debug = "objective(0) is NA"
    ))
  }

  root_tolerance <- 1e-12
  max_expansions <- 60L
  c_hat <- NA_real_
  bracket <- c(NA_real_, NA_real_)
  bracket_values <- c(NA_real_, NA_real_)

  if (is.finite(f_zero) && abs(f_zero) <= root_tolerance) {
    c_hat <- 0
  } else if (f_zero > 0) {
    ## The target lies to the right because tau(c) decreases with c.
    lower <- 0
    f_lower <- f_zero
    upper <- 1

    for (i in seq_len(max_expansions)) {
      f_upper <- objective(upper)
      if (is.finite(f_upper) && f_upper <= 0) {
        bracket <- c(lower, upper)
        bracket_values <- c(f_lower, f_upper)
        break
      }
      upper <- upper * 2
      if (!is.finite(upper)) break
    }
  } else {
    ## The target lies to the left because tau(c) decreases with c.
    upper <- 0
    f_upper <- f_zero
    lower <- -1

    for (i in seq_len(max_expansions)) {
      f_lower <- objective(lower)
      if (is.finite(f_lower) && f_lower >= 0) {
        bracket <- c(lower, upper)
        bracket_values <- c(f_lower, f_upper)
        break
      }
      lower <- lower * 2
      if (!is.finite(lower)) break
    }
  }

  if (is.na(c_hat)) {
    if (anyNA(bracket)) {
      return(warn_result(
        paste0(
          "Failed to bracket a numerically stable solution for the requested ",
          "'lambda', 'xi', and 'tau'."
        ),
        debug = paste0(
          "lambda=", signif(lambda, 6),
          ", xi=", signif(xi, 6),
          ", tau=", signif(tau, 6),
          ", objective(0)=", signif(f_zero, 6)
        )
      ))
    }

    root_fit <- tryCatch(
      stats::uniroot(
        objective,
        interval = bracket,
        f.lower = bracket_values[1L],
        f.upper = bracket_values[2L],
        tol = root_tolerance
      ),
      error = function(e) e
    )

    if (inherits(root_fit, "error")) {
      return(warn_result(
        "Root finding failed for the requested reparameterization.",
        debug = conditionMessage(root_fit)
      ))
    }
    c_hat <- root_fit$root
  }

  ## ---- recover mu and sigma ----
  z_hat <- z_star(c_hat)
  denominator <- c_hat + lambda * z_hat[2L]

  if (!is.finite(denominator) || denominator <= 0) {
    return(warn_result(
      "The root produced an invalid scale parameter.",
      debug = paste0("c=", signif(c_hat, 8),
                     ", denominator=", signif(denominator, 8))
    ))
  }

  log_sigma <- lambda * log_xi - log(denominator)
  sigma <- exp(log_sigma)
  mu <- expm1(lambda * log_xi) / lambda - sigma * z_hat[2L]

  if (!is.finite(mu) || !is.finite(sigma) || sigma <= 0) {
    return(warn_result(
      "The requested reparameterization produced non-finite parameters.",
      debug = paste0("c=", signif(c_hat, 8),
                     ", mu=", signif(mu, 8),
                     ", sigma=", signif(sigma, 8))
    ))
  }

  ## ---- numerical verification ----
  log_median <- log1p(lambda * (mu + sigma * z_hat[2L])) / lambda
  median_relative_error <- abs(expm1(log_median - log_xi))
  tau_relative_error <- abs(expm1(
    log_tau_from_c(c_hat) - log_tau_target
  ))

  verification_tolerance <- 1e-8
  if (!is.finite(median_relative_error) ||
      !is.finite(tau_relative_error) ||
      median_relative_error > verification_tolerance ||
      tau_relative_error > verification_tolerance) {
    return(warn_result(
      "The reparameterization could not be verified to the required accuracy.",
      debug = paste0(
        "c=", signif(c_hat, 8),
        ", median relative error=", signif(median_relative_error, 6),
        ", tau relative error=", signif(tau_relative_error, 6)
      )
    ))
  }

  return(list(mu = mu, sigma = sigma))
}
