test_that("lambda = 0 agrees with the log-normal closed form", {
  xi <- 100
  tau <- 1
  fit <- pnd_reparm(lambda = 0, xi = xi, tau = tau)

  expect_equal(fit$mu, log(xi), tolerance = 1e-12)
  expect_equal(
    fit$sigma,
    asinh(tau / 2) / qnorm(0.75),
    tolerance = 1e-12
  )
})


test_that("the requested median and relative IQR are recovered", {
  settings <- expand.grid(
    lambda = c(-2, -1, -0.5, -0.01, -0.001,
               0, 0.001, 0.01, 0.5, 1, 2),
    xi = c(1, 100),
    tau = c(0.1, 0.5),
    KEEP.OUT.ATTRS = FALSE
  )

  for (i in seq_len(nrow(settings))) {
    s <- settings[i, ]
    fit <- pnd_reparm(s$lambda, s$xi, s$tau)
    q <- qpnd(
      c(0.25, 0.50, 0.75),
      s$lambda,
      fit$mu,
      fit$sigma
    )

    expect_equal(
      unname(q[2]),
      s$xi,
      tolerance = 1e-7,
      info = paste("lambda =", s$lambda,
                   "xi =", s$xi,
                   "tau =", s$tau)
    )
    expect_equal(
      unname((q[3] - q[1]) / q[2]),
      s$tau,
      tolerance = 1e-7,
      info = paste("lambda =", s$lambda,
                   "xi =", s$xi,
                   "tau =", s$tau)
    )
  }
})


test_that("the reparameterization is continuous around lambda = 0", {
  xi <- 100
  tau <- 1

  fit_negative <- pnd_reparm(-1e-6, xi, tau)
  fit_zero <- pnd_reparm(0, xi, tau)
  fit_positive <- pnd_reparm(1e-6, xi, tau)

  expect_equal(fit_negative$mu, fit_zero$mu, tolerance = 2e-5)
  expect_equal(fit_positive$mu, fit_zero$mu, tolerance = 2e-5)
  expect_equal(fit_negative$sigma, fit_zero$sigma, tolerance = 2e-5)
  expect_equal(fit_positive$sigma, fit_zero$sigma, tolerance = 2e-5)
})


test_that("infeasible relative IQR values return NA with a warning", {
  expect_warning(
    fit <- pnd_reparm(lambda = 2, xi = 100, tau = 1),
    "No finite parameterization exists"
  )
  expect_true(is.na(fit$mu))
  expect_true(is.na(fit$sigma))
})


test_that("invalid inputs are rejected", {
  expect_error(pnd_reparm(NA_real_, 100, 1), "must be finite")
  expect_error(pnd_reparm(0.5, 0, 1), "'xi' must be positive")
  expect_error(pnd_reparm(0.5, 100, 0), "'tau' must be positive")
  expect_error(pnd_reparm(c(0.2, 0.3), 100, 1), "scalar numeric")
})
