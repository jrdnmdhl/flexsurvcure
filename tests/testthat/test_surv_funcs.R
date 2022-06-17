test_that("Mean survival works", {

  # MIXTURE MODELS
  # Mean should be infinite for models with positive cure fraction
  mix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="llogis")
  mix_some_cured_res <- summary(mix_some_cured, type = "mean", tidy=T)
  expect_equal(as.numeric(mix_some_cured_res), c(Inf, Inf, Inf))

  # For cure fraction = 0, mean should equal that of base distribution
  mix_none_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="identity", dist="exp")
  mix_none_cured_res <- summary(mix_none_cured, type = "mean", tidy=T)
  mix_none_cured_base_res <- mean_exp(mix_none_cured$res[2,1])
  expect_equal(mix_none_cured_res[1,1], mix_none_cured_base_res)

  # NON-MIXTURE MODELS
  # Mean should be infinite for models with positive cure fraction
  nmix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="llogis", mixture = F)
  nmix_some_cured_res <- summary(nmix_some_cured, type = "mean", tidy=T)
  expect_equal(as.numeric(nmix_some_cured_res), c(Inf, Inf, Inf))


  # Test case where theta is zero
  expect_equal(
    mean_nmixsurv(pgenf, 0, mu = 1.2, sigma = 0.8, Q = 0.2, P = 0.3),
    0
  )
  expect_equal(
    mean_mixsurv(pgenf, 0, mu = 1.2, sigma = 0.8, Q = 0.2, P = 0.3),
    mean_genf(1.2, 0.8, 0.2, 0.3)
  )

  # Test with vector arguments
  expect_equal(
    mean_nmixsurv(pweibull, c(0.1, 0.1, 0.1), shape=c(1.2, 1.3, 1.5), scale=c(20, 21, 50)),
    c(Inf, Inf, Inf)
  )
  expect_equal(
    mean_nmixsurv(pweibull, c(0.1), shape=c(1.2, 1.3, 1.5), scale=c(20, 21, 50)),
    c(Inf, Inf, Inf)
  )
  expect_equal(
    mean_nmixsurv(pweibull, c(0), shape=c(1.2, 1.3, 1.5), scale=c(20, 21, 50)),
    c(0, 0, 0)
  )
  expect_equal(
    mean_nmixsurv(pweibull, c(0, 1, 1), shape=c(1.2), scale=c(20)),
    c(0, Inf, Inf)
  )
  expect_error(
    mean_nmixsurv(pweibull, c(0, 1, 1), shape=c(1.2, 1.3), scale=c(20, 21)),
    'Parameter values provided were of incompatible length'
  )
  expect_equal(
    mean_mixsurv(pweibull, c(0.1, 0.1, 0.1), shape=c(1.2, 1.3, 1.5), scale=c(20, 21, 50)),
    c(Inf, Inf, Inf)
  )
  expect_equal(
    mean_mixsurv(pweibull, c(0.1), shape=c(1.2, 1.3, 1.5), scale=c(20, 21, 50)),
    c(Inf, Inf, Inf)
  )
  expect_equal(
    mean_mixsurv(pweibull, c(0, 0, 0), shape=c(1.2, 1.3, 1.5), scale=c(20, 21, 50)),
    c(18.81311716, 19.39511074, 45.13726464)
  )
  expect_equal(
    mean_mixsurv(pweibull, c(0), shape=c(1.2, 1.3, 1.5), scale=c(20, 21, 50)),
    c(18.81311716, 19.39511074, 45.13726464)
  )
  expect_equal(
    mean_mixsurv(pweibull, c(0, 1, 1), shape=c(1.2), scale=c(20)),
    c(18.81311716, Inf, Inf)
  )
  expect_error(
    mean_mixsurv(pweibull, c(0, 1, 1), shape=c(1.2, 1.3), scale=c(20, 21)),
    'Parameter values provided were of incompatible length'
  )

})

test_that("RMST Works", {
  # MIXTURE MODELS
  # RMST should be equal to duration * theta + (1-theta) * uncured_rmst
  t_rmst <- 10000
  mix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="gompertz")
  mix_some_cured_res <- summary(mix_some_cured, t=t_rmst, type = "rmst", tidy=T)
  mix_some_cured_res_u <- rmst_gompertz(
    t = t_rmst,
    shape=mix_some_cured$res[2,1],
    rate=mix_some_cured$res[3,1]
  )
  expect_equal(
    mix_some_cured_res$est,
    (1 - mix_some_cured$res[1,1]) * mix_some_cured_res_u + mix_some_cured$res[1,1] * t_rmst
  )

  nmix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="gompertz", mix = F)
  nmix_some_cured_res <- summary(nmix_some_cured, t=t_rmst, type = "rmst", tidy=T)
  expect_equal(
    nmix_some_cured_res$est,
    integrate(function(x) pnmixsurv(
      pgompertz,
      x,
      nmix_some_cured$res[1,1],
      shape = nmix_some_cured$res[2,1],
      rate = nmix_some_cured$res[3,1],
      lower.tail = F), 0, t_rmst
    )$value
  )
})

test_that("Survival projections", {

  # MIXTURE MODELS
  # Survival should be equal to cure fraction for large values of t
  mix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="loglog", dist="lnorm")
  mix_some_cured_res <- summary(mix_some_cured, t=1e99, type = "survival", tidy=T)
  expect_equal(as.numeric(mix_some_cured_res)[2], as.numeric(mix_some_cured$res[1]))

  # NON-MIXTURE MODELS
  # Survival should be equal to cure fraction for large values of t
  nmix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="loglog", dist="lnorm", mixture=F)
  nmix_some_cured_res <- summary(nmix_some_cured, t=1e99, type = "survival", tidy=T)
  expect_equal(as.numeric(nmix_some_cured_res)[2], as.numeric(nmix_some_cured$res[1]))
})

test_that("Cumulative hazard projections", {

  # MIXTURE MODELS
  # Cumulative hazard should equal -log(cure fraction) for large values of t
  mix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="identity", dist="weibull")
  mix_some_cured_res <- summary(mix_some_cured, t=1e99, type = "cumhaz", tidy=T)
  expect_equal(as.numeric(mix_some_cured_res)[2], -log(as.numeric(mix_some_cured$res[1])))

  # NON-MIXTURE MODELS
  # Cumulative hazard should equal -log(cure fraction) for large values of t
  nmix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="identity", dist="weibull", mixture=F)
  nmix_some_cured_res <- summary(nmix_some_cured, t=1e99, type = "cumhaz", tidy=T)
  expect_equal(as.numeric(nmix_some_cured_res)[2], -log(as.numeric(nmix_some_cured$res[1])))

  # MIXTURE MODELS
  # Cumulative hazard should flatten
  mix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="exp")
  mix_some_cured_res <- summary(mix_some_cured, t=c(9999999, 1e99), type = "cumhaz", tidy=T)
  expect_equal(mix_some_cured_res$est[1], mix_some_cured_res$est[2])

  # NON-MIXTURE MODELS
  # Cumulative hazard should flatten
  nmix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="exp", mixture=F)
  nmix_some_cured_res <- summary(nmix_some_cured, t=c(9999999, 1e99), type = "cumhaz", tidy=T)
  expect_equal(nmix_some_cured_res$est[1], nmix_some_cured_res$est[2])
})

test_that("Hazard rate projections", {

  # MIXTURE MODELS
  # Hazard at t = Inf should be zero
  mix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="exp")
  mix_some_cured_res <- summary(mix_some_cured, t=Inf, type = "hazard", tidy=T)
  expect_equal(as.numeric(mix_some_cured_res)[2], 0)

  # NON-MIXTURE MODELS
  # Hazard at t = Inf should be zero
  nmix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="exp", mixture=F)
  nmix_some_cured_res <- summary(nmix_some_cured, t=Inf, type = "hazard", tidy=T)
  expect_equal(as.numeric(nmix_some_cured_res)[2], 0)
})

test_that("Random sampling", {

  # MIXTURE MODELS
  expect_equal(mean(rmixsurv(qexp, n = 1000000, theta = 0.0, rate = 1/50)), 50, tolerance = 1e-2)
  expect_equal(median(rmixsurv(qexp, n = 10000000, theta = 0.20, rate = 1/50)), qexp(0.625, rate = 1/50), tolerance = 1e-1)

  # NON-MIXTURE MODELS
  expect_equal(mean(rnmixsurv(qexp, n = 1000000, theta = 0.0, rate = 1/50)), 50, tolerance = 1e-2)

})

test_that("P function works with infinite input", {

  # MIXTURE MODELS
  # Hazard at t = Inf should be zero
  mix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="exp")
  mix_some_cured_res <- summary(mix_some_cured, t=Inf, type = "survival", tidy=T)
  expect_equal(as.numeric(mix_some_cured_res)[2], 0)

  # NON-MIXTURE MODELS
  # Hazard at t = Inf should be zero
  nmix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="exp", mixture=F)
  nmix_some_cured_res <- summary(nmix_some_cured, t=Inf, type = "survival", tidy=T)
  expect_equal(as.numeric(nmix_some_cured_res)[2], 0)
})

test_that("Quantile functions", {

  # MIXTURE MODELS
  expect_equal(
    qmixsurv(qexp, c(0.25, 0.5, 0.75), theta = 0.2, rate = 1/50),
    qgeneric(function(...) pmixsurv(pexp, ...), c(0.25, 0.5, 0.75), theta = 0.2, rate = 1/50),
    tolerance = 1e-4
  )

  expect_equal(
    qmixsurv(qweibull, c(0.25, 0.5, 0.75), theta = 0.15, shape = 1.2, scale = 50),
    qgeneric(function(...) pmixsurv(pweibull, ...), c(0.25, 0.5, 0.75), theta = 0.15, shape = 1.2, scale = 50),
    tolerance = 1e-4
  )

  # NON-MIXTURE MODELS
  expect_equal(
    qnmixsurv(qexp, c(0.25, 0.5, 0.75), theta = 0.2, rate = 1/50),
    qgeneric(function(...) pnmixsurv(pexp, ...), c(0.25, 0.5, 0.75), theta = 0.2, rate = 1/50),
    tolerance = 1e-4
  )

  expect_equal(
    qnmixsurv(qweibull, c(0.25, 0.5, 0.75), theta = 0.15, shape = 1.2, scale = 50),
    qgeneric(function(...) pnmixsurv(pweibull, ...), c(0.25, 0.5, 0.75), theta = 0.15, shape = 1.2, scale = 50),
    tolerance = 1e-4
  )

})

test_that("Probit link works", {
  probit_model <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="probit", dist="llogis")
  expect_equal(probit_model$res.t[1,1], qnorm(probit_model$res[1,1]))
  expect_equal(probit_model$res[1,1], pnorm(probit_model$res.t[1,1]))
})


