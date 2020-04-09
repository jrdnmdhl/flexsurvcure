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
})

test_that("Survival projections", {

  # MIXTURE MODELS
  # Survival at t = Inf should equal cure fraction
  mix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="loglog", dist="lnorm")
  mix_some_cured_res <- summary(mix_some_cured, t=Inf, type = "survival", tidy=T)
  expect_equal(as.numeric(mix_some_cured_res)[2], as.numeric(mix_some_cured$res[1]))

  # NON-MIXTURE MODELS
  # Mean should be infinite for models with positive cure fraction
  nmix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="loglog", dist="lnorm", mixture=F)
  nmix_some_cured_res <- summary(nmix_some_cured, t=Inf, type = "survival", tidy=T)
  expect_equal(as.numeric(nmix_some_cured_res)[2], as.numeric(nmix_some_cured$res[1]))
})

test_that("Cumulative hazard projections", {

  # MIXTURE MODELS
  # Cumulaitve hazard at t = Inf should equal -log(cure fraction0
  mix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="identity", dist="weibull")
  mix_some_cured_res <- summary(mix_some_cured, t=Inf, type = "cumhaz", tidy=T)
  expect_equal(as.numeric(mix_some_cured_res)[2], -log(as.numeric(mix_some_cured$res[1])))

  # NON-MIXTURE MODELS
  # Cumulaitve hazard at t = Inf should equal -log(cure fraction0
  nmix_some_cured <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="identity", dist="weibull", mixture=F)
  nmix_some_cured_res <- summary(nmix_some_cured, t=Inf, type = "cumhaz", tidy=T)
  expect_equal(as.numeric(nmix_some_cured_res)[2], -log(as.numeric(nmix_some_cured$res[1])))
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


