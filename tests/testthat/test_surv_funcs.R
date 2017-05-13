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
  t_rmst <- seq(from=0, to=10000, by=2500)
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
