test_that("Weibull Mixture matches stata", {

  # Weibull, Logistic
  # ------------------------------------------------------------------------------
  #   _t |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
  # -------------+----------------------------------------------------------------
  #   pi           |
  #   _cons |  -.4731367   .1641441    -2.88   0.004    -.7948533   -.1514202
  # -------------+----------------------------------------------------------------
  #   ln_lambda    |
  #   _cons |  -11.10519   .5901206   -18.82   0.000     -12.2618   -9.948572
  # -------------+----------------------------------------------------------------
  #   ln_gamma     |
  #   _cons |    .448163   .0576425     7.77   0.000     .3351858    .5611403
  # ------------------------------------------------------------------------------
  logistic_params = c( -.4731367, .448163, -11.10519)
  logistic_lower = c( -.7948533, .3351858, -12.2618)
  logistic_upper = c(-.1514202, .5611403,  -9.948572)
  logistic_se = c(.1641441 , .0576425 ,  .5901206)
  logistic_null = flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="weibullPH")

  # Tolerance here is very permissive
  expect_equal(logistic_params, unname(logistic_null$res.t[ ,1]), tolerance=1e-2)
  expect_equal(logistic_lower, unname(logistic_null$res.t[ ,2]), tolerance=1e-2)
  expect_equal(logistic_upper, unname(logistic_null$res.t[ ,3]), tolerance=1e-2)
  expect_equal(logistic_se, unname(logistic_null$res.t[ ,4]), tolerance=1e-2)

  # Weibull, Identity Link
  # _t |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
  # -------------+----------------------------------------------------------------
  #   pi           |
  #   _cons |   .3838745   .0388224     9.89   0.000     .3077839    .4599651
  # -------------+----------------------------------------------------------------
  #   ln_lambda    |
  #   _cons |  -11.10519   .5901226   -18.82   0.000    -12.26181   -9.948572
  # -------------+----------------------------------------------------------------
  #   ln_gamma     |
  #   _cons |   .4481635   .0576427     7.77   0.000     .3351859     .561141
  # ------------------------------------------------------------------------------
  ident_params = c(.3838745, .4481635, -11.10519)
  ident_lower = c(.3077839, .3351859, -12.26181)
  ident_upper = c(.4599651, .561141,  -9.948572)
  ident_se = c(.0388224, .0576427,  .5901226)
  ident_null = flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="identity", dist="weibullPH")

  # Tolerance here is very permissive
  expect_equal(ident_params, unname(ident_null$res.t[ ,1]), tolerance=1e-2)
  expect_equal(ident_lower, unname(ident_null$res.t[ ,2]), tolerance=1e-2)
  expect_equal(ident_upper, unname(ident_null$res.t[ ,3]), tolerance=1e-2)
  expect_equal(ident_se, unname(ident_null$res.t[ ,4]), tolerance=1e-2)



  #loglog_null = flexsurvcure(survival::Surv(rectime, censrec)~1,data=flexsurv::bc,link="loglog", dist="weibull")
})
test_that("Lognormal Mixture matches stata", {


  # Lognormal, Loglog
  # ------------------------------------------------------------------------------
  #   _t |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
  # -------------+----------------------------------------------------------------
  #   pi           |
  #   _cons |   .2484884   .1977661     1.26   0.209    -.1391261    .6361029
  # -------------+----------------------------------------------------------------
  #   mu           |
  #   _cons |   6.982434   .1304159    53.54   0.000     6.726824    7.238045
  # -------------+----------------------------------------------------------------
  #   ln_sigma     |
  #   _cons |  -.0933944   .0800745    -1.17   0.243    -.2503376    .0635487
  # ------------------------------------------------------------------------------
  loglog_params = c(.2484884, 6.982434, -.0933944)
  loglog_lower = c( -.1391261, 6.726824, -.2503376)
  loglog_upper = c(.6361029, 7.238045, .0635487)
  loglog_se = c(.1977661, .1304159, .0800745)
  loglog_null = flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="loglog", dist="lnorm")

  # Tolerance here is very permissive
  expect_equal(loglog_params, unname(loglog_null$res.t[ ,1]), tolerance=1e-2)
  expect_equal(loglog_lower, unname(loglog_null$res.t[ ,2]), tolerance=1e-2)
  expect_equal(loglog_upper, unname(loglog_null$res.t[ ,3]), tolerance=1e-2)
  expect_equal(loglog_se, unname(loglog_null$res.t[ ,4]), tolerance=1e-2)

  # Lognormal, Identity Link
  # ------------------------------------------------------------------------------
  #   _t |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
  # -------------+----------------------------------------------------------------
  #   pi           |
  #   _cons |   .2774585   .0703503     3.94   0.000     .1395744    .4153425
  # -------------+----------------------------------------------------------------
  #   mu           |
  #   _cons |   6.982433   .1304159    53.54   0.000     6.726823    7.238044
  # -------------+----------------------------------------------------------------
  #   ln_sigma     |
  #   _cons |   -.093395   .0800746    -1.17   0.243    -.2503383    .0635482
  # ------------------------------------------------------------------------------
  ident_params = c(.2774585, 6.982433, -.093395)
  ident_lower = c(.1395744, 6.726823, -.2503383)
  ident_upper = c(.4153425, 7.238044,  .0635482)
  ident_se = c(.0703503, .1304159,  .0800746)
  ident_null = flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="identity", dist="lnorm")

  # Tolerance here is very permissive
  expect_equal(ident_params, unname(ident_null$res.t[ ,1]), tolerance=1e-2)
  expect_equal(ident_lower, unname(ident_null$res.t[ ,2]), tolerance=1e-2)
  expect_equal(ident_upper, unname(ident_null$res.t[ ,3]), tolerance=1e-2)
  expect_equal(ident_se, unname(ident_null$res.t[ ,4]), tolerance=1e-2)

})
