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
  logistic_params <- c( -.4731367, .448163, -11.10519)
  logistic_lower <- c( -.7948533, .3351858, -12.2618)
  logistic_upper <- c(-.1514202, .5611403,  -9.948572)
  logistic_se <- c(.1641441 , .0576425 ,  .5901206)
  logistic_null <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="weibullPH")

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
  ident_params <- c(.3838745, .4481635, -11.10519)
  ident_lower <- c(.3077839, .3351859, -12.26181)
  ident_upper <- c(.4599651, .561141,  -9.948572)
  ident_se <- c(.0388224, .0576427,  .5901226)
  ident_null <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="identity", dist="weibullPH")

  # Tolerance here is very permissive
  expect_equal(ident_params, unname(ident_null$res.t[ ,1]), tolerance=1e-2)
  expect_equal(ident_lower, unname(ident_null$res.t[ ,2]), tolerance=1e-2)
  expect_equal(ident_upper, unname(ident_null$res.t[ ,3]), tolerance=1e-2)
  expect_equal(ident_se, unname(ident_null$res.t[ ,4]), tolerance=1e-2)

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
  loglog_params <- c(.2484884, 6.982434, -.0933944)
  loglog_lower <- c( -.1391261, 6.726824, -.2503376)
  loglog_upper <- c(.6361029, 7.238045, .0635487)
  loglog_se <- c(.1977661, .1304159, .0800745)
  loglog_null <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="loglog", dist="lnorm")

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
  ident_params <- c(.2774585, 6.982433, -.093395)
  ident_lower <- c(.1395744, 6.726823, -.2503383)
  ident_upper <- c(.4153425, 7.238044,  .0635482)
  ident_se <- c(.0703503, .1304159,  .0800746)
  ident_null <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="identity", dist="lnorm")

  # Tolerance here is very permissive
  expect_equal(ident_params, unname(ident_null$res.t[ ,1]), tolerance=1e-2)
  expect_equal(ident_lower, unname(ident_null$res.t[ ,2]), tolerance=1e-2)
  expect_equal(ident_upper, unname(ident_null$res.t[ ,3]), tolerance=1e-2)
  expect_equal(ident_se, unname(ident_null$res.t[ ,4]), tolerance=1e-2)

})
test_that("Weibull Non-Mixture matches stata", {

  # Weibull, Logistic
  # ------------------------------------------------------------------------------
  #   _t |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
  # -------------+----------------------------------------------------------------
  #   pi           |
  #   _cons |  -.5151878   .1806544    -2.85   0.004    -.8692639   -.1611116
  # -------------+----------------------------------------------------------------
  #   ln_lambda    |
  #   _cons |  -12.16467   .6119361   -19.88   0.000    -13.36405    -10.9653
  # -------------+----------------------------------------------------------------
  #   ln_gamma     |
  #   _cons |   .5110451   .0585201     8.73   0.000     .3963477    .6257425
  # ------------------------------------------------------------------------------
  logistic_params <- c( -.5151878, .5110451, -12.16467)
  logistic_lower <- c( -.8692639, .3963477, -13.36405)
  logistic_upper <- c(-.1611116, .6257425,  -10.9653)
  logistic_se <- c(.1806544 , .0585201 ,  .6119361)
  logistic_null <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="logistic", dist="weibullPH", mixture=F)

  # Tolerance here is very permissive
  expect_equal(logistic_params, unname(logistic_null$res.t[ ,1]), tolerance=1e-2)
  expect_equal(logistic_lower, unname(logistic_null$res.t[ ,2]), tolerance=1e-2)
  expect_equal(logistic_upper, unname(logistic_null$res.t[ ,3]), tolerance=1e-2)
  expect_equal(logistic_se, unname(logistic_null$res.t[ ,4]), tolerance=1e-2)

  # Weibull, Loglog
  # ------------------------------------------------------------------------------
  #   _t |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
  # -------------+----------------------------------------------------------------
  #   pi           |
  #   _cons |  -.0165781   .1149836    -0.14   0.885    -.2419419    .2087857
  # -------------+----------------------------------------------------------------
  #   ln_lambda    |
  #   _cons |  -12.16465   .6119389   -19.88   0.000    -13.36402   -10.96527
  # -------------+----------------------------------------------------------------
  #   ln_gamma     |
  #   _cons |   .5110431   .0585205     8.73   0.000     .3963451    .6257412
  # ------------------------------------------------------------------------------

  loglog_params <- c(-.0165781, .5110431, -12.16465)
  loglog_lower <- c( -.2419419, .3963451, -13.36402)
  loglog_upper <- c(.2087857, .6257412, -10.96527)
  loglog_se <- c(.1149836, .0585205, .6119389)
  loglog_null <- flexsurvcure(Surv(rectime, censrec)~1,data=bc,link="loglog", dist="weibullPH", mixture=F)

  # Tolerance here is very permissive
  expect_equal(loglog_params, unname(loglog_null$res.t[ ,1]), tolerance=1e-2)
  expect_equal(loglog_lower, unname(loglog_null$res.t[ ,2]), tolerance=1e-2)
  expect_equal(loglog_upper, unname(loglog_null$res.t[ ,3]), tolerance=1e-2)
  expect_equal(loglog_se, unname(loglog_null$res.t[ ,4]), tolerance=1e-2)
})
test_that("Weibull Mixture w/ covariate matches stata", {

  good_data <- bc
  good_data$good <- ifelse(good_data$group == "Good", 1, 0)

  # Weibull, loglog
  # ------------------------------------------------------------------------------
  #   _t |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
  # -------------+----------------------------------------------------------------
  #   pi           |
  #   good |  -1.304919   .1775218    -7.35   0.000    -1.652856   -.9569831
  # _cons |    .360375   .1384624     2.60   0.009     .0889936    .6317564
  # -------------+----------------------------------------------------------------
  #   ln_lambda    |
  #   _cons |  -11.05986   .5894734   -18.76   0.000     -12.2152    -9.90451
  # -------------+----------------------------------------------------------------
  #   ln_gamma     |
  #   _cons |   .4434316   .0580233     7.64   0.000      .329708    .5571553
  # ------------------------------------------------------------------------------

  loglog_params <- c(.360375, .4434316, -11.05986, -1.304919)
  loglog_lower <- c(.0889936, .329708, -12.2152, -1.652856)
  loglog_upper <- c(.6317564, .5571553,-9.90451, -.9569831)
  loglog_se <- c(.1384624, .0580233, .5894734, .1775218)
  loglog_cov <- flexsurvcure(Surv(rectime, censrec)~good,data=good_data,link="loglog", dist="weibullPH")

  # Tolerance here is very permissive
  expect_equal(loglog_params, unname(loglog_cov$res.t[ ,1]), tolerance=1e-2)
  expect_equal(loglog_lower, unname(loglog_cov$res.t[ ,2]), tolerance=1e-2)
  expect_equal(loglog_upper, unname(loglog_cov$res.t[ ,3]), tolerance=1e-2)
  expect_equal(loglog_se, unname(loglog_cov$res.t[ ,4]), tolerance=1e-2)

})


test_that("Weibull Mixture w/ baseline hazard", {
  # Weibull, Logistic
  # ------------------------------------------------------------------------------
  #   _t |      Coef.   Std. Err.      z    P>|z|     [95% Conf. Interval]
  # -------------+----------------------------------------------------------------
  #   pi           |
  #   _cons |  -.4589729   .1644206    -2.79   0.005    -.7812314   -.1367144
  # -------------+----------------------------------------------------------------
  #   ln_lambda    |
  #   _cons |  -1.871874   .1175981   -15.92   0.000    -2.102362   -1.641386
  # -------------+----------------------------------------------------------------
  #   ln_gamma     |
  #   _cons |   .4501334   .0579905     7.76   0.000     .3364742    .5637926
  # ------------------------------------------------------------------------------
  #
  logistic_params <- c( -.4589729, .4501334, -1.871874)
  logistic_lower <- c( -.7812314, .3364742, -2.102362)
  logistic_upper <- c(-.1367144, .5637926,  -1.641386)
  logistic_se <- c(.1644206 , .0579905 ,  .1175981)
  logistic_null <- flexsurvcure(Surv(recyrs, censrec)~1,data=bc,link="logistic", dist="weibullPH", bhazard=rep(0.001,686))

  # Tolerance here is very permissive
  expect_equal(logistic_params, unname(logistic_null$res.t[ ,1]), tolerance=1e-2)
  expect_equal(logistic_lower, unname(logistic_null$res.t[ ,2]), tolerance=1e-2)
  expect_equal(logistic_upper, unname(logistic_null$res.t[ ,3]), tolerance=1e-2)
  expect_equal(logistic_se, unname(logistic_null$res.t[ ,4]), tolerance=1e-2)

})
