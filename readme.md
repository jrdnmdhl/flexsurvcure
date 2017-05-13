flexsurvcure
============

The development repository for flexsurvcure, an R package for parametric mixture and non-mixture cure models.  Flexsurvcure is based on [flexsurv](http://cran.r-project.org/package=flexsurv), the R package for parametric survival modelling.

## Installation (development version)

```r
install.packages("flexsurv")
install.packages("devtools")
devtools::install_github('jrdnmdhl/flexsurvcure')
```

## Supported Distributions

All of the built-in distribution with flexsurvreg are supported, though some currently have issues with convergence and numerical instability.  The following distributions currently seem reliable:

- Exponential (exp)
- Weibull (weibull, weibullPH)
- Lognormal (lnorm)
- Log-Logistic (llogis)

The following distributions are supported, but may not be reliable:

- Gompertz (gompertz)
- Generalized Gamma (gengamma, gengamma.orig)
- Generalized F (genf, genf.orig)

Custom distributions can also be used by passing a distribution list (see [flexsurv examples](https://cran.r-project.org/web/packages/flexsurv/vignettes/flexsurv-examples.pdf)).

## Fitting a mixture cure model
```r
mixture = flexsurvcure(Surv(rectime,censrec)~group, data=bc, dist="weibullPH", link="logistic", mixture = T)
plot(mixture)
```


## Fitting a non-mixture cure model
```r
non_mixture = flexsurvcure(Surv(rectime,censrec)~group, data=bc, dist="weibullPH", link="loglog", mixture = F)
plot(non_mixture)
```

## Covariates on parameters other than cure fraction
```r
non_mixture_covarite_scale = flexsurvcure(Surv(rectime,censrec)~group, data=bc, anc=list(scale=~group), dist="weibullPH", link="loglog", mixture = F)
plot(non_mixture_covarite_scale)
```
