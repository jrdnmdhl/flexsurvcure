flexsurvcure
============

The development repository for flexsurvcure, an R package for parametric mixture and non-mixture cure models.  Flexsurvcure is based on [flexsurv](http://cran.r-project.org/package=flexsurv), the R package for parametric survival modelling.

## Installation (development version)

```r
install.packages("flexsurv")
install.packages("devtools") # if devtools not already installed
install.packages("flexsurv")
devtools::install_github('jrdnmdhl/flexsurvcure')
```

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
