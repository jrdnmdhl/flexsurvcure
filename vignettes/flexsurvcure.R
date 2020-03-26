## ---- warning=FALSE, message=FALSE---------------------------------------
library(flexsurvcure)
cure_model <- flexsurvcure(Surv(rectime, censrec)~group, data=bc, link="logistic", dist="weibullPH", mixture=T)
print(cure_model)

## ---- warning=FALSE, message=FALSE---------------------------------------
plot(cure_model)

## ---- warning=FALSE, message=FALSE---------------------------------------
summary(cure_model, t=seq(from=0,to=3000,by=1000), type="survival", tidy=T)

## ---- warning=FALSE, message=FALSE---------------------------------------
cure_model_complex <- flexsurvcure(Surv(rectime, censrec)~group, data=bc, link="logistic", dist="weibullPH", mixture=T, anc=list(scale=~group))
print(cure_model_complex)
plot(cure_model_complex)

## ---- warning=FALSE, message=FALSE---------------------------------------
library(flexsurvcure)
cure_model_nmix <- flexsurvcure(Surv(rectime, censrec)~group, data=bc, link="loglog", dist="weibullPH", mixture=F)
print(cure_model_nmix)

