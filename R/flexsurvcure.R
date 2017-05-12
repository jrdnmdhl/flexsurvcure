
# Taken from flexsurv, needed to wrap init functions of base distributions

expand.inits.args <- function(inits) {
  inits2 <- inits
  formals(inits2) <- alist(t=,mf=,mml=,aux=)
  body(inits2) <- body(inits)
  inits2
}

#' Mixture and Non-Mixture Parametric Cure Models
##'
##' Mixture and non-mixture cure models using flexible base distribtuions
##' from the flexsurv package.
##'
##' This function works as a wrapper around \code{flexsurv::flexsurvreg} by
##' dynamically constructing a custom distribution using
##' \code{\link{dsurvmix}}, \code{\link{psurvmix}}}.
##'
##' In a parametric mixture model, it is assume that there is a group of individuals
##' who experience no excess mortality, with the proportion of patients being given
##' by the estimated cure fraction, and a parametric distribution representing excess
##' mortality for the remaining individuals.
##'
##' By contrast, a parametric non-mixture model simply rescales an existing parametric
##' distribution such that survival reaches an asyptote >0.
##'
##' @param formula A formula expression in conventional R linear modelling
##' syntax. The response must be a survival object as returned by the
##' \code{\link{Surv}} function, and any covariates are given on the right-hand
##' side.  For example,
##'
##' \code{Surv(time, dead) ~ age + sex}
##'
##' \code{Surv} objects of \code{type="right"},\code{"counting"},
##' \code{"interval1"} or \code{"interval2"} are supported, corresponding to
##' right-censored, left-truncated or interval-censored observations.
##'
##' If there are no covariates, specify \code{1} on the right hand side, for
##' example \code{Surv(time, dead) ~ 1}.
##'
##' By default, covariates are placed on the ``theta'' parameter of the
##' distribution, representing the cure fraction, through a linear
##' model with the selected link function.
##'
##' Covariates can be placed on parameters of the base distribution by using the
##' name of the parameter as a ``function'' in the formula.  For example, in a
##' Weibull model, the following expresses the scale parameter in terms of age
##' and a treatment variable \code{treat}, and the shape parameter in terms of
##' sex and treatment.
##'
##' \code{Surv(time, dead) ~ age + treat + shape(sex) + shape(treat)}
##'
##' However, if the names of the ancillary parameters clash with any real
##' functions that might be used in formulae (such as \code{I()}, or
##' \code{factor()}), then those functions will not work in the formula.  A
##' safer way to model covariates on ancillary parameters is through the
##' \code{anc} argument to \code{flexsurv::flexsurvreg}.
##'
##' \code{\link{survreg}} users should also note that the function
##' \code{strata()} is ignored, so that any covariates surrounded by
##' \code{strata()} are applied to the location parameter.
##' @param anc An alternative and safer way to model covariates on ancillary
##' parameters, that is, parameters other than the cure fraction.  This is a
##' named list of formulae, with the name of each component giving the parameter
##' to be modelled.  The model above can also be defined as:
##'
##' \code{Surv(time, dead) ~ age + treat, anc = list(shape = ~ sex + treat)}
##' @param data A data frame in which to find variables supplied in
##' \code{formula}.  If not given, the variables should be in the working
##' environment.
##' @param weights Optional variable giving case weights.
##' @param bhazard Optional variable giving expected hazards for relative
##' survival models.
##' @param subset Vector of integers or logicals specifying the subset of the
##' observations to be used in the fit.
##' @param na.action a missing-data filter function, applied after any 'subset'
##' argument has been used. Default is \code{options()$na.action}.
##' @param dist A string representing one of the built-in distributions of flexsurv.
##' @param link A string representing the link function to use for estimation of the
##' cure fraction.  Defaults to logistic.
##' @param mixture optional TRUE/FALSE to specify whether a mixture model should be fitted.  Defaults to TRUE.
##' @export
flexsurvcure <- function(formula, data, weights, bhazard, subset, dist, link = "logistic", mixture = T, ...) {
  call <- match.call()
  indx <- match(c("formula", "data", "weights", "bhazard", "subset", "na.action"), names(call), nomatch = 0)
  if (indx[1] == 0)
    stop("A \"formula\" argument is required")
  temp <- call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")
  if (missing(data)) temp[["data"]] <- environment(formula)
  if (missing(data)) data <- environment(formula)
  if (missing(dist)) stop("Must provide dist")
  optim = list()

  # Patch the transformations based on link argument
  dist_list <- flexsurv.dists[[dist]]
  dist_list$name <- paste0(dist_list$name, "_mix")
  n_base_par <- length(dist_list$pars)
  dist_list$pars <- c("theta", dist_list$pars)
  dist_list$location <- "theta"
  if(is.null(dist_list)) stop("Distribution not found")
  if (link == "logistic") {
    dist_list$transforms <- append(list(logit), dist_list$transforms)
    dist_list$inv.transforms <- append(list(inv.logit), dist_list$inv.transforms)
  } else if(link == "loglog") {
    dist_list$transforms <- append(list(function(x) log(-log(x))), dist_list$transforms)
    dist_list$inv.transforms <- append(list(function(x) exp(-exp(x))), dist_list$inv.transforms)
  } else if(link == "identity") {
    dist_list$transforms <- append(list(identity), dist_list$transforms)
    dist_list$inv.transforms <- append(list(identity), dist_list$inv.transforms)
    optim$method <- "L-BFGS-B"
    optim$lower = c(0, rep(-Inf, n_base_par))
    optim$upper = c(1, rep(Inf, n_base_par))
  } else {
    stop("Link must be 'logistic', 'loglog', or 'identity'")
  }

  base_init <- expand.inits.args(dist_list$inits)

  dist_list$inits <- function(t, mf, mml, aux) {
    # To estimate initial values:
    # -Cure fraction based on minimum KM survival
    # -Other parameters based on normal initial values
    #  run only on events.
    surv <- as.matrix(mf[ ,1])
    weights <- mf[ ,ncol(mf)]
    selector <- surv[ ,2] == 1
    aux_sf <- list(
      formula = aux$forms[[1]],
      data = aux$data,
      weights = aux$weights
    )
    sf <- do.call(survfit, aux_sf)
    # Can't allow value of 0
    theta = max(min(sf$surv), 0.01)
    aux_events <- aux
    aux_events$data = aux$data[selector, ]
    out <- c(theta, base_init(t=t[selector], mf=mf[selector, ], mml=mml[selector, ], aux=aux_events))
    return(out)
  }

  # Build function list
  pfun = get(paste0("p", dist))
  dfun = get(paste0("d", dist))
  if(mixture) {
    dfns_list = list(
      p = function(q, ...) {
        pmixsurv(pfun, q, ...)
      },
      d = function(x, ...) {
        dmixsurv(dfun, pfun, x, ...)
      },
      h = function(x, ...) hmixsurv(dfun, pfun, x, ...),
      q = function(p, ...) qmixsurv(pfun, p, ...),
      mean = function(t, ...) rmst_mixsurv(pfun, t, ...),
      rmst = function(t, ...) mean_mixsurv(pfun, t, ...)
    )
  } else {
    dfns_list = list(
      p = function(q, ...) {
        pnmixsurv(pfun, q, ...)
      },
      d = function(x, ...) {
        dnmixsurv(dfun, pfun, x, ...)
      },
      h = function(x, ...) hnmixsurv(dfun, x, ...),
      q = function(p, ...) qnmixsurv(pfun, p, ...),
      mean = function(t, ...) rmst_nmixsurv(pfun, t, ...),
      rmst = function(t, ...) mean_nmixsurv(pfun, t, ...)
    )
  }

  # Generate fit
  out <- do.call(
    "flexsurvreg",
    append(
      list(
        formula,
        data = temp$data,
        weights = temp$weights,
        subset = temp$subset,
        bhazard = temp$bhazard,
        dist = dist_list,
        dfns = dfns_list,
        ...
      ),
      optim
    )
  )

  # Use top-level call and set additional properties/attributes
  out$call <- call
  class(out) <- c("flexsurvcure", class(out))
  out$link <- link
  out$mixture <- mixture
  out
}

