expand.inits.args <- function(inits){
  inits2 <- inits
  formals(inits2) <- alist(t=,mf=,mml=,aux=)
  body(inits2) <- body(inits)
  inits2
}

#' @export
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

  # Patch the transformations based on link argument
  dist_list <- flexsurv.dists[[dist]]
  dist_list$name <- paste0(dist_list$name, "_mix")
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
  } else {
    stop("Link must be 'logistic', 'loglog', or 'identity'")
  }

  base_init <- expand.inits.args(dist_list$inits)
  if(dist == "weibullPHdfdf") {
    dist_list$inits <- function(t, mf) {
      surv <- as.matrix(mf[ ,1])
      weights <- mf[ ,ncol(mf)]
      events <- surv[surv[ ,2] == 1, ]
      theta <- 1 - mean(surv[ ,2] * weights)
      shape <- 1
      scale <- 1 / mean(events[ ,1])
      out <- c(theta, shape, scale)
      print(out)
      return(out)
    }
  } else {
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
  }

  # Build function list
  pfun = get(paste0("p", dist))
  dfun = get(paste0("d", dist))
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

  # Generate fit
  out <- do.call(
    "flexsurvreg",
    list(
      formula,
      data = temp$data,
      weights = temp$weights,
      subset = temp$subset,
      bhazard = temp$bhazard,
      dist = dist_list,
      dfns = dfns_list,
      optim = list(
        maxit = 1000
      ),
      ...
    )
  )

  # Use top-level call and set additional properties/attributes
  out$call <- call
  class(out) <- c("flexsurvcure", class(out))
  out$link <- link
  out
}

