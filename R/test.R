pcuremix <- function(theta, surv) {
  theta + (1 - theta) * surv
}

#' @export
pweibullmix <- function(q, theta, shape, scale = 1, lower.tail = T, log.p = F) {
  out <- pcuremix(
    theta,
    flexsurv::pweibullPH(q, shape = shape, scale = scale, lower.tail = F)
  )
  if (lower.tail) {
    out <- 1 - out
  }
  if (log.p) {
    out <- log(out)
  }
  return(out)
}

#' @export
pweibullmix.quiet <- function(q, theta, shape, scale = 1, lower.tail = T, log.p = F) {
  out <- suppressWarnings(
    pweibullmix(q, theta, shape, scale = scale, lower.tail = lower.tail, log.p = log.p)
  )
  return(out)
}


pcurenmix <- function(theta, surv) {
  theta ^ (1 - surv)
}

hcuremix <- function(theta, dens, surv) {
  ((1 - theta) * dens) / (theta + (1 - theta) * surv)
}

#' @export
hweibullmix <- function(x, theta, shape, scale = 1, log = F) {
  out <- hcuremix(
    theta,
    flexsurv::dweibullPH(x, shape = shape, scale = scale),
    flexsurv::pweibullPH(x, shape = shape, scale = scale, lower.tail = F)
  )
  if (log) {
    out <- log(out)
  }
  return(out)
}

#' @export
hweibullmix.quiet <- function(x, theta, shape, scale = 1, log = F) {
  out <- suppressWarnings(
    hweibullmix(x, theta, shape, scale = scale, log = log)
  )
  return(out)
}

hcurnemix <- function(theta, haz, surv) {
  -log(theta) * surv * haz
}

dcuremix <- function(theta, haz, surv) {
  pcuremix(theta, surv) * hcuremix(theta, haz, surv)
}

#' @export
dweibullmix <- function(x, theta, shape, scale = 1, log = F) {
  out <- dcuremix(
    theta,
    flexsurv::dweibullPH(x, shape = shape, scale = scale),
    flexsurv::pweibullPH(x, shape = shape, scale = scale, lower.tail = F)
  )
  if (log) {
    out <- log(out)
  }
  return(out)
}

#' @export
dweibullmix.quiet <- function(x, theta, shape, scale = 1, log = F) {
  out <- suppressWarnings(
    dweibullmix(x, theta, shape, scale = scale, log = log)
  )
  return(out)
}



dcurenmix <- function(theta, haz, surv) {
  out <- pcurenmix(theta, surv) * hcurenmix(theta, haz, surv)
  return(out)
}

#' @export
qweibullmix <- function(p, theta, shape, scale = 1, lower.tail = T, log.p = F) {
  if (log.p) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  qgeneric(
    pweibullmix,
    p = p,
    theta = theta,
    shape = shape,
    scale = scale
  )
}

#' @export
qweibullmix.quiet <- function(p, theta, shape, scale = 1, lower.tail = T, log.p = F) {
  out <- suppressWarnings(
    pweibullmix.quiet(p, theta, shape, scale = scale, lower.tail = lower.tail, log.p = log.p)
  )
  return(out)
}


flexsurvcure.dists = list(
  weibull = list(
    name = "weibullmix.quiet",
    pars = c("theta", "shape", "scale"),
    location = "theta",
    transforms = c(gtools::logit, log, log),
    inv.transforms = c(gtools::inv.logit, exp, exp),
    inits = function(t, mf) {
      surv <- as.matrix(mf[ ,1])
      weights <- mf[ ,ncol(mf)]
      events <- surv[surv[ ,2] == 1, ]
      theta <- 1 - mean(surv[ ,2] * weights)
      shape <- 1
      scale <- 1 / mean(events[ ,1])
      out <- c(theta, shape, scale)
      return(out)
    }
  )
)

