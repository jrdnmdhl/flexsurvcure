##' Non-Mixture Cure Models
##'
##' Probability density, distribution, quantile, random generation, hazard
##' cumulative hazard, mean, and restricted mean functions for generic
##' non-mixture cure models.  These distribution functions take as arguments
##' the corresponding functions of the base distribution used.
##'
##' es dnmixsurv pnmixsurv qnmixsurv rnmixsurv
##' hnmixsurv Hnmixsurv mean_nmixsurv rmst_nmixsurv
##' @param pfun The base distribution's cumulative distribution function.
##' @param dfun The base distribution's probability density function.
##' @param x,q,t Vector of times.
##' @param p Vector of probabilities.
##' @param n Number of random numbers to simulate.
##' @param theta The estimated cure fraction.
##' @param ... Parameters to be passed to the pdf or cdf of the base
##' distribution.
##' @return \code{dnmixsurv} gives the density, \code{pnmixsurv} gives the
##' distribution function, \code{hnmixsurv} gives the hazard and
##' \code{Hnmixsurv} gives the cumulative hazard.
##'
##' \code{qnmixsurv} gives the quantile function, which is computed by crude
##' numerical inversion.
##'
##' \code{rnmixsurv} generates random survival times by using \code{qnmixsurv}
##' on a sample of uniform random numbers.  Due to the numerical root-finding
##' involved in \code{qnmixsurv}, it is slow compared to typical random number
##' generation functions.
##'
##' \code{mean_nmixsurv} and \code{rmst_nmixsurv} give the mean and restricted
##' mean survival times, respectively.
##' @author Jordan Amdahl <jrdnmdhl@gmail.com>
##' @keywords distribution
##' @name nmixsurv
NULL

##' @export
##' @rdname nmixsurv
pnmixsurv = function(pfun, q, theta, ...) {
  dots <- list(...)
  args <- dots
  args$lower.tail <- T
  args$log.p <- F
  out <- theta ^ do.call(pfun, append(list(q), args))
  if (is.null(dots$lower.tail) || dots$lower.tail) {
    out <- 1 - out
  }
  if (!is.null(dots$log.p) && dots$log.p) {
    out <- log(out)
  }
  return(out)
}

##' @export
##' @rdname nmixsurv
hnmixsurv = function(dfun,x, theta, ...) {
  dots <- list(...)
  args <- dots
  args$log <- F
  u_pdf <-
  out <- -log(theta) * do.call(dfun, append(list(x), args))
  if (!is.null(dots$log) && dots$log) {
    out <- log(out)
  }
  return(out)
}

##' @export
##' @rdname nmixsurv
Hnmixsurv = function(pfun, x, theta, ...) {
  dots <- list(...)
  pargs <- dots
  pargs$lower.tail <- F
  pargs$log.p <- F
  pargs$log <- NULL
  surv <- do.call(pnmixsurv, append(list(pfun, x), pargs))
  out <- -log(surv)
  if (!is.null(dots$log) && dots$log) {
    out <- log(out)
  }
  return(out)
}

##' @export
##' @rdname nmixsurv
dnmixsurv = function(dfun, pfun, x, theta, ...) {
  dots <- list(...)
  pargs <- dots
  pargs$lower.tail <- F
  pargs$log.p <- F
  pargs$log <- NULL
  hargs <- dots
  hargs$log <- F
  u_surv <- do.call(pnmixsurv, append(list(pfun, x, theta), pargs))
  u_haz <- do.call(hnmixsurv, append(list(dfun, x, theta), hargs))
  out <- u_surv * u_haz
  if (!is.null(dots$log) && dots$log) {
    out <- log(out)
  }
  return(out)
}

##' @export
##' @rdname nmixsurv
qnmixsurv = function(pfun, p, theta, ...) {
  dots <- list(...)
  args <- dots
  args$lower.tail <- F
  args$log.p <- F
  if (dots$log.p) p <- exp(p)
  if (!dots$lower.tail) p <- 1 - p
  out <- do.call(
    qgeneric,
    append(
      list(
        function(...) pnmixsurv(pfun, ...),
        p = p,
        theta = theta
      ),
      args
    )
  )
  return(out)
}


##' @export
##' @rdname nmixsurv
rnmixsurv = function(pfun, n, theta, ...) {
  dots <- list(...)
  args <- dots
  args$lower.tail <- F
  args$log.p <- F
  if (dots$log.p) p <- exp(p)
  if (!dots$lower.tail) p <- 1 - p
  out <- do.call(
    qgeneric,
    append(
      list(
        function(...) pnmixsurv(pfun, ...),
        p = runif(n),
        theta = theta
      ),
      args
    )
  )
  return(out)
}

##' @export
##' @rdname nmixsurv
rmst_nmixsurv = function(pfun, t, theta, ...) {
  args <- list(...)
  out <- do.call(
    rmst_generic,
    append(
      list(
        function(...) pnmixsurv(pfun, ...),
        t = t,
        theta = theta
      ),
      args
    )
  )
  return(out)
}

##' @export
##' @rdname nmixsurv
mean_nmixsurv = function(pfun, theta, ...) {
  if(theta > 0) {
    out <- Inf
  }else {
    args <- list(...)
    out <- do.call(
      rmst_generic,
      append(
        list(
          pfun,
          t = Inf,
          start = 0
        ),
        args
      )
    )
  }
  return(out)
}

