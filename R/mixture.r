##' Mixture cure models
##'
##' Probability density, distribution, quantile, random generation, hazard
##' cumulative hazard, mean, and restricted mean functions for generic
##' mixture cure models.  These distribution functions take as arguments
##' the corresponding functions of the base distribution used.
##'
##' @aliases dmixsurv pmixsurv qmixsurv rmixsurv
##' hmixsurv Hmixsurv mean_mixsurv rmst_mixsurv
##' @param pfun The base distribution's cumulative distribution function.
##' @param dfun The base distribution's probability density function.
##' @param qfun The base distribution's quantile function.
##' @param x,q,t Vector of times.
##' @param p Vector of probabilities.
##' @param n Number of random numbers to simulate.
##' @param theta The estimated cure fraction.
##' @param ... additional parameters to be passed to the pdf or cdf of the base
##' distribution.
##' @return \code{dmixsurv} gives the density, \code{pmixsurv} gives the
##' distribution function, \code{hmixsurv} gives the hazard and
##' \code{Hmixsurv} gives the cumulative hazard.
##'
##' \code{qmixsurv} gives the quantile function, which is computed by crude
##' numerical inversion.
##'
##' \code{rmixsurv} generates random survival times by using \code{qmixsurv}
##' on a sample of uniform random numbers.  Due to the numerical root-finding
##' involved in \code{qmixsurv}, it is slow compared to typical random number
##' generation functions.
##'
##' \code{mean_mixsurv} and \code{rmst_mixsurv} give the mean and restricted
##' mean survival times, respectively.
##' @author Jordan Amdahl <jrdnmdhl@gmail.com>
##' @keywords distribution
##' @name mixsurv
NULL

##' @export
##' @rdname mixsurv
pmixsurv = function(pfun, q, theta, ...) {
  dots <- list(...)
  args <- dots
  args$lower.tail <- F
  args$log.p <- F
  out <- theta + (1 - theta) * do.call(pfun, append(list(q), args))
  if (is.null(dots$lower.tail) || dots$lower.tail) {
    pos_inf <- is.infinite(q) & (q > 0)
    out[pos_inf] <- 0
    out <- 1 - out
  }
  if (!is.null(dots$log.p) && dots$log.p) {
    out <- log(out)
  }
  return(out)
}

##' @export
##' @rdname mixsurv
hmixsurv = function(dfun, pfun, x, theta, ...) {
  dots <- list(...)
  pargs <- dots
  pargs$lower.tail <- F
  pargs$log.p <- F
  pargs$log <- NULL
  dargs <- dots
  dargs$log <- F
  u_surv <- do.call(pfun, append(list(x), pargs))
  u_pdf <- do.call(dfun, append(list(x), dargs))
  out <- ((1 - theta) * u_pdf) / (theta + (1 - theta) * u_surv)
  if (!is.null(dots$log) && dots$log) {
    out <- log(out)
  }
  return(out)
}

##' @export
##' @rdname mixsurv
Hmixsurv = function(pfun, x, theta, ...) {
  dots <- list(...)
  pargs <- dots
  pargs$lower.tail <- F
  pargs$log.p <- F
  pargs$log <- NULL
  surv <- do.call(pmixsurv, append(list(pfun, x, theta), pargs))
  out <- -log(surv)
  if (!is.null(dots$log) && dots$log) {
    out <- log(out)
  }
  return(out)
}

##' @export
##' @rdname mixsurv
dmixsurv = function(dfun, pfun, x, theta, ...) {
  dots <- list(...)
  pargs <- dots
  pargs$lower.tail <- F
  pargs$log.p <- F
  pargs$log <- NULL
  hargs <- dots
  hargs$log <- F
  u_surv <- do.call(pmixsurv, append(list(pfun, x, theta), pargs))
  u_haz <- do.call(hmixsurv, append(list(dfun, pfun, x, theta), hargs))
  out <- u_surv * u_haz
  if (!is.null(dots$log) && dots$log) {
    out <- log(out)
  }
  return(out)
}

##' @export
##' @rdname mixsurv
qmixsurv = function(qfun, p, theta, ...) {
  inv_p <- 1 - p
  dots <- list(...)
  args <- dots
  args$lower.tail <- F
  args$log.p <- F
  if (theta == 0) {
    out <- do.call(qfun, append(list(inv_p), args))
  } else {
    uncured <- inv_p > theta
    out <- rep(Inf, length(inv_p))
    out[uncured] <- do.call(qfun, append(list((inv_p[uncured] - theta) / (1 - theta)), args))
  }
  return(out)
}


##' @export
##' @rdname mixsurv
rmixsurv = function(qfun, n, theta, ...) {

  # Plug random uniform into quantile function
  out <- qmixsurv(qfun, runif(n = n), theta, ...)

  return(out)
}

##' @export
##' @rdname mixsurv
rmst_mixsurv = function(pfun, t, theta, ...) {
  args <- list(...)
  out <- do.call(
    rmst_generic,
    append(
      list(
        function(q, ...) pmixsurv(pfun, q, ...),
        t = t,
        theta = theta
      ),
      args
    )
  )
  return(out)
}

##' @export
##' @rdname mixsurv
mean_mixsurv = function(pfun, theta, ...) {

  # Put together arguments for call to rmst_generic
  args <- append(
    list(
      pfun,
      t = Inf,
      start = 0
    ),
    list(...)
  )

  # Figure out what length the output should be and create
  # a vector to store result
  out_length <- get_param_length_and_check(theta, args)
  out <- numeric(length(out_length))


  # Identify indices where mean survival will be infinite
  # cure fraction is > 0.
  inf_indices <- (theta > 0) & rep(T, out_length)
  out[inf_indices] <- Inf

  # Create arguments to call rmst_generic to estimate mean
  # for indices where cure fraction is zero.
  non_inf_args <- lapply(args, function(x) {

    # Handle x is length 1 and shouldn't be indexed on
    if (length(x) == 1) {
      return(x)
    }

    return(x[!inf_indices])
  })

  # Set output for indices where theta is zero
  non_inf_res <- do.call(rmst_generic, non_inf_args)

  out[!inf_indices] <- non_inf_res

  out
}

