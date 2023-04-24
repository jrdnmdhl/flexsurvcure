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
##' @param qfun The base distribution's quantile function.
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
##' @rdname nmixsurv
hnmixsurv = function(dfun,x, theta, ...) {
  dots <- list(...)
  args <- dots
  args$log <- F
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
  surv <- do.call(pnmixsurv, append(list(pfun, x, theta), pargs))
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
qnmixsurv = function(qfun, p, theta, ...) {
  inv_p <- 1 - p
  dots <- list(...)
  args <- dots
  args$lower.tail <- F
  args$log.p <- F
  uncured <- inv_p > theta
  out <- rep(Inf, length(inv_p))
  if (length(theta)==1) {
    theta <- rep(theta, length(out))
  }

  zeroThetaInd = uncured & theta == 0

  if (any(zeroThetaInd)) {
    # If no cure then just use base qfun
    out[zeroThetaInd] <- do.call(qfun, append(list(inv_p[zeroThetaInd]), args))
  } else {
    # Calculations below are meant to map the quantile distribution of the base
    # distribution to the quantile distribution of the cure model using the
    # following algebra:
    #
    # S(t) = theta ^ (1 - Su(t))
    # ln[S(t)] = ln[theta ^ (1 - Su(t))]
    # ln[S(t)] = (1 - Su(t)) * ln[theta]
    # ln[S(t)] / ln[theta] = 1 - Su(t)
    # Su(t) = 1 - ln[S(t)] / ln[theta]
    #
    # Where Su(t) is the baseline survival distribution
    p_surv_to_lookup <- 1 - (log(inv_p) / log(theta))
    out[uncured] <- do.call(qfun, append(list(p_surv_to_lookup), args))[uncured]
  }
  return(out)
}


##' @export
##' @rdname nmixsurv
rnmixsurv = function(qfun, n, theta, ...) {

  # Plug random uniform into quantile function
  out <- qnmixsurv(qfun, runif(n = n), theta, ...)

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
        function(q, ...) pnmixsurv(pfun, q, ...),
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

  # This is a very silly function because if theta is greater
  # than zero then the mean is infinite and if theta is zero
  # then the mean is zero. Still need to have it since all
  # flexsurv models should support getting means through
  # summary.flexsurv.

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

  out[!inf_indices] <- 0

  out

}

