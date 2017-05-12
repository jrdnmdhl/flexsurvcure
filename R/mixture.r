#' @export
pmixsurv = function(pfun, q, theta, ...) {
  dots <- list(...)
  args <- dots
  args$lower.tail <- F
  args$log.p <- F
  out <- theta + (1 - theta) * do.call(pfun, append(list(q), args))
  if (is.null(dots$lower.tail) || dots$lower.tail) {
    out <- 1 - out
  }
  if (!is.null(dots$log.p) && dots$log.p) {
    out <- log(out)
  }
  return(out)
}

#' @export
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

#' @export
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

#' @export
qmixsurv = function(pfun, p, theta, ...) {
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
        function(...) pmixsurv(pfun, ...),
        p = p,
        theta = theta
      ),
      args
    )
  )
  return(out)
}

rmst_mixsurv = function(pfun, t, theta, ...) {
  args <- list(...)
  out <- do.call(
    rmst_generic,
    append(
      list(
        function(...) pmixsurv(pfun, ...),
        t = t,
        theta = theta
      ),
      args
    )
  )
  return(out)
}

mean_mixsurv = function(pfun, t, theta, ...) {
  args <- list(...)
  args$start <- 0
  out <- do.call(
    rmst_generic,
    append(
      list(
        function(...) pmixsurv(pfun, ...),
        t = t,
        theta = theta
      ),
      args
    )
  )
  return(out)
}

