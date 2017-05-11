#' @export
flexsurvcure <- function(formula, data, link = "logistic", mixture = T, ...) {
  dist = flexsurvcure.dists$weibull
  if (link == "logistic") {
    dist$transforms[[1]] <- gtools::logit
    dist$inv.transforms[[1]] <- gtools::inv.logit
  } else if(link == "loglog") {
    dist$transforms[[1]] <- function(x) log(-log(x))
    dist$inv.transforms[[1]] <- function(x) exp(-exp(x))
  } else {
    dist$transforms[[1]] <- identity
    dist$inv.transforms[[1]] <- identity
  }
  print(list(...))
  out <- do.call(
    flexsurv::flexsurvreg,
    list(
      formula,
      data = data,
      dist = dist,
      ...
    )
  )
  class(out) <- c("flexsurvcure", "flexsurv")
  out
}
