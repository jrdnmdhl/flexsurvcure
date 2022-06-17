get_param_length_and_check <- function(theta, args) {
  param_lengths <- c(length(theta), sapply(args, length))
  max_length <- max(param_lengths)
  valid_lengths <- (param_lengths <= 1) | (param_lengths == max_length)

  # If any two parameters with length > 1 have different lengths
  # then throw an error.
  if (!all(valid_lengths)) {
    stop('Parameter values provided were of incompatible length')
  }

  max_length
}
