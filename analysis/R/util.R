
logspace.f <- function(x, f, ...) {
  exp(f(log(x), ...))
}

geometric.mean <- function(x, ...) {
  logspace.f(x, mean, ...)
}

geometric.sd <- function(x, ...) {
  logspace.f(x, sd, ...)
}
