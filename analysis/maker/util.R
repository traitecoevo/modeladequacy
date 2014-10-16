category_to_logical <- function(x, trueval) {
  x[x == ""] <- NA
  x == trueval
}

strstrip <- function(x) {
  gsub("(^ +| +$)", "", x)
}

logspace_f <- function(x, f, ...) {
  exp(f(log(x), ...))
}

geometric_mean <- function(x, ...) {
  logspace_f(x, mean, ...)
}

geometric_sd <- function(x, ...) {
  logspace_f(x, sd, ...)
}

## Quiet version of treedata:
treedata_q <- function(...) {
  geiger::treedata(..., warnings=FALSE)
}

unlist_with_names <- function(x) {
  ret <- unlist(x, recursive=FALSE, use.names=FALSE)
  names(ret) <- unlist(lapply(x, names))
  ret
}
