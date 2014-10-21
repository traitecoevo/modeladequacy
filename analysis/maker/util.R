category_to_logical <- function(x, trueval) {
  x[x == ""] <- NA
  x == trueval
}

strstrip <- function(x) {
  gsub("(^ +| +$)", "", x)
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
