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

capitalise_first <- function(x) {
  sub("^([a-z])", "\\U\\1", x, perl=TRUE)
}

rename_traits <- function(x, as_factor=TRUE) {
  tr <- c(sla="SLA",
          seed_mass="SeedMass",
          leaf_n="LeafN")
  rename(x, tr, as_factor)
}

rename_variables <- function(x, as_factor=TRUE) {
  from <- pvalue_names_arbutus()
  tr <- sub("([a-z])\\.([a-z]+)", "italic(\\U\\1[\\U\\2])", from,
            perl=TRUE)
  names(tr) <- from
  rename(x, tr, as_factor)
}

rename <- function(x, tr, as_factor=TRUE) {
  ret <- unname(tr[match(as.character(x), names(tr))])
  if (as_factor) {
    ret <- factor(ret, unname(tr))
  }
  ret
}

write_csv <- function(x, ...) {
  write.csv(x, ..., row.names=FALSE)
}
