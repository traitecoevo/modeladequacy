## Creates the split data sets
make_data_subsets <- function(data) {
  ret <- c(make_data_times(data),
           make_data_families(data),
           make_data_orders(data))
  names(ret) <- NULL
  ret
}

## Converts a raw data set into one with the tree associated with it.
build_data <- function(raw, tree) {
  ## NOTE: base 10 log
  dat <- structure(log10(raw$mean), names=raw$gs)

  ## Drop species from tree not in data, and v.v.
  ## This step is the slow point.
  t <- extract.clade(tree, node="Angiospermae")
  phy <- geiger:::.drop.tip(phy=t, tip=setdiff(t$tip.label, names(dat)))
  dat <- dat[phy$tip.label]
  trait <- attr(raw, "trait")
  se <- attr(raw, "sd")

  list(phy=phy, states=dat, se=se, trait=trait)
}

## Used by both make_data_times and make_data_rank
add_metadata <- function(d, rank, taxa, trait, se) {
  d$rank <- rank
  d$taxa <- taxa
  d$trait <- trait
  d$se <- se
  d
}
