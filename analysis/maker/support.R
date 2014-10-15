build_data <- function(raw, tree) {
  ## NOTE: base 10 log
  dat <- structure(log10(raw$mean), names=raw$gs)

  ## Drop species from tree not in data, and v.v.
  ## This step is the slow point.
  t <- extract.clade(tree, node="Angiospermae")
  phy <- geiger:::.drop.tip(phy=t, tip=setdiff(t$tip.label, names(dat)))
  dat <- dat[phy$tip.label]
  se <- attr(raw, "se")

  list(phy=phy, states=dat, SE=se)
}
