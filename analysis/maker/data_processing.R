## Converts a raw data set into one with the tree associated with it.
build_data <- function(dat_spp, tree) {
  ## Note that these are *not* the log of one another for the case
  ## where there is >1 observation per species: the raw value will
  ## tend to be larger than the log value because of the difference
  ## between an arithmetic and geometric mean.
  states_log <- structure(dat_spp$mean_log, names=dat_spp$gs)
  states_raw <- structure(dat_spp$mean_raw, names=dat_spp$gs)

  ## From the paper:
  ## > we estimated a single SE for each trait by calculating the mean
  ## > standard deviation for all species for which we had multiple
  ## > measurements
  ## NOTE: Do this *before* filtering by presence in the tree.
  ok <- dat_spp$n_obs > 1L
  se_raw <- mean(dat_spp$sd_raw[ok])
  se_log <- mean(dat_spp$sd_log[ok])

  ## Drop species from tree not in data, and v.v.
  t <- extract.clade(tree, node="Angiospermae")
  phy <- geiger:::.drop.tip(phy=t, tip=setdiff(t$tip.label, dat_spp$gs))
  states_raw <- states_raw[phy$tip.label]
  states_log <- states_log[phy$tip.label]

  list(trait=attr(dat_spp, "trait"), phy=phy,
       states_raw=states_raw, states_log=states_log,
       se_raw=se_raw, se_log=se_log)
}
