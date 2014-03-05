#!/usr/bin/env Rscript
source("R/build-angio-data.R")

## TODO: Be more specific perhaps.
## read in and process data

## TODO: Where did these come from?
## time slices
tt <- c(0.2687, 50.2697, 100.2697, 150.2697, 200.2697)

sla <- get.sla.data()
sla.fam <- build.angio.data.clade(sla, rank="family", trait="SLA", min.size=20)
sla.ord <- build.angio.data.clade(sla, rank="order", trait="SLA", min.size=20)
sla.tt <- lapply(tt, function(t)
                 build.angio.data.time(sla, age=t, trait="SLA", min.size=20))

sdm <- get.seedmass.data()
sdm.fam <- build.angio.data.clade(sdm, rank="family", trait="seedmass", min.size=20)
sdm.ord <- build.angio.data.clade(sdm, rank="order", trait="seedmass", min.size=20)
sdm.tt <- lapply(tt, function(t)
                 build.angio.data.time(sdm, age=t, trait="seedmass", min.size=20))

lfn <- get.leafn.data()
lfn.fam <- build.angio.data.clade(lfn, rank="family", trait="leafn", min.size=20)
lfn.ord <- build.angio.data.clade(lfn, rank="order", trait="leafn", min.size=20)
lfn.tt <- lapply(tt, function(t)
                 build.angio.data.time(lfn, age=t, trait="leafn", min.size=20))

all.trait.data <- c(sla.fam, sla.ord, unlist(sla.tt, FALSE),
                    sdm.fam, sdm.ord, unlist(sdm.tt, FALSE),
                    lfn.fam, lfn.ord, unlist(lfn.tt, FALSE))

## write to rds
saveRDS(all.trait.data, "output/angio-trait-data-all.rds")
