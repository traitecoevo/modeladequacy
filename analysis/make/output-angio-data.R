#!/usr/bin/env Rscript
source("R/build-angio-data.R")

## Create the directory if it does not exist, but remove all files
## from it if it does:
path <- "output/angio-data"
dir.create(path, FALSE)
invisible(file.remove(dir(path, full.names=TRUE)))

## TODO: Where did these times come from?
tt <- c(0.2687, 50.2697, 100.2697, 150.2697, 200.2697)

## Two nested for loops.  Not the cleanest but it does the trick.
##
## For each taxonomic rank (family, order) and for each time slice
## (elements in tt) we will extract  number of sub trees.  Each of
## these needs saving to its own file.
make.data <- function(trait, min.size=20) {
  dat <- get.data(trait)
  for (rank in c("family", "order")) {
    dat.rank <- build.angio.data.clade(dat, rank, trait, min.size)
    for (i in dat.rank) {
      saveRDS(i, file.path(path, filename.analysis(trait, i$taxa)))
    }
  }

  for (age in tt) {
    dat.time <- build.angio.data.time(dat, age, trait, min.size)
    for (i in seq_along(dat.time)) {
      saveRDS(dat.time[[i]],
              file.path(path, filename.analysis(trait, "timeslice",
                                                age, i)))
    }
  }
}

## Prepare files for each of the three traits.
for (trait in c("SLA", "leafn", "seedmass")) {
  make.data(trait)
}
