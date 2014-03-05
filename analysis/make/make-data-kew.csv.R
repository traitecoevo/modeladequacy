#!/usr/bin/env Rscript
source("R/kew.R")

clades <- kew.clades.tr(kew.clades())

for (clade in clades) {
  message(sprintf("Processing: %s", clade))
  kew.html.to.csv(sprintf("data/kew/kew_apg_%s.html", clade))
}

files <- sprintf("data/kew/kew_apg_%s.csv", clades)
dat <- lapply(files, read.csv, stringsAsFactors=FALSE)
res <- cbind(clade=I(rep(clades, sapply(dat, nrow))),
             do.call(rbind, dat))
write.csv(res, "data/kew.csv", row.names=FALSE)
