#!/usr/bin/env Rscript
prefix <- 'http://datadryad.org/bitstream/handle/10255/'
suffix <- 'dryad.55548/PhylogeneticResources.zip?sequence=1'
dest <- "data/zae"
dir.create(dest, FALSE)
download.file(paste0(prefix, suffix),
              file.path(dest, "PhylogeneticResources.zip"))
