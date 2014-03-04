#!/usr/bin/env Rscript
prefix <- 'http://datadryad.org/bitstream/handle/10255/'
suffix <- 'dryad.55304/Spermatophyta_Genera.csv?sequence=2'
dest <- "data/zae"
dir.create(dest, FALSE)
download.file(paste0(prefix, suffix),
              file.path(dest, "Spermatophyta_Genera.csv"))
