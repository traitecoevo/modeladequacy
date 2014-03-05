#!/usr/bin/env Rscript
tmp <- tempdir()
unzip(file.path('data/zae/PhylogeneticResources.zip'),
      'PhylogeneticResources/Vascular_Plants_rooted.dated.tre',
      junkpaths=TRUE, exdir=tmp)
file.rename(file.path(tmp, "Vascular_Plants_rooted.dated.tre"),
            "data/vascular_plant_phylogeny.tre")
invisible(NULL)
