## Script for evaluating model adequacy of angiosperm clades
## Fitting models using maximum likelihood
source("R/model-adequacy-fit.R")

dir.create(path.ml(), FALSE)

options(mc.cores=2) # Adjust accordingly
files <- dir(path.data())
ok <- mclapply(files, run.model.ad, "ml", mc.preschedule=FALSE)

# All the fits can be read in like this, but it takes a while and
# results in objects that are large.
fits <- lapply(dir(path.ml(), full.names=TRUE), readRDS)

# Combine everything together by stripping out some summary statistics
results <- combine_fits(fits)
write.csv(results, "output/results-ml.csv", row.names=FALSE)
