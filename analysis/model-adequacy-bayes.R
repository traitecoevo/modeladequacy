## Script for evaluating model adequacy of angiosperm clades
## Fitting models using mcmc
source("R/model-adequacy-fit.R")

dir.create(path.bayes(), FALSE)

options(mc.cores=2) # Adjust accordingly
files <- dir(path.data())
ok <- mclapply(files, run.model.ad, "bayes", mc.preschedule=FALSE)

# All the fits can be read in like this, but it takes a while and
# results in objects that are large.
fits <- lapply(dir(path.bayes(), full.names=TRUE), readRDS)

# Combine everything together by stripping out some summary statistics
results <- combine_fits(fits)
write.csv(results, "output/results-bayes.csv", row.names=FALSE)
