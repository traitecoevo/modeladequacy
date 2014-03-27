## Script for evaluating model adequacy of angiosperm clades
## Fitting models using mcmc

## Read in all the helper files
source("R/modelfit-helper-fxns.R")

dir.create(path.bayes(), FALSE)

options(mc.cores=2) # Adjust accordingly
files <- dir(path.data())
# ok <- mclapply(files, run.model.ad, "bayes", mc.preschedule=FALSE)
ok <- lapply(files, run.model.ad, "bayes")


tmp <- run.model.ad(files[[1]], "bayes")

# Next, process the output to create a file with summarised results.
fits <- lapply(dir(path.bayes(), full.names=TRUE), readRDS)
results <- combine(fits)
write.csv(results, "output/results-bayes.csv", row.names=FALSE)
