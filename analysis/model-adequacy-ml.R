## Script for evaluating model adequacy of angiosperm clades
## Fitting models using maximum likelihood
source("R/modelfit-helper-fxns.R")

dir.create(path.ml(), FALSE)
options(mc.cores=2) # Adjust accordingly
files <- dir(path.data())
# ok <- mclapply(files, run.model.ad, "ml", mc.preschedule=FALSE)
ok <- lapply(files, run.model.ad, "ml")

# Next, process the output to create a file with summarised results.
fits <- lapply(dir(path.ml(), full.names=TRUE), readRDS)
results <- combine(fits)
write.csv(results, "output/results-ml.csv", row.names=FALSE)
