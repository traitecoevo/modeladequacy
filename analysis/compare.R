d.new <- read.csv("output/results-ml.csv", stringsAsFactors=FALSE)
d.old <- read.csv("output/ml-results.csv", stringsAsFactors=FALSE,
                  row.names=1)
rownames(d.old) <- NULL

# Check that all columns are present:
all(names(d.new) %in% names(d.old))
all(names(d.old) %in% names(d.new))

# Reorder rows because we go through the analyses in a slightly
# different order:
d.new <- d.new[names(d.old)]

dd.new <-
  d.new[order(d.new$trait, d.new$rank, d.new$rank, d.new$size, d.new$age),]
dd.old <-
  d.old[order(d.old$trait, d.old$rank, d.old$rank, d.old$size, d.old$age),]
rownames(dd.new) <- rownames(dd.old) <- NULL

info <- c("taxa", "rank", "trait", "size", "age")
all.equal(dd.new[info], dd.old[info])

## Looks like we very slightly disagree on sizes, but I'm not freaking
## out about this.  Could be treatment of rounding issues.
plot(dd.new$size ~ dd.old$size, log="xy")
plot(dd.new$size - dd.old$size ~ dd.old$size)

## This is more concerning:
plot(dd.new$age ~ dd.old$age, log="xy")
plot(dd.new$age - dd.old$age ~ dd.old$size)

ok <- abs(dd.new$age - dd.old$age) < 1e-4 & dd.new$size == dd.old$size

# Focussing just on the cases where we have the same basic
# information:
plot(dd.new$aic.bm[ok] - dd.old$aic.bm[ok])
plot(dd.new$aic.ou[ok] - dd.old$aic.ou[ok])
plot(dd.new$aic.eb[ok] - dd.old$aic.eb[ok])

## My guess is that the BM and EB cases are fine (all the differences
## are probably due to a change in the kew data).  But some of the OU
## cases look to have failed to converge.

## These look quite different to me: worryingl so, quite frankly.
plot(dd.new$mv.ml.bm[ok] - dd.old$mv.ml.bm[ok])

i <- which(ok)[order(abs(dd.new$mv.ml.bm[ok] - dd.old$mv.ml.bm[ok]))]
dd.new[i,1:5]

## Script for evaluating model adequacy of angiosperm clades
## Fitting models using maximum likelihood
source("R/modelfit-helper-fxns.R")

## Look at the error around the distance estimate:
dat <- readRDS(file.path(path.data(), "SLA_clade_Poales.rds"))
res <- lapply(1:10, function(i) model.ad(dat, "BM", "ml", i))
dist <- sapply(res, function(x) mean(mv.summ.stats(x$ma)))
# About 5% CI
sd(dist) / mean(dist)
