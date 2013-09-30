\name{phy.model.check}
\alias{phy.model.check}
\title{Assess the adequacy of continuous trait models}
\usage{
  phy.model.check(x, nsim = 1000, stats = NULL, ...)
}
\arguments{
  \item{x}{a fitted model object. See
  \code{\link{as.unit.tree}} for a description of what type
  of model types are implemented}

  \item{nsim}{the number of datasets to simulate. This is
  passed as an argument to the function
  \code{\link{sim.char.unit}}.}

  \item{stats}{a named list of summary statistics to
  calculate on observed and simulated datasets. See
  \code{\link{summ.stats}} for details.}

  \item{...}{additional arguments to be passed to
  \code{\link{as.unit.tree}}}
}
\description{
  Use summary statistics to assess model adequacy for
  phylogenetic models of continuous character evolution.
  This function is a simply a wrapper which does all the
  steps involved in evaluating model adequacy using the
  approach outlined in Pennell et al. (citation)
}
\details{
  lots of details here. should be a mini-tutorial on the
  method and borrow heavily from the manuscript.
}
\examples{
## finch data
data(geospiza)
td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
phy <- td$phy
data <- td$data[,"wingL"]

## fit Brownian motion model
## using geiger's fitContinuous function
fit.bm <- fitContinuous(phy=phy, dat=data, model="BM",
                                 control=list(niter=10))

## check adequacy of BM model
## get p-values for default summary statistics
modelad.bm <- phy.model.check(fit.bm, nsim=10)
modelad.bm


## fit Ornstein-Uhlenbeck model
## again, using geiger's fitContinuous function
fit.ou <- suppressWarnings(fitContinuous(phy=phy, dat=data, model="OU",
                                                  control=list(niter=10)))

## check adequacy of OU model
modelad.ou <- phy.model.check(fit.ou, nsim=10)
}
\seealso{
  \code{\link{as.unit.tree}}, \code{\link{summ.stats}},
  \code{\link{sim.char.unit}},
  \code{\link{compare.summ.stats}}
}

