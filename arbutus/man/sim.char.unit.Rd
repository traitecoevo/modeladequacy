\name{sim.char.unit}
\alias{sim.char.unit}
\title{Simulate character evolution on unit tree}
\usage{
  sim.char.unit(unit.tree, nsim = 1000)
}
\arguments{
  \item{unit.tree}{a 'unit.tree' object}

  \item{nsim}{the number of datasets to simulate (if single
  unit.tree is provided; see below)}
}
\value{
  A list of unit.trees. In each unit.tree, the phylogeny
  will be the same as that of the input unit.tree(s) but
  the data and the contrasts will be replaced by the
  simulated datasets.
}
\description{
  Simulates continuous characters along a unit tree
  according to a Brownian motion process with rate 1.
}
\details{
  If a single unit.tree is supplied, the function will
  simulate \code{nsim} datasets under a Brownian motion
  process with a rate of 1. If a lits of unit.trees
  supplied (such as those derived from a Bayesian analysis)
  \code{sim.char.unit} will simulate a single data set on
  each tree in the list and the \code{nsim} argument will
  be ignored.
}
\examples{
data(geospiza)
td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
phy <- td$phy
dat <- td$data[,"wingL"]

## build unit.tree object
unit.tree <- as.unit.tree(phy, dat)

## simulate 2 datasets on unit tree
sims <- sim.char.unit(unit.tree, nsim=2)

sims
}
\author{
  Matt Pennell
}
\seealso{
  \code{\link{geiger::sim.char}} which this function wraps
}

