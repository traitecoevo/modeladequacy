\name{pval.summ.stats}
\alias{pval.summ.stats}
\title{Extract p--values for summary statistics}
\usage{
  pval.summ.stats(x)
}
\arguments{
  \item{x}{a 'phy.ss' object from the function
  \code{\link{compare.summ.stats}}}
}
\value{
  a named vector of two-tailed p-values
}
\description{
  Utility function for extracting p-values from the output
  of \code{\link{compare.summ.stats}}
}
\examples{
data(finch)
phy <- finch$phy
dat <- finch$data[,"wingL"]
unit.tree <- as.unit.tree(phy, data=dat)

## calculate default summary stats on observed data
ss.obs <- summ.stats(unit.tree, stats=NULL)

## simulate data on unit.tree
sims <- sim.char.unit(unit.tree, nsim=10)

## calculate default summary stats on simulated data
ss.sim <- summ.stats(sims, stats=NULL)

## compare simulated to observed summary statistics
res <- compare.summ.stats(ss.obs, ss.sim)

## extract p-values
pval.summ.stats(res)
}
\seealso{
  \code{\link{compare.summ.stats}}
}

