\name{print.phy.ss}
\alias{print.phy.ss}
\title{Print summary statistics}
\usage{
  \method{print}{phy.ss} (x)
}
\arguments{
  \item{x}{a 'phy.ss' object}
}
\value{
  prints p.values to screen
}
\description{
  Generic function for printing 'phy.ss' objects from
  \code{\link{compare.summ.stats}}
}
\examples{
data(geospiza)
td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
phy <- td$phy
dat <- td$data[,"wingL"]

## fit BM model using fitContinuous
f <- fitContinuous(phy, dat, model="BM", control=list(niter=10))

## use phy.model.check to assess model adequacy
p <- phy.model.check(f, nsim=10)

## print output
p
}

