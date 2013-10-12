\name{print.phy.ss}
\alias{print.phy.ss}
\title{Print summary statistics}
\usage{
  \method{print}{phy.ss} (x, ...)
}
\arguments{
  \item{x}{a 'phy.ss' object}

  \item{...}{additional arguments to be passed to print}
}
\value{
  prints p.values to screen
}
\description{
  Generic function for printing 'phy.ss' objects from
  \code{\link{compare.summ.stats}}
}
\examples{
data(finch)
phy <- finch$phy
dat <- finch$data[,"wingL"]

## use phy.model.check to assess model adequacy
## pretend as if tree has already been rescaled
p <- phy.model.check(phy, data=dat, nsim=10)

## print output
p
}

