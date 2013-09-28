\name{var.pic}
\alias{var.pic}
\title{Variance of contrasts}
\usage{
  var.pic(unit.tree)
}
\arguments{
  \item{unit.tree}{a 'unit.tree' object}
}
\value{
  the estimate of the variance of the absolute value of the
  contrasts
}
\description{
  Calculates the variance of the absolute value of the
  contrasts
}
\details{
  This summary is used to evaluate whether the model is
  sufficiently capturing the variation in rates across the
  tree.

  The variance of contrasts is included as a default
  summary statistic in the function
  \code{\link{summ.stats}}
}
\examples{
data(geospiza)
td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
phy <- td$phy
dat <- td$data[,"wingL"]
unit.tree <- as.unit.tree(phy, dat)

## estimate variance of contrasts
var.pic(unit.tree)
}
\author{
  Matt Pennell
}
\seealso{
  \code{\link{summ.stats}}, \code{\link{def.summ.stats}}
}

