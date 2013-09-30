\name{slope.pic.var}
\alias{slope.pic.var}
\title{Relationship between contrasts and their variances}
\usage{
  slope.pic.var(unit.tree)
}
\arguments{
  \item{unit.tree}{a unit.tree object}
}
\value{
  the estimated slope paramter
}
\description{
  Fits a linear model between the absolute value of
  constrasts and their expected variances
}
\details{
  This summary statistic is the estimate of the slope from
  fitting a linear model between the contrasts and their
  expected variance. It is used to evaluate whether the
  model is adequately capturing variation relative to
  branch lengths. If the generating model is correct, we
  expect there not to be a relationship between the
  contrasts and their variances (i.e. slope ~ 0).

  The slope of the contrasts and their variances is
  included as a default summary statistic in the function
  \code{\link{summ.stats}}
}
\examples{
data(geospiza)
td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
phy <- td$phy
dat <- td$data[,"wingL"]
unit.tree <- as.unit.tree(phy, dat)

## estimate slope
slope.pic.var(unit.tree)
}
\seealso{
  \code{\link{summ.stats}}, \code{\link{def.summ.stats}},
  \code{\link{stats::lm}}
}

