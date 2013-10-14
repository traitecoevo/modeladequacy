\name{cor.contrast.nh}
\alias{cor.contrast.nh}
\title{Relationship between contrasts and node height}
\usage{
  cor.contrast.nh(unit.tree)
}
\arguments{
  \item{unit.tree}{a 'unit.tree' object}
}
\value{
  the estimated slope parameter
}
\description{
  Fits a linear model between the absolute value of the
  contrasts and the node height at which they were
  calculated
}
\details{
  This summary statistic is the estimate of the slope from
  fitting a linear model between the contrasts and the node
  height (i.e. tree depth) at which they were calculated.
  It is used to evaluate whether the model is adequately
  capturing variation relative to tree depth. If the
  generating model is correct, we expect there not to be a
  relationship between the contrasts and their node height
  (i.e. slope ~ 0).

  The slope of the contrasts and their variances is
  included as a default summary statistic in the function
  \code{\link{summ.stats}}
}
\examples{
data(finch)
phy <- finch$phy
dat <- finch$data[,"wingL"]
unit.tree <- as.unit.tree(phy, data=dat)

## estimate slope
cor.contrast.nh(unit.tree)
}
\seealso{
  \code{\link{summ.stats}}, \code{\link{def.summ.stats}},
  \code{\link[stats]{lm}}
}

