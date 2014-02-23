\name{cor.contrast.var}
\alias{cor.contrast.var}
\title{Relationship between contrasts and their variances}
\usage{
  cor.contrast.var(unit.tree)
}
\arguments{
  \item{unit.tree}{a \code{unit.tree} object}
}
\value{
  \code{s.var} the estimated slope paramter
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
  branch lengths. If the model is adequate, we expect there
  not to be a relationship between the contrasts and their
  variances (i.e. slope ~ 0).

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
cor.contrast.var(unit.tree)
}
\seealso{
  \code{\link{summ.stats}}, \code{\link{def.summ.stats}},
  \code{\link[stats]{lm}}
}

