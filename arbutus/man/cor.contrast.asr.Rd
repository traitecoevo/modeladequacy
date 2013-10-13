\name{cor.contrast.asr}
\alias{cor.contrast.asr}
\title{Relationship between contrasts and inferred ancestral state}
\usage{
  cor.contrast.asr(unit.tree)
}
\arguments{
  \item{unit.tree}{a 'unit.tree' object}
}
\value{
  the estimated slope parameter
}
\description{
  Fits a linear model between the absolute value of the
  contrasts and the ancestral state inferred at the node at
  which the contrast was calculated
}
\details{
  This summary statistic is the estimate of the slope from
  fitting a linear model between the contrasts and the
  inferred ancestral state at the node at which the
  contrasts were calculated. It is used to evaluate whether
  the model is adequately capturing variation relative to
  ancestral state. The ancestral states are calculated
  using the contrast based approach and are not technically
  ancestral state estimates (see Felsenstein 1985 --
  American Naturalist); this is conceptually and
  statistically different from the maximum likelihood
  approach to reconstructing ancestral states. The contrast
  based approach is used here as there is a clear
  expectation of the relationship between the the contrasts
  and trait value. If the generating model is correct, we
  expect there not to be a relationship between the
  contrasts and their node height (i.e. slope ~ 0).

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
cor.contrast.asr(unit.tree)
}
\seealso{
  \code{\link{summ.stats}}, \code{\link{def.summ.stats}},
  \code{\link[stats]{lm}}, \code{\link[ape]{ace}}
}

