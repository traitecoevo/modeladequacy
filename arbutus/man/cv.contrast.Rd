\name{cv.contrast}
\alias{cv.contrast}
\title{Coefficient of variation of contrasts}
\usage{
  cv.contrast(unit.tree)
}
\arguments{
  \item{unit.tree}{a \code{unit.tree} object}
}
\value{
  \code{v.pic} the coefficient of variation of the
  contrasts
}
\description{
  Calculates the coefficient of variation on a set of
  contrasts
}
\details{
  This summary is used to evaluate whether the model is
  sufficiently capturing the variation in rates across the
  tree.

  The coefficient of variation is the sd/mean. The
  coefficient of variation of the contrasts is included as
  a default summary statistic in the function
  \code{\link{summ.stats}}
}
\examples{
data(finch)
phy <- finch$phy
dat <- finch$data[,"wingL"]
unit.tree <- as.unit.tree(phy, data=dat)

## estimate variance of contrasts
cv.contrast(unit.tree)
}
\seealso{
  \code{\link{summ.stats}}, \code{\link{def.summ.stats}}
}

