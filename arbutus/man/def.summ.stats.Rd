\name{def.summ.stats}
\alias{def.summ.stats}
\title{Internal function for getting default summary statistics}
\usage{
  def.summ.stats()
}
\value{
  a named list of functions
}
\description{
  Creates a list of the default summary statistics to be
  used to asses model adequacy
}
\details{
  The following summary statistics are produced by this
  function: \enumerate{ \item{m.pic: }{The mean of the
  squared contrasts. This is equivalent to the REML
  estimate of sigsq.}

  \item{v.pic: }{The coefficient of variation of the
  absolute value of the contrasts.}

  \item{s.var: }{The slope of a linear model fit between
  the contrasts and their expected variances.}

  \item{s.anc: }{The slope of a linear model fit between
  the contrasts and their inferred ancestral state.}

  \item{s.hgt: }{The slope of a linear model fit between
  the contrasts and the node height at which they were
  calculated.}

  \item{d.ks: }{The D-statistic from a KS test comparing
  the distribution of the contrasts to a normal
  distribution with mean 0 and variance equal to the square
  root of the squared mean of the contrasts.} }
}
\examples{
## get default summary stats
stats <- def.summ.stats
stats
}
\seealso{
  \code{\link{summ.stats}}, \code{\link{sigsq.est}},
  \code{\link{cv.contrast}},
  \code{\link{cor.contrast.var}},
  \code{\link{cor.contrast.nh}},
  \code{\link{cor.contrast.nh}}, \code{\link{ks.contrast}}
}
\keyword{internal}

