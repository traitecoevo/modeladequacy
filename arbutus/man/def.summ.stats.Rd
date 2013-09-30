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
  function: \enumerate{ \item{reml.sigsq: }{The mean of the
  squared contrasts. This is equivalent to the REML
  estimate of sigsq.}

  \item{var.con: }{The variance of the absolute value of
  the contrasts.}

  \item{slope.con.var: }{The slope of a linear model fit
  between the contrasts and their expected variances.}

  \item{slope.con.asr: }{The slope of a linear model fit
  between the contrasts and their inferred ancestral
  state.}

  \item{slope.con.nh: }{The slope of a linear model fit
  between the contrasts and the node height at which they
  were calculated.}

  \item{ks.pic: }{The D-statistic from a KS test comparing
  the distribution of the contrasts to a normal
  distribution with mean 0 and variance equal to the square
  root of the squared mean of the contrasts.} }
}
\examples{
## get default summary stats
stats <- def.summ.stats
stats
}
\author{
  Matt Pennell
}
\seealso{
  \code{\link{summ.stats}}, \code{\link{sigsq.reml}},
  \code{\link{var.pic}}, \code{\link{slope.con.bl}},
  \code{\link{slope.con.nh}}, \code{\link{slope.con.asr}},
  \code{\link{ks.pic}}
}
\keyword{internal}

