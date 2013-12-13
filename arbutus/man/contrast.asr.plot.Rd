\name{contrast.asr.plot}
\alias{contrast.asr.plot}
\title{Plot contrasts versus their inferred ancestral states}
\usage{
  contrast.asr.plot(x, col = c("dodgerblue4", "darkblue"),
    ...)
}
\arguments{
  \item{x}{a model fitted object or \code{phylo} object
  that can be passed to \code{\link{as.unit.tree}}.}

  \item{col}{colours used for plotting}

  \item{...}{additional arguments to be passed to
  \code{\link{as.unit.tree}}}
}
\description{
  Plot the contrasts derived from a rescaled 'unit.tree'
  against their inferred ancestral state.
}
\details{
  This function generates a plot of the observed contrasts
  versus their inferred ancestral state (see
  \code{\link[ape]{pic}}, \code{\link[ape]{ace}}). Note
  that the ancestral states are estimated using the 'pic'
  method, which differs somewhat from the ML estimate of
  the ancestral states. The slope of this plot is equal to
  the summary statistic \code{\link{cor.contrast.asr}} and
  is used as a default summary statistic in
  \code{\link{summ.stats}}.

  If the model is a good fit to the data (i.e. an adequate
  model), the slope ~ 0. The purpose of the plot is to get
  a quick visual glance at the data and to assess whether
  or not the model inadequacy may be driven by outliers.

  The function can currently only take a single dataset
  (does not integrate across multiple trees).

  \code{contrast.asr.plot} uses the
  \code{\link[ggplot2]{ggplot2}} library.
}
\examples{
## finch data
data(finch)
phy <- finch$phy
data <- finch$data[,"wingL"]

\dontrun{
require(geiger)
## fit Brownian motion model
## using geiger's fitContinuous function
fit.bm <- fitContinuous(phy=phy, dat=data, model="BM",
                                 control=list(niter=10))

## plot contrasts versus their variances
contrast.asr.plot(fit.bm)
}
}
\seealso{
  \code{\link{as.unit.tree}}, \code{\link{summ.stats}},
  \code{\link{cor.contrast.asr}}
}

