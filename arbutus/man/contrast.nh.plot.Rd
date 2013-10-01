\name{contrast.nh.plot}
\alias{contrast.nh.plot}
\title{Plot contrasts versus their node height}
\usage{
  contrast.nh.plot(x, col = c("dodgerblue4", "darkblue"),
    ...)
}
\arguments{
  \item{x}{a model fitted object or 'phylo' object that can
  be passed to \code{\link{as.unit.tree}}.}

  \item{col}{colours used for plotting}

  \item{...}{additional arguments to be passed to
  \code{\link{as.unit.tree}}}
}
\description{
  Plot the contrasts derived from a rescaled 'unit.tree'
  against the node height (i.e. depth) at which they are
  computed
}
\details{
  This function generates a plot of the observed contrasts
  versus their the node height at which they are computed
  (see \code{\link{ape::pic}}). The slope of this plot is
  equal to the summary statistic 'slope.pic.nh', calculated
  by \code{\link{slope.pic.nh}} and is used as a default
  summary statistic in \code{\link{summ.stats}}.

  If the model is a good fit to the data (i.e. an adequate
  model), the slope ~ 0. The purpose of the plot is to get
  a quick visual glance at the data and to assess whether
  or not the model inadequacy may be driven by outliers.

  The function can currently only take a single dataset
  (does not integrate across multiple trees).

  \code{contrast.nh.plot} uses the
  \code{\link{ggplot2::ggplot2}} library.
}
\examples{
## finch data
data(geospiza)
td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
phy <- td$phy
data <- td$data[,"wingL"]

## fit Brownian motion model
## using geiger's fitContinuous function
fit.bm <- fitContinuous(phy=phy, dat=data, model="BM",
                                 control=list(niter=10))

## plot contrasts versus their variances
contrast.nh.plot(fit.bm)
}
\seealso{
  \code{\link{as.unit.tree}}, \code{\link{summ.stats}},
  \code{\link{slope.pic.nh}}
}

