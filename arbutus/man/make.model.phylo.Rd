\name{make.model.phylo}
\alias{make.model.phylo}
\title{Rescale phylogeny based on fitted model parameters}
\usage{
  make.model.phylo(x, ...)
}
\arguments{
  \item{x}{an object inherited from a \code{model.info}
  function}

  \item{...}{additional arguments to be passed to
  \code{make.model.phylo}}
}
\value{
  a 'phylo' object that can be used to from a 'unit.tree'
  object
}
\description{
  Rescales phylogeny to form a 'unit.tree' using parameters
  from fitted model
}
\details{
  This is a generic function which rescales the phylogeny
  based on the model specific information. While the class
  and information may differ, the object must include the
  'phylo' object to be rescaled and model specific
  information such as model type and parameter values.

  To include additional types of models, researchers will
  need to build a new \code{\link{model.info}} function for
  the given model type. The output should be assigned a
  novel class. A \code{make.model.phylo.classX} object will
  need to be created to perform the rescaling.
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

## get model info using internal arbutus function
info.bm <- model.info(fit.bm)

## rescale phylogeny based on model info
make.model.phylo(info.bm)
}
}
\seealso{
  \code{\link{as.unit.tree}}
}

