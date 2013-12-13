\name{as.unit.tree}
\alias{as.unit.tree}
\title{Create unit.tree object}
\usage{
  as.unit.tree(x, ...)
}
\arguments{
  \item{x}{fitted model object or
  \code{phylo}/\code{multiPhylo} object}

  \item{...}{additional arguments}
}
\value{
  a \code{unit.tree} object containing (or a list of
  \code{unit.tree} objects, each containing) the following
  elements: \describe{ \item{phy}{a \code{phylo} object. If
  a model fitted object has been supplied (see details),
  the \code{phylo} object will be rescaled based on fitted
  model parameters.} \item{data}{the original comparative
  data} \item{pics}{a matrix consisting of the contrasts
  calculated using the original data on the (rescaled)
  phylogeny. The matrix has two columns: the "contrasts"
  and a "variance" for each contrast.} }
}
\description{
  Generic method for creating a 'unit.tree' object
}
\details{
  This function is a generic function which takes a fitted
  model object (from a number of different packages),
  rescales the phylogeny based on the model fitted
  parameters, computes the contrasts on the rescaled
  phylogeny and returns an object of class
  \code{unit.tree}. The \code{unit.tree} object can then be
  used to assess the adequacy of phylogenetic models of
  continuous character evolution. Alternatively, the
  function can take a \code{phylo} object or a
  \code{multiPhylo} object and compute the contrasts on the
  phylogeny without rescaling. This is only meaningful for
  downstream analyses if the phylogenies have been rescaled
  beforehand. Currently, the following object types have
  been implemented: \itemize{ \item a \code{gfit} object
  returned from fitting a model of continuous character
  evolution using \code{fitContinuous} in the \code{geiger}
  package.

  \item a \code{fit.mle} object returned from fitting a
  model of continuous character evolution using
  \code{find.mle} in the \code{diversitree} package. As the
  \code{fit.mle} object does not include all of the
  information required for creating a \code{unit.tree}, a
  second argument \code{lik} needs to be supplied,
  providing the likelihood function used in
  \code{find.mle}.

  \item a \code{mcmcsamples} object returned from fitting a
  model of continuous character evolution using MCMC
  methods in the \code{diversitree} package.
  \code{as.unit.tree} will apply the same trait dataset to
  a set of unit trees based on sampled parameters. By
  default this will create a unit tree for every sample in
  the mcmc chain. To modify this, additional arguments can
  be use.  \code{burnin} specifies how many samples to
  remove from the beginning of the chain.  \code{thin}
  specifies the thinning interval (e.g. if \code{thin=5},
  the function will create a unit tree from every fifth
  parameter set sampled.  \code{sample} specifies how many
  samples to draw from the MCMC run.  As the
  \code{mcmcsamples} object does not include all of the
  information required for creating a \code{unit.tree}, a
  second argument \code{lik} needs to be supplied,
  providing the likelihood function used in
  \code{mcmcsamples}.

  \item a \code{gls} object returned from fitting a
  phylogenetic generalized least squares model of character
  correlation using \code{gls} in the \code{nlme} package.

  \item a \code{pgls} object returned from fitting a
  phylogenetic generalized least squares model of character
  correlation using \code{pgls} in the \code{caper}
  package.

  \item a \code{phylolm} object returned from fitting a
  phylogenetic generalized linear model of character
  correlation using \code{phylolm} in the \code{phylolm}
  package. As the phylogeny is not returned with the
  \code{phylolm} object, a \code{phy} argument must also be
  specified.

  \item a \code{phylo} object. If a \code{phylo} object is
  supplied, the tree is assumed to have been rescaled
  previously. A second argument \code{data} must also be
  provided included the trait data as a named vector with
  names equal to the tip.labels of the phylogeny.

  \item a \code{multiPhylo object}. If a \code{multiPhylo}
  object is supplied, the tree is assumed to have been
  rescaled previously. A second argument \code{data} must
  also be provided included the trait data as a named
  vector with names equal to the tip.labels of the
  phylogenies. Note that this function will append the same
  data set to every tree in the \code{multiPhylo} object. }
}
\examples{
## finch data
data(finch)
phy <- finch$phy
dat <- finch$data[,"wingL"]

## using just the given phylogeny
unit.tree.phy <- as.unit.tree(phy, data=dat)

\dontrun{
require(geiger)
## fit Brownian motion model
## using geiger's fitContinuous function
fit.bm <- fitContinuous(phy=phy, dat=data, model="BM",
                                 control=list(niter=10))

## this creates a 'gfit' object which can be used
## in 'as.unit.tree()'
unit.tree.geiger <- as.unit.tree(fit.bm)
unit.tree.geiger

require(diversitree)
## fit Brownian motion model using ML
## using diversitree's find.mle function

bmlik <- make.bm(phy, data)
fit.bm.dt <- find.mle(bmlik, 0.1)

## this creates a 'fit.mle' object which can be used
## in 'as.unit.tree()'
unit.tree.dt <- as.unit.tree(fit.bm.dt, lik=bmlik)

## fit Brownian motion model using MCMC
mcmc.bm.dt <- mcmc(bmlik, x.init=1, nsteps=1000, w=1)

## construct a unit tree object from mcmcsamples
## removing 100 samples as burnin
## and sampling 10 parameter sets
unit.tree.mcmc <- as.unit.tree(mcmc.bm.dt, lik=bmlik,
                       burnin=100, samples=10)


require(nlme)
## Use pgls to look for a correlation between two traits

t1 <- data
t2 <- td$data[,"tarsusL"]
dd <- cbind.data.frame(t1, t2)

## fit gls model with corPagel correlation structure
fit.gls <- gls(t1~t2, data=dd, correlation=corPagel(phy=phy, value=1))

## this creates a 'gls' object which can be used
## in 'as.unit.tree()'
unit.tree.gls <- as.unit.tree(fit.gls)

}
}
\seealso{
  \code{\link[ape]{pic}}, \code{\link{phy.model.check}}
}

