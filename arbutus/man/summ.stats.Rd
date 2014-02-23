\name{summ.stats}
\alias{summ.stats}
\title{Calculate summary statistics on a unit tree}
\usage{
  summ.stats(unit.tree, stats = NULL)
}
\arguments{
  \item{unit.tree}{a \code{unit.tree} object or list of
  unit.trees}

  \item{stats}{a named list of summary statistics to
  calculate on the unit.tree.  If no \code{stats} argument
  supplied, default summary statistics are used.}
}
\value{
  A data.frame with the calculated summary statistics
  across all unit.trees provided.
}
\description{
  Calculates a set of summary statistics on the contrasts
  included with a 'unit.tree' object.
}
\details{
  This function can be applied to either a single
  'unit.tree' of object or a list of 'unit.tree' objects.
  If \code{stats=NULL} default summary statistics are used.
  The default summary statistics are the following:
  \enumerate{ \item{m.pic: }{The mean of the squared
  contrasts. This is equivalent to the REML estimate of
  sigsq.}

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
  root of the squared mean of the contrasts.} } User
  defined summary statistics can be supplied as a named
  list of functions (see examples). The functions supplied
  must take a unit.tree as argument and perform an
  operation on at least one of the elements of the object
  (see \code{\link{as.unit.tree}} for details).
}
\examples{
data(finch)
phy <- finch$phy
dat <- finch$data[,"wingL"]
unit.tree <- as.unit.tree(phy, data=dat)

## use default statistics
summ.stat <- summ.stats(unit.tree, stats=NULL)
summ.stat

## user defined statistics
mean.con <- function(x) mean(x$pics[,"contrasts"])
max.con <- function(x) max(x$pics[,"contrasts"])

summ.stat.user <- summ.stats(unit.tree,
                    stats=list(mean = mean.con, max = max.con))
summ.stat.user
}
\seealso{
  \code{\link{def.summ.stats}}, \code{\link{sigsq.est}},
  \code{\link{var.contrast}},
  \code{\link{cor.contrast.var}},
  \code{\link{cor.contrast.nh}},
  \code{\link{cor.contrast.asr}}, \code{\link{ks.contrast}}
}

