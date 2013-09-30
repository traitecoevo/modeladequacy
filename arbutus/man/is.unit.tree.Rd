\name{is.unit.tree}
\alias{is.unit.tree}
\title{Check if object is a unit.tree}
\usage{
  is.unit.tree(x)
}
\arguments{
  \item{x}{object}
}
\value{
  logical
}
\description{
  Utility function for checking if object is of class
  'unit.tree'
}
\examples{
## finch data
data(geospiza)
td <- suppressWarnings(treedata(geospiza$phy, geospiza$dat))
phy <- td$phy
data <- td$data[,"wingL"]


## using just the given phylogeny
unit.tree.phy <- as.unit.tree(phy, data)

is.unit.tree(unit.tree.phy)
}
\seealso{
  \code{\link{as.unit.tree}}
}

