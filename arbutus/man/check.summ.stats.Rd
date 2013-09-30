\name{check.summ.stats}
\alias{check.summ.stats}
\title{Internal function to check summary stats}
\usage{
  check.summ.stats(stats)
}
\arguments{
  \item{stats}{a named list of summary statistics}
}
\value{
  named list of summary statistics after being checked.
}
\description{
  Makes sure summary statistics sent to
  \code{\link{summ.stats}} are in the correct format
}
\examples{
## produce list of summary stats (using default summary statistics)
check.summ.stats(def.summ.stats())

## use custom list of summary statistics
foo <- function(x) mean(x$pics[,"contrasts"])

my.stats <- list(mean.of.contrasts=foo)
check.summ.stats(my.stats)
}
\author{
  Rich FitzJohn
}
\seealso{
  \code{\link{summ.stats}}
}
\keyword{internal}

