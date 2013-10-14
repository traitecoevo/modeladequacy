#' Arbutus: Evaluate model adequacy for phylogenetic models of continuous character evolution
#'
#' @description Arbutus was developed to evaluate the adequacy (or absolute goodness-of-fit)
#' of phylogenetic models of continuous character evolution. The basic principle
#' underlying the approach is that if the generating model is a Brownian motion
#' process, the contrasts (sensu Felsenstein, 1985) will be independent and
#' identically distributed (I.I.D.). We can evaluate this condition by calculating
#' a set of summary statistics on our observed data, then simulating datasets under
#' a Brownian motion process along the phylogeny and calculating the summary statistics
#' on each simulated data set. We can then compare our observed summary statistics to the
#' simulated distribution of summary statistics.
#'
#' The I.I.D. property of the contrasts will hold only for the case of Brownian motion.
#' However, the approach can be extended to any arbitrarily complex model of continuous character
#' evolution (as long as the model is based on a multivariate normal distribution) by first
#' rescaling the phylogeny based on the fitted parameter estimates, creating what we refer to
#' as a 'unit.tree'. Once the rescaling is done, the contrasts will again be I.I.D. if the model
#' that was fit is a good fit to the data; which again we can evaluate with the use of
#' summary statistics.
#'
#' The basic functionality of arbutus can be employed using the function \code{\link{phy.model.check}}.
#' This function wraps several other functions which can be used separately:
#' \enumerate{
#'  \item{\code{\link{as.unit.tree}}: }{rescales phylogeny based on fitted parameter values}
#'  \item{\code{\link{summ.stats}}: }{calculates summary statistics on observed data}
#'  \item{\code{\link{sim.char.unit}}: }{simulates datasets under Brownian motion}
#'  \item{\code{\link{summ.stats}}: }{calculates summary statistics on simulated data}
#'  \item{\code{\link{compare.summ.stats}}: }{compares observed to simulated statistics}
#' }
#'
#' Arbutus is designed to interact with a number of other packages for fitting trait models. The
#' architecture is flexible such that additional packages can be incorporated relatively easily.
#'
#' @references Matthew W. Pennell, Richard G. FitzJohn, William K. Cornwell and Luke J. Harmon
#' (in prep.).  On the adequacy of phylogenetic trait models.
#'
#' @references Joseph Felsenstein (1985). Phylogenies and the comparative method. The American
#' Naturalist 125:1-15.
#' 
#' @docType package
#' @import ape ggplot2
#' @name arbutus
#' @aliases arbutus-arbutus arbutus
#'
NULL







#' Galapagos finch dataset
#'
#' A dataset including a phylogeny of 13 galapagos finches with
#' trait data matched to the phylogeny.
#'
#' \itemize{
#'  \item phy. the finch phylogeny
#'  \item data. measured traits for each species in the phylogeny
#' }
#'
#' @docType data
#' @keywords datasets
#' @format a list with two elements
#' @name finch
NULL
