## These functions are designed to read and process the data
require(geiger, quietly=TRUE)

path.data <- function() {
  "output/angio-data"
}
path.ml <- function() {
  "output/results-ml"
}
path.bayes <- function() {
  "output/results-bayes"
}

get.tree <- function() {
  readRDS("output/vascular_plant_phylogeny.rds")
}

get.synonyms <- function() {
  readRDS("output/synonyms.rds")
}

get.corrections <- function() {
  readRDS("output/corrections.rds")
}

get.errors <- function() {
  read.csv("data/errors.csv", stringsAsFactors=FALSE)
}

get.data <- function(trait) {
  trait <- match.arg(trait, c("SLA", "leafn", "seedmass"))
  # TODO: This needs dealing with throughout.  One capitalisation for
  # each please.
  tr <- c(SLA="sla", leafn="leafN", seedmass="seedMass")
  readRDS(sprintf("output/data-%s.rds", tr[[trait]]))
}

filename.analysis <- function(trait, type, age=NULL, idx=NULL) {
  if (type == "timeslice") {
    if (is.null(age) || is.null(idx))
      stop('age and idx must be provided for type="timeslice"')
    type <- paste(type, age, idx, sep="_")
  } else {
    type <- paste("clade", type, sep="_")
  }
  sprintf("%s_%s.rds", trait, type)
}
