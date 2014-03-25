## These functions are designed to read and process the data
require(geiger, quietly=TRUE)

build.data <- function(dataset) {
  dataset <- match.arg(dataset, c("sla", "seedMass", "leafN"))

  t <- get.tree()
  raw <- read.csv(sprintf("output/species-mean-%s.csv", dataset),
                  stringsAsFactors=FALSE)
  t <- extract.clade(t, node="Angiospermae")

  ## NOTE: base 10 log
  dat <- structure(log10(raw$mean), names=raw$gs)

  ## Drop species from tree not in data, and v.v.
  ## This step is the slow point.
  phy <- geiger:::.drop.tip(phy=t, tip=setdiff(t$tip.label, names(dat)))
  dat <- dat[phy$tip.label]

  ## TODO: recompute SLA from data?
  ##   mean(log10(raw$sd[log10(raw$sd) > 0.0001]), na.rm=TRUE)
  se <- c(sla=0.1039405, seedMass=0.1551108, leafN=0.07626127)[[dataset]]

  list(phy=phy, states=dat, SE=se)
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

get.sla.v.leafn.data <-function(){
    sla <- get.data("SLA")
    ln <- get.data("leafn")

    tmp <- sla$phy$tip.label[!(sla$phy$tip.label %in% ln$phy$tip.label)]
    phy <- geiger:::.drop.tip(phy=sla$phy, tip=tmp)

    states.sla <- sla$states[phy$tip.label] 
    states.ln <- ln$states[phy$tip.label]
    states <- cbind.data.frame(states.sla, states.ln)
    colnames(states) <- c("sla", "leafN")

    SE <- c(sla$SE, ln$SE)
    names(SE) <- c("sla", "leafN")

    list(phy=phy, states=states, SE=SE)
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
