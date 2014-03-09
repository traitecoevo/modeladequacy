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

get.sla.data <- function() {
  readRDS("output/data-sla.rds")
}

get.seedmass.data <- function() {
  readRDS("output/data-seedMass.rds")
}

get.leafn.data <- function() {
  readRDS("output/data-leafN.rds")
}

get.sla.v.leafn.data <-function(){
    sla <- get.sla.data()
    ln <- get.leafn.data()

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
