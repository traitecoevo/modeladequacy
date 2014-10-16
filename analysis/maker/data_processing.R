## Converts a raw data set into one with the tree associated with it.
build_data <- function(raw, tree) {
  ## NOTE: base 10 log
  dat <- structure(log10(raw$mean), names=raw$gs)

  ## Drop species from tree not in data, and v.v.
  ## This step is the slow point.
  t <- extract.clade(tree, node="Angiospermae")
  phy <- geiger:::.drop.tip(phy=t, tip=setdiff(t$tip.label, names(dat)))
  dat <- dat[phy$tip.label]
  trait <- attr(raw, "trait")
  se <- attr(raw, "sd")

  list(phy=phy, states=dat, se=se, trait=trait)
}


## The more serious data munging: building the leafN, seed_mass and
## sla data sets:
make_species_leafN <- function(wright_2004, synonyms, corrections) {
  dat <- wright_2004[c("Species", "N.mass")]
  names(dat) <- c("gs", "N.mass")

  dat$gs <- scrub_wrapper(dat$gs, corrections)

  ## Correct errors:
  ##
  ## [no known errors in glopnet, so nothing here]

  dat$gs <- update_synonomy(dat$gs, synonyms)

  ## Build a little data frame with the species names and geometric mean
  ## of the trait:
  dat_spp <-
    dat[complete.cases(dat),] %>%
      group_by(gs)            %>%
        summarise(mean = geometric_mean(N.mass),
                  sd   = geometric_sd(N.mass))

  sd <- mean(log10(dat_spp$sd[log10(dat_spp$sd) > 0.0001]), na.rm=TRUE)
  attr(dat_spp, "sd") <- sd
  attr(dat_spp, "trait") <- "leafN"
  
  dat_spp
}

make_species_seed_mass <- function(kew, synonyms, corrections) {
  dat <- data.frame(gs=kew$species,
                    seedMass=kew$value,
                    stringsAsFactors=FALSE)

  dat$gs <- scrub_wrapper(kew$species, corrections)

  ## Correct errors:
  errors <- read.csv("data/errors.csv", stringsAsFactors=FALSE)
  kew_errors <- errors[errors$Dataset=="kewSeed" &
                       errors$trait=="seedMass",]
  kew_errors$Original <- as.numeric(kew_errors$Original)
  kew_errors$Changed  <- as.numeric(kew_errors$Changed)
  kew_errors$genus_species <- gsub(" ", "_", kew_errors$genus_species)
  #double column matching

  replace_matrix <- which(
    outer(kew_errors$genus_species, dat$gs, "==") &
    outer(round(kew_errors$Original), round(dat$seedMass), "=="),
    arr.ind=TRUE)
  ## TODO: Again, only 21 of 60 match.
  dat$seedMass[replace_matrix[,2]] <-
    kew_errors$Changed[replace_matrix[,1]]

  #using modified plant list synonmy
  dat$gs <- update_synonomy(dat$gs, synonyms)

  ## Build a little data frame with the species names and geometric mean
  ## of the trait:
  dat_spp <-
    dat[complete.cases(dat),] %>%
      group_by(gs)            %>%
        summarise(mean = geometric_mean(seedMass),
                  sd   = geometric_sd(seedMass))

  sd <- mean(log10(dat_spp$sd[log10(dat_spp$sd) > 0.0001]), na.rm=TRUE)
  attr(dat_spp, "sd") <- sd
  attr(dat_spp, "trait") <- "seed_size"

  dat_spp
}

make_species_sla <- function(wright_2004, leda, synonyms, corrections) {
  ## TODO: Check $logLMA is not character
  ## TODO: We have unlogged LMA here...
  wright_2004$LogLMA <- as.numeric(wright_2004$LogLMA)
  wright_2004 <- data.frame(gs=wright_2004$Species,
                            sla=10000/10^(wright_2004$LogLMA),
                            dataset="glop",
                            stringsAsFactors=FALSE)
  leda <- data.frame(gs=leda$SBS.name,
                     sla=10 * leda$SLA.mean,
                     dataset="leda",
                     stringsAsFactors=FALSE)
  dat <- rbind(wright_2004, leda)

  ## Sort out synonomy and mispellings in species names.
  dat$gs <- scrub_wrapper(dat$gs, corrections)

  ## Correct errors:
  ##
  errors <- read.csv("data/errors.csv", stringsAsFactors=FALSE)
  leda_errors <- errors[errors$Dataset=="LEDA" & errors$trait=="sla",]
  leda_errors$Original <- as.numeric(leda_errors$Original)
  leda_errors$Changed  <- as.numeric(leda_errors$Changed)
  leda_errors$gs <- sub(" ", "_", leda_errors$genus_species, fixed=TRUE)

  # double column matching

  ## Hmm, looks like this was not working.  I needed to do the
  ## underscore substitute, but even then I don't get all the matches I
  ## need.  I'd expect 21 here, but only get 10.
  ##
  ## TODO: talk through this with Matt, replace with something more
  ## robust?
  replace_matrix <- which(
    outer(leda_errors$gs, dat$gs, "==") &
    outer(round(leda_errors$Original), round(dat$sla), "=="),
    arr.ind=TRUE)
  dat$sla[replace_matrix[,2]] <-
    leda_errors$Changed[replace_matrix[,1]]

  dat$gs <- update_synonomy(dat$gs, synonyms)

  ## Build a little data frame with the species names and geometric mean
  ## of the trait:
  dat_spp <-
    dat[complete.cases(dat),] %>%
      group_by(gs)            %>%
        summarise(mean = geometric_mean(sla),
                  sd   = geometric_sd(sla))

  sd <- mean(log10(dat_spp$sd[log10(dat_spp$sd) > 0.0001]), na.rm=TRUE)
  attr(dat_spp, "sd") <- sd
  attr(dat_spp, "trait") <- "sla"

  dat_spp
}

## Used by both make_data_times and make_data_rank
add_metadata <- function(d, rank, taxa, trait, se) {
  d$rank <- rank
  d$taxa <- taxa
  d$trait <- trait
  d$se <- se
  d
}
