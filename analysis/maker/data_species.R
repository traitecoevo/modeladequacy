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

  ## Here, what we want is four possible things:
  ##   * Mean and sd of log(trait)
  ##   * Mean and sd of trait
  dat_spp <-
    dat[complete.cases(dat),] %>%
      group_by(gs)            %>%
        summarise(mean_log = mean(log(N.mass)),
                  sd_log   = sd(log(N.mass)),
                  mean_raw = mean(N.mass),
                  sd_raw   = sd(N.mass),
                  n_obs    = length(N.mass))

  ## From the paper:
  ## > we estimated a single SE for each trait by calculating the mean
  ## > standard deviation for all species for which we had multiple
  ## > measurements
  ok <- dat_spp$n_obs > 1L
  attr(dat_spp, "sd_raw") <- mean(dat_spp$sd_raw[ok])
  attr(dat_spp, "sd_log") <- mean(dat_spp$sd_log[ok])

  ## Remember the trait for later:
  attr(dat_spp, "trait") <- "leafN"
  
  dat_spp
}

make_species_seed_mass <- function(kew, synonyms, corrections) {
  dat <- data.frame(gs=kew$species,
                    seed_mass=kew$value,
                    stringsAsFactors=FALSE)

  dat$gs <- scrub_wrapper(kew$species, corrections)

  browser(skipCalls=2L)

  ## Correct errors:
  errors <- read.csv("data/errors.csv", stringsAsFactors=FALSE)
  kew_errors <- errors[errors$Dataset=="kewSeed" &
                       errors$trait=="seed_mass",]
  kew_errors$Original <- as.numeric(kew_errors$Original)
  kew_errors$Changed  <- as.numeric(kew_errors$Changed)
  kew_errors$genus_species <- gsub(" ", "_", kew_errors$genus_species)
  #double column matching

  replace_matrix <- which(
    outer(kew_errors$genus_species, dat$gs, "==") &
    outer(round(kew_errors$Original), round(dat$seed_mass), "=="),
    arr.ind=TRUE)
  ## TODO: Again, only 21 of 60 match.
  dat$seed_mass[replace_matrix[,2]] <-
    kew_errors$Changed[replace_matrix[,1]]

  #using modified plant list synonmy
  dat$gs <- update_synonomy(dat$gs, synonyms)

  ## Build a little data frame with the species names and geometric mean
  ## of the trait:
  dat_spp <-
    dat[complete.cases(dat),] %>%
      group_by(gs)            %>%
        summarise(mean_log = mean(log(seed_mass)),
                  sd_log   = sd(log(seed_mass)),
                  mean_raw = mean(seed_mass),
                  sd_raw   = sd(seed_mass),
                  n_obs    = length(seed_mass))

  ## From the paper:
  ## > we estimated a single SE for each trait by calculating the mean
  ## > standard deviation for all species for which we had multiple
  ## > measurements
  ok <- dat_spp$n_obs > 1L
  attr(dat_spp, "sd_raw") <- mean(dat_spp$sd_raw[ok])
  attr(dat_spp, "sd_log") <- mean(dat_spp$sd_log[ok])

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
        summarise(mean_log = mean(log(sla)),
                  sd_log   = sd(log(sla)),
                  mean_raw = mean(sla),
                  sd_raw   = sd(sla),
                  n_obs    = length(sla))

  ## From the paper:
  ## > we estimated a single SE for each trait by calculating the mean
  ## > standard deviation for all species for which we had multiple
  ## > measurements
  ok <- dat_spp$n_obs > 1L
  attr(dat_spp, "sd_raw") <- mean(dat_spp$sd_raw[ok])
  attr(dat_spp, "sd_log") <- mean(dat_spp$sd_log[ok])
  
  attr(dat_spp, "trait") <- "sla"

  dat_spp
}
