make_synonyms <- function(plantlist, tree) {
  dat <- read.csv(plantlist, stringsAsFactors=FALSE)
  ## Drop an extra column that has turned up for some reason:
  dat <- dat[c("synonym", "species", "genus")]
  dat <- dat[dat$synonym != dat$species,]
  dat$synonym[dat$synonym %in% tree$tip.label] <- NA
  dat
}

make_corrections <- function(filename, tree) {
  corrections <- read.delim(filename, stringsAsFactors=FALSE)
  corrections[!(corrections$originalName %in%
                sub("_"," ",tree$tip.label)),]
}

update_synonomy <- function(names, syn) {
  ## TODO: What is going on here?
  temp <- syn$species[match(names, syn$synonym)]
  names[names %in% syn$synonym] <- temp[!is.na(temp)]
  names
}

scrub_wrapper<-function(names, corrections) {
  ## saves time by only looking for unique species names
  sp.list <- unique(names)
  ## sp_list <- gsub(" ", "_", sp.list) # TODO: this was never used!
  ## correcting names and dropping to binomial
  out_spp <- scrub(sp.list, corrections)
  out_spp <- gsub(" ", "_", out_spp, fixed=TRUE)
  temp <- out_spp[match(names, sp.list)] #working out the many-to-one mapping
  names[names %in% sp.list] <- as.character(temp)
  names
}

## This code was written by Ginger Jui -- original header below.
## Slight modifications by Rich and Will

## #-------------------------------------------------------------------#
## # ScrubNames
## #
## # Name scrubbing function to clean up species names
## # from plant trait databases
## # 
## # Written by Ginger Jui
## # 20110304
## # 
## # Args:
## #	sp.names: vector of species names as strings
## #	lookup: dataframe of manually scrubbed names with two columns,
## #		"Original" for the original names in the datasets and 
## #		"Changed" for the manually cleaned binomial
## #	
## # Returns:
## #	vector of species names (Genus species separated by a space)
## # 	in same order as original vector; junk names are replaced with NA	
## #-------------------------------------------------------------------#

scrub <- function(sp.names, lookup) {
  ## Replace known problem names by matching to lookup file
  #for (i in 1:nrow(lookup)) {
  #  a <- which(sp.names %in% lookup$Original[i])
  #  sp.names[a] <- lookup$Changed[i]
  #}
  ## This may be faster:
  i <- match(sp.names,lookup$originalName)
  sp.names[!is.na(i)] <- lookup$newName[i[!is.na(i)]]

  ## ###################################################
  ## FILTER 1: Identify common mistakes and            #
  ## fix names in sp.names                             #
  ## ###################################################

  ##---------------------------------------------------#
  ## Turn all underscores and >1 space between names
  ## into single space
  ##---------------------------------------------------#
  re.underscore <- "[_ \n]+"
  fix.underscore <- grepl(re.underscore, sp.names)
  sp.names[fix.underscore] <-
    gsub(re.underscore," ", sp.names[fix.underscore])

  ##---------------------------------------------------#
  ## Delete cf collapsing to binomial
  ## Delete aff. or sp. aff., collapsing to binomial
  ##---------------------------------------------------#
  re.aff <- "(.*) (cf\\.?|aff\\.?|sp\\.? aff\\.?) (.*)"
  fix.aff <- grepl(re.aff, sp.names)
  sp.names[fix.aff] <- gsub(re.aff, "\\1 \\3", sp.names[fix.aff])

  re.cf2 <- "(.*)(\\.cf.*\\.)(.*)"
  fix.cf2 <- grepl(re.cf2, sp.names)
  sp.names[fix.cf2] <- gsub(re.cf2, "\\1 \\3", sp.names[fix.cf2])

  ##---------------------------------------------------#
  ## Remove everything following a var. subsp. ssp. subvar. 
  ##---------------------------------------------------#
  rm.str1 <- paste("\\.*",c("var", "subsp", "ssp", "sp", "subvar",
                            "agg", "cv"), " ", sep="")
  rm.str2 <- paste("\\.*",c("var", "subsp", "ssp", "sp", "subvar",
                            "agg", "cv"), "\\.", sep="")
  rm.str <- c(rm.str1, rm.str2)
  re.ssp <- paste("^(.*) (", paste(rm.str, collapse="|"), ").*$", sep="")
  fix.ssp <- grepl(re.ssp, sp.names)
  sp.names[fix.ssp] <- sub(re.ssp, "\\1", sp.names[fix.ssp])

  ##---------------------------------------------------#
  ## Bionomials with alternate genus in parentheses
  ## Dump things in parentheses, collapsing to binomial
  ##---------------------------------------------------#

  ## ending in parentheses
  re.end.in.par <- " ?\\(.*\\)$"
  fix.end.in.par <- grepl(re.end.in.par, sp.names)
  sp.names[fix.end.in.par] <-
    sub(re.end.in.par, "", sp.names[grepl(re.end.in.par, sp.names)])

  ## parentheses embedded in name
  re.syn <- "(.*)(\\(.+\\))(.+)$"
  fix.syn <- grepl(re.syn, sp.names)
  sp.names[fix.syn] <- gsub(re.syn, "\\1 \\3", sp.names[fix.syn])

  ## Turn all spaces between names into single space
  sp.names <- gsub("( )+"," ", sp.names)
 
  ##---------------------------------------------------#
  ## Extract potential binomials, trinomials, etc.
  ## and collapse to binomails
  ##---------------------------------------------------#

  ## All binomials or trinomials
  re.bi.tri <- "(^[A-Z][a-z]+) ([-a-z]+)( .*)+$"
  fix.bi.tri <- grepl(re.bi.tri, sp.names)
  sp.names[fix.bi.tri] <-
    gsub(re.bi.tri, "\\1 \\2", sp.names[fix.bi.tri])

  ## Nuke trailing whitespace:
  sp.names <- sub(" +$", "", sp.names)
  ## Nuke trailing periods from possible species names (not
  ## abbreviations)
  sp.names <- sub("([a-z]{4,})\\.+$", "\\1", sp.names)
  
  ## ###################################################
  ## FILTER 2: Identify names to exclude               #
  ## ###################################################
  ## Logical vector for tracking which names are good binomials
  good.names <- rep(FALSE, length(sp.names))
  bad.names <- rep(FALSE, length(sp.names))

  re.potentials <- "^[A-Z][a-z]+(-[A-Z][a-z]+)? [-a-z]+( .*)*$"
  good.names[grep(re.potentials, sp.names)] <- TRUE

  ## Exclude species name "indet"
  ex.indet <- grepl("[A-Za-z]+ indet", sp.names)
  bad.names[ex.indet] <- TRUE
  
      ## Exclude species name "nondet"
  ex.nondet <- grepl("[A-Za-z]+ nondet", sp.names)
  bad.names[ex.nondet] <- TRUE
  
    ## Exclude species name "sect"
  ex.sect <- grepl("[A-Za-z]+ sect", sp.names)
  bad.names[ex.sect] <- TRUE
  
      ## Exclude species name "x"
  ex.hybrid <- grepl("[A-Za-z]+ x$", sp.names)
  bad.names[ex.hybrid] <- TRUE

  ## Exclude species that start with parentheses
  re.start.in.par <- "^\\((.*)$"
  fix.start.in.par <- grepl(re.start.in.par, sp.names)
  bad.names[fix.start.in.par] <- TRUE

  ## Exclude rows with species name that end in sp. or sp
  ## and followed by numbers
  re.sp <- "(.+) (sp.|sp|spp.|sp[0-9]+|sp.[0-9]+)( |$)"
  fix.sp <- grepl(re.sp, sp.names)
  bad.names[fix.sp] <- TRUE

  ## Exclude species names that include question marks
  re.quest <- "\\?"
  fix.quest <- grepl(re.quest, sp.names) 
  bad.names[fix.quest] <- TRUE

  ## Exclude hybrids
  re.x <- " (x|X) "
  fix.x <- grepl(re.x, sp.names)
  bad.names[fix.x] <- TRUE

  ## Exclude mono names:
  bad.names[!grepl(" ", sp.names)] <- TRUE
  ## Exclude anything with a number
  bad.names[grepl("[0-9]", sp.names)] <- TRUE
  
 

  sort(unique(sp.names[!good.names & !bad.names]))

  sp.names[!(good.names & !bad.names)] <- NA
  
   ##match with lookup table again
 i <- match(sp.names,lookup$originalName)
  sp.names[!is.na(i)] <- lookup$newName[i[!is.na(i)]]
  
  sp.names
}
