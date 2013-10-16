source("R/load-scrubbing-tools.R")

#load_raw_data
wright<-read.csv("data/wright-2004.csv",as.is=T)


options(stringsAsFactors = FALSE)
trait.raw<-data.frame(gs=c(wright$Species),Nmass=c(wright$N.mass))

#cleaning names
trait.raw$gs<-scrub.wrapper(trait.raw$gs)

#correcting errors (could be a seprate file, but it's only a few lines)
#no known errors in glopnet

#using modified plant list synonmy
temp<-pl.mod$species[match(trait.raw$gs,pl.mod$synonym)]
trait.raw$gs[trait.raw$gs%in%pl.mod$synonym]<-temp[!is.na(temp)]


leafN.subset<-subset(trait.raw,!is.na(trait.raw$Nmass))
#geometric mean for species
leafN.spp.mean<-10^(tapply(log10(leafN.subset$Nmass),as.factor(leafN.subset$gs),FUN=mean))
leafN.spp.sd<-tapply(log10(leafN.subset$Nmass),as.factor(leafN.subset$gs),FUN=sd,na.rm=T)
sd.real<-leafN.spp.sd[leafN.spp.sd>0.0001]
mean(sd.real,na.rm=T)

write.csv(leafN.spp.mean,"output/species-mean-leafN.csv")

#Code below here need manual supervision
#False positive rate is too high without checking

# 
#  full.names<-gsub(" ","_",names(leafN.spp.mean))
#   still.missing<-full.names[!full.names%in%tree$tip.label]
#   hello<-sapply(still.missing,agrep.for.names,good.names=tree$tip.label)
# this found about 15-20 more legit misspellings + many false positives.

