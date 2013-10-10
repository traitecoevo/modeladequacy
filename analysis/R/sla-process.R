source("R/load-scrubbing-tools.R")

#load_raw_data
wright<-read.csv("data/wright-2004.csv",as.is=T)
leda<-read.csv("data/leda.csv",as.is=T)

options(stringsAsFactors = FALSE)
trait.raw<-data.frame(gs=c(wright$Species,leda$SBS.name),sla=c(10000/10^(wright$LogLMA),10*leda$SLA.mean),Nmass=c(wright$N.mass,rep(NA,length(leda$SLA.mean))),dataset=c(rep("glop",length(wright$N.mass)),rep("leda",length(leda$SLA.mean))))

#cleaning names
trait.raw$gs<-scrub.wrapper(trait.raw$gs)


#correcting errors (could be a seprate file, but it's only a few lines)
#only errors in LEDA, none in glopnet
leda.errors<-subset(errors,errors$Dataset=="LEDA"&errors$trait=="sla")
leda.errors$Original<-as.numeric(leda.errors$Original)
leda.errors$Changed<-as.numeric(leda.errors$Changed)
#double column matching
replace.matrix<-which( outer(leda.errors$genus_species, trait.raw$gs, "==") & 
       outer(round(leda.errors$Original), round(trait.raw$sla), "=="), 
       arr.ind=TRUE)
#replace bad data points
trait.raw$sla[replace.matrix[,2]]<-leda.errors$Changed[replace.matrix[,1]]

#using modified plant list synonmy
temp<-pl.mod$species[match(trait.raw$gs,pl.mod$synonym)]
trait.raw$gs[trait.raw$gs%in%pl.mod$synonym]<-temp[!is.na(temp)]


sla.subset<-subset(trait.raw,!is.na(trait.raw$sla))
#geometric mean for species
sla.spp.mean<-10^(tapply(log10(sla.subset$sla),as.factor(sla.subset$gs),FUN=mean))
sla.spp.sd<-tapply(log10(sla.subset$sla),as.factor(sla.subset$gs),FUN=sd,na.rm=T)
sd.real<-sla.spp.sd[sla.spp.sd>0.0001]
mean(sd.real,na.rm=T)

write.csv(sla.spp.mean,"output/species-mean-sla.csv")

#Code below here need manual supervision
#False positive rate is too high without checking

# 
# trait.raw$gs<-gsub(" ","_",trait.raw$gs)
# sla.subset<-subset(trait.raw,!is.na(trait.raw$sla))
#  still.missing<-sla.subset$gs[!sla.subset$gs%in%tree$tip.label]
#  hello<-sapply(still.missing,agrep.for.names,good.names=tree$tip.label)
#this found about 15-20 more legit misspellings + many false positives.