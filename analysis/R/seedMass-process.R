source("R/load-scrubbing-tools.R")

kew<-read.csv("data/kew.csv",as.is=T)

#cleaning names
kew$species<-scrub.wrapper(kew$species)



#correcting errors (could be a seprate file, but it's only a few lines)
kew.errors<-subset(errors,errors$Dataset=="kewSeed"&errors$trait=="seedMass")
kew.errors$Original<-as.numeric(kew.errors$Original)
kew.errors$Changed<-as.numeric(kew.errors$Changed)
kew.errors$genus_species<-gsub(" ","_",kew.errors$genus_species)
#double column matching
replace.matrix<-which( outer(kew.errors$genus_species, kew$species, "==") & 
                         outer(round(kew.errors$Original), round(kew$value), "=="), 
                       arr.ind=TRUE)
#replace bad data points
kew$value[replace.matrix[,2]]<-kew.errors$Changed[replace.matrix[,1]]

#using modified plant list synonmy
temp<-pl.mod$species[match(kew$species,pl.mod$synonym)]
kew$species[kew$species%in%pl.mod$synonym]<-temp[!is.na(temp)]


kew<-subset(kew,!is.na(kew$value))
seed.mean<-10^(tapply(log10(kew$value),as.factor(kew$species),FUN=mean))
seed.spp.sd<-tapply(log10(kew$value),as.factor(kew$species),FUN=sd,na.rm=T)
sd.real<-seed.spp.sd[seed.spp.sd>0.000001]
mean(sd.real,na.rm=T)
#hist(log10(seed.mean))

write.csv(seed.mean,"output/species-mean-seedMass.csv")


#Code below here need manual supervision
#False positive rate is too high without checking


# still.missing<-names(seed.mean)[!names(seed.mean)%in%tree$tip.label]
# hello<-sapply(still.missing[10001:11765],agrep.for.names,good.names=tree$tip.label)
# 

