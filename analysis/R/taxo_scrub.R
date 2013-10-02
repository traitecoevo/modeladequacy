library(ape)
#wright<-read.csv("../data/wright-2004.csv",as.is=T)
#leda<-read.csv("../data/leda.csv",as.is=T)

#trait.raw<-data.frame(gs=c(wright$Species,leda$SBS.name),sla=c(10000/10^(wright$LogLMA),10*leda$SLA.mean),
 #                     Nmass=c(wright$N.mass,rep(NA,length(leda$SLA.mean))),dataset=c(rep("glop",length(wright$N.mass)),rep("leda",length(leda$SLA.mean))))
setwd("/Users/willcornwell/Documents/modeladequacy/analysis/R")
trait.raw<-read.csv("species_mean_sla.csv")


#loading tools from forest
path.forest <- readLines("~/.forest_path")
plantList<-read.csv(file.path(path.forest,"taxonomic/spermatophyta_synonyms_PLANTLIST.csv"),as.is=TRUE)
corrections<-read.delim(file.path(path.forest,"db/names-tr.txt"),as.is=TRUE)
tree<-read.tree(file.path(path.forest,"/trees/tempomode_trees_02112013/tempo_scrubbed_CONSTRAINT_rooted.dated.tre"))
source(file.path(path.forest,"R/import-scrub.R"))

#
#correcting errors (could be a seprate file, but it's only a few lines)
#only errors in LEDA, none in glopnet
errors<-read.csv(file.path(path.forest,"db/errors.csv"),as.is=TRUE)
leda.errors<-subset(errors,errors$Dataset=="LEDA"&errors$trait=="sla")
leda.errors$Original<-as.numeric(leda.errors$Original)
leda.errors$Changed<-as.numeric(leda.errors$Changed)
#double column matching
replace.matrix<-which( outer(leda.errors$genus_species, trait.raw$gs, "==") & 
       outer(round(leda.errors$Original), round(trait.raw$sla), "=="), 
       arr.ind=TRUE)
#replace bad data points
trait.raw$sla[replace.matrix[,2]]<-leda.errors$Changed[replace.matrix[,1]]


#now taxon scrubbing
sp.list<-unique(gsub(" ","_",trait.raw$gs))
good.names<-sp.list[sp.list%in%tree$tip.label]
missing.names<-sp.list[!sp.list%in%tree$tip.label]

#trying ginger's script
out.spp<-scrub(missing.names,corrections)
still.missing<-missing.names[missing.names==out.spp]
out.spp<-gsub(" ","_",out.spp)
fixed.names1<-out.spp[out.spp%in%tree$tip.label]
sum(fixed.names1%in%tree$tip.label)
#fixes 162 species

#using plant list synonmy
new.missing<-out.spp[!out.spp%in%fixed.names1]
potential.fixes<-plantList$species[match(new.missing,plantList$synonym)]
fixed.names2<-potential.fixes[potential.fixes%in%tree$tip.label]
#fixes 101 species

look.up.table<-data.frame(original.names=c(missing.names[which(out.spp%in%tree$tip.label)]
                                         ,new.missing[which(potential.fixes%in%tree$tip.label)])
                                         ,fixed.names=c(fixed.names1,
                                         fixed.names2))
trait.raw$gs<-gsub(" ","_",trait.raw$gs)

temp<-look.up.table$fixed.names[match(trait.raw$gs,look.up.table$original.names)]
temp<-as.character(temp[!is.na(temp)])
trait.raw$gs[trait.raw$gs%in%look.up.table$original.names]<-temp

sla.subset<-subset(trait.raw,!is.na(trait.raw$sla))
#geometric mean for species
sla.spp.mean<-exp(tapply(log(sla.subset$sla),as.factor(sla.subset$gs),FUN=mean))
write.csv(sla.spp.mean,"species_mean_sla.csv")

#Code below here need manual supervision
#False positive rate is too high without checking

# agrep.for.names<-function(good.names,bad.name){
#   fix<-good.names[agrep(bad.name,good.names,max.distance=0.08)]
#  if (length(fix>0)){
#   print(paste("bad name:",bad.name))
#   print("fixes:")
#   print(fix)
#  }
# }
# 
# still.missing<-sort(look.up.table$missing.names[is.na(look.up.table$fixed.names)])
# hello<-sapply(still.missing,agrep.for.names,good.names=tree$tip.label)
#this found about 15-20 more legit misspellings + many false positives.