library(ape)
wright<-read.csv("../data/wright-2004.csv",as.is=T)
leda<-read.csv("../data/leda.csv",as.is=T)

trait.raw<-data.frame(gs=c(wright$Species,leda$SBS.name),sla=c(1000/10^(wright$LogLMA),leda$SLA.mean),Nmass=c(wright$N.mass,rep(NA,length(leda$SLA.mean))))

path.forest <- readLines("~/.forest_path")

#loading tools from forest
plantList<-read.csv(file.path(path.forest,"taxonomic/spermatophyta_synonyms_PLANTLIST.csv"),as.is=TRUE)
corrections<-read.delim(file.path(path.forest,"db/names-tr.txt"),as.is=TRUE)
errors<-read.csv(file.path(path.forest,"db/errors.csv"),as.is=TRUE)
errors<-read.csv(file.path(path.forest,"db/errors.csv"),as.is=TRUE)
tree<-read.tree(file.path(path.forest,"/trees/tempomode_trees_02112013/tempo_scrubbed_CONSTRAINT_rooted.dated.tre"))
source(file.path(path.forest,"R/import-scrub.R"))


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
still.missing<-gsub(" ","_",still.missing)
fixed.names2<-plantList$species[match(still.missing,plantList$synonym)]
still.missing2<-still.missing[!fixed.names2%in%tree$tip.label]
sum(fixed.names2%in%tree$tip.label)
#fixes 97 species



still.missing.3<-sort(still.missing2[!still.missing2%in%plantList$synonym&!is.na(still.missing2)])
agrep.for.names<-function(good.names,bad.name){
  fix<-good.names[agrep(bad.name,good.names,max.distance=0.1)]
 if (length(fix>0)){
  print(paste("bad name:",bad.name))
  print("fixes:")
  print(fix)
 }
}

hello<-sapply(still.missing.3,agrep.for.names,good.names=tree$tip.label)
#this found about 15-20 more legit misspellings + many false positives.