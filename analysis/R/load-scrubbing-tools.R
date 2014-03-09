require(ape, quietly=TRUE)
source("R/import-scrub.R")
source("R/read-data-functions.R")

scrub.wrapper<-function(names){
  corrections <- get.corrections()
  #saves time by only looking for unique species names
  sp.list<-unique(names)
  sp_list<-gsub(" ","_",sp.list)
  out.spp<-scrub(sp.list,corrections) #correcting names and dropping to binomial
  out_spp<-gsub(" ","_",out.spp) 
  temp<-out_spp[match(names,sp.list)] #working out the many-to-one mapping
  names[names%in%sp.list]<-as.character(temp)
  names
}

agrep.for.names<-function(good.names,bad.name){
  fix<-good.names[agrep(bad.name,good.names,max.distance=0.03)]
  if (length(fix>0)&length(fix)<5){
    print(paste("bad name:",sub("_"," ",bad.name)))
    print("fixes:")
    print(sub("_"," ",fix))
    #out.df<-data.frame(good=fix,bad=bad.name)
    #return(out.df)
  }
}
# 

logspace.f <- function(x, f, ...) {
  exp(f(log(x), ...))
}

geometric.mean <- function(x, ...) {
  logspace.f(x, mean, ...)
}

geometric.sd <- function(x, ...) {
  logspace.f(x, sd, ...)
}
