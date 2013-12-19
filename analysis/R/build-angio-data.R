## functions for building rds file for all data and all traits

require(geiger)
require(phytools)

## get the data for the SLA
get.sla.data <- function(){
    t <- read.tree(file="data/tempo_scrubbed_CONSTRAINT_rooted.dated.tre")
    sla.raw <- read.csv(file="output/species-mean-sla.csv", header=TRUE)

    ## only angiosperms
    t <- extract.clade(t, node="Angiospermae")

    ## base 10 log
    sla <- log10(sla.raw[,"x"])
    names(sla) <- sla.raw[,"X"]

    ## drop extra tips
    tmp <- t$tip.label[!t$tip.label %in% names(sla)]
    phy <- geiger:::.drop.tip(phy=t, tip = tmp)

    sla <- sla[phy$tip.label]

    ## return tree and data
    list(phy=phy, states=sla, SE=0.1039405)
}

## get the data for seedmass
get.seedmass.data <- function(){
    t <- read.tree(file="data/tempo_scrubbed_CONSTRAINT_rooted.dated.tre")
    sm.raw <- read.csv(file="output/species-mean-seedMass.csv", header=TRUE)

    ## only angiosperms
    t <- extract.clade(t, node="Angiospermae")

    ## base 10 log
    sm <- log10(sm.raw[,"x"])
    names(sm) <- sm.raw[,"X"]

    ## drop all for which the value is -Inf
    # sm <- sm[-which(sm == -Inf)]

    ## drop extra tips
    tmp <- t$tip.label[!t$tip.label %in% names(sm)]
    phy <- geiger:::.drop.tip(phy=t, tip = tmp)

    sm <- sm[phy$tip.label]
        
    ## return tree and data
    list(phy=phy, states=sm, SE=0.1551108)
}


## get the data for leafn
get.leafn.data <- function(){
    t <- read.tree("data/tempo_scrubbed_CONSTRAINT_rooted.dated.tre")
    ln.raw <- read.csv("output/species-mean-leafN.csv")

    ## only angiosperms
    t <- extract.clade(t, node="Angiospermae")

    ## log base 10
    ln <- log10(ln.raw[,"x"])
    names(ln) <- ln.raw[,"X"]

    ## drop extra tips
    tmp <- t$tip.label[!t$tip.label %in% names(ln)]
    phy <- geiger:::.drop.tip(phy=t, tip = tmp)

    ln <- ln[phy$tip.label]

    ## return tree and data
    list(phy=phy, states=ln, SE=0.07626127) ## need to update SE
}


## functions for extracting subtrees
## TODO: NEED TO SEE WHICH OF THESE I ACTUALLY NEED
extract.all.nodes<-function(tree,node.height.min,node.height.max){
  #note that this returns a matrix of the same dimensions as tree$edge 
  #containing the height above the root of each node in edge
  tree$node.label<-rep("",length(tree$node.label))
  a<-nodeHeights(tree)
  interesting.nodes<-tree$edge[which(a[,2]>node.height.min&a[,2]<node.height.max)] #these are edges which pass through a certain age
  non.terminal.interesting.nodes<-interesting.nodes[interesting.nodes>length(tree$node.label)+2]
  return(non.terminal.interesting.nodes)
}

where.to.cut<-function(tree,age){
  #note that this returns a matrix of the same dimensions as tree$edge 
  #containing the height above the root of each node in edge
  tree$node.label<-rep("",length(tree$node.label))
  a<-nodeHeights(tree)
  interesting.nodes<-tree$edge[which(a[,1]<age&a[,2]>age),2] #these are edges which pass through a certain age
  non.terminal.interesting.nodes<-interesting.nodes[interesting.nodes>length(tree$node.label)+2]
  return(non.terminal.interesting.nodes)
}

find.node.label<-function(tree,edge.matrix.id,plot.tree=TRUE){
  new.node.label<-rep("",length(tree$node.label))
  node.id<-edge.matrix.id-length(tree$tip.label)
  new.node.label[node.id]<-1:length(node.id)
  tree$node.label<-new.node.label
  if(plot.tree) plot(tree,show.tip.label=F,show.node.label=T)
  return(new.node.label)
}

extract.sub.trees<-function(tree,species.richness,poss.nodes){
  # Extracts subclades with a species richness criterion
  # Args:
  #   tree: the big tree.
  #   species.richness: the cut off for which sub-clade to return
  #   poss.nodes: the possible nodes to check
  # Returns:
  #   list of sub-clades that have enough species
  tree.list<-list()  
  tree$node.label<-poss.nodes
  pn<-poss.nodes[poss.nodes!=""]
  tree.list<-lapply(pn,FUN=extract.clade,phy=tree)
  sr<-sapply(tree.list,FUN=function(x)length(x$tip.label))
#   throw error here?
   if (sum(sr>species.richness)<1) {
     stop("no clades within the specified age range have enough species")
   }
  out.tree.list<-tree.list[which(sr>species.richness)]
  return(out.tree.list)
}

time.slice.tree<-function(time.slice,temp.tree,sr){
edge.matrix.id<-where.to.cut(temp.tree,time.slice)
node.label<-find.node.label(temp.tree,edge.matrix.id,FALSE)
extract.sub.trees(tree=temp.tree,species.richness=sr,poss.nodes=node.label)
}

extract.relevant.nodes<-function(tree,sr,node.height.min,node.height.max){
  #
  # Extracts subclades with a species richness criterion
  # Args:
  #   tree: the big tree.
  #   species.richness: the cut off for which sub-clade to return
  #   node.height.min/max: the age range of interest
  # Returns:
  #   list of sub-clades within the age range that have enough species
  #
  #note can extract the species richness later with: sapply(out.tree.list,FUN=function(x) length(x$tip.label))
  #can extract the age of the clade with: sapply(out.tree.list,FUN=function(x) max(nodeHeights(x)))
  #
  edge.matrix.id<-extract.all.nodes(tree,node.height.min,node.height.max)
  node.label<-find.node.label(tree,edge.matrix.id,TRUE)
  out.tree.list<-extract.sub.trees(tree=tree,species.richness=sr,poss.nodes=node.label)
  return(out.tree.list)
}



## function which pulls out data sets by clade from
## full tree and dataset
## assume rownames of dataset are taxon labels
## rank can be 'family' or 'order' (will add genus later)
## min.size is the minimum clade size
treedata.taxon <- function(phy, data, rank="family", min.size){
    
    if (rank == "family"){
        ## get all family nodes in tree
        tax <- phy$node.label[grep("[A-z]+ceae", phy$node.label, perl=TRUE)]
        ## extract subtrees
        trees <- lapply(tax, function(x) extract.clade(phy, node=x))
        ## use treedata to match to family level
        td <- lapply(trees, function(x) treedata(phy=x, data=data))
        names(td) <- tax
        ## get number of tips
        tips <- lapply(td, function(x) Ntip(x$phy))
        ## extract those which meet threshold
        dd <- td[tips >= min.size]
    }
    if (rank == "order"){
        ## get all ordinal nodes in tree
        tax <- phy$node.label[grep("[A-z]+ales", phy$node.label, perl=TRUE)]
        ## extract subtrees
        trees <- lapply(tax, function(x) extract.clade(phy, node=x))
        ## use treedata to match to family level
        td <- lapply(trees, function(x) treedata(phy=x, data=data))
        names(td) <- tax
        ## get number of tips
        tips <- lapply(td, function(x) Ntip(x$phy))
        ## extract those which meet threshold
        dd <- td[tips >= min.size]
    }

    dd
}


treedata.time <- function(phy, data, age, min.size) {
    tr <- time.slice.tree(time.slice=age, temp.tree=phy, sr=min.size)

    td <- lapply(tr, function(x) treedata(phy=x, data=data))

    td
}




## function for making data objects with the following elements
## phy
## states
## SE
## trait
## taxa
## rank
build.angio.data.clade <- function(tree.states, rank, trait, min.size){
    td <- treedata.taxon(phy=tree.states$phy, data=tree.states$states,
                         rank=rank, min.size=min.size)

    all.trees <- lapply(seq_len(length(td)), function(x)
                        return(list(phy=td[[x]]$phy, states=td[[x]]$data,
                                    SE=tree.states$SE, taxa=names(td[x]),
                                    rank=rank, trait=trait)))

    all.trees
}


build.angio.data.time <- function(tree.states, age, trait, min.size){
    td <- treedata.time(phy=tree.states$phy, data=tree.states$states,
                        age=age, min.size=min.size)
    
    all.trees <- lapply(seq_len(length(td)), function(x)
                        return(list(phy=td[[x]]$phy, states=td[[x]]$data,
                                    SE=tree.states$SE, taxa="random",
                                    rank="timeslice", trait=trait)))

    all.trees
}



## read in and process data

## time slices
tt <- c(0.2687, 50.2697, 100.2697, 150.2697, 200.2697)

sla <- get.sla.data()
sla.fam <- build.angio.data.clade(sla, rank="family", trait="SLA", min.size=20)
sla.ord <- build.angio.data.clade(sla, rank="order", trait="SLA", min.size=20)
sla.tt1 <- build.angio.data.time(sla, age=tt[1], trait="SLA", min.size=20)
sla.tt2 <- build.angio.data.time(sla, age=tt[2], trait="SLA", min.size=20)
sla.tt3 <- build.angio.data.time(sla, age=tt[3], trait="SLA", min.size=20)
sla.tt4 <- build.angio.data.time(sla, age=tt[4], trait="SLA", min.size=20)
sla.tt5 <- build.angio.data.time(sla, age=tt[5], trait="SLA", min.size=20) 


sdm <- get.seedmass.data()
sdm.fam <- build.angio.data.clade(sdm, rank="family", trait="seedmass", min.size=20)
sdm.ord <- build.angio.data.clade(sdm, rank="order", trait="seedmass", min.size=20)
sdm.tt1 <- build.angio.data.time(sdm, age=tt[1], trait="seedmass", min.size=20)
sdm.tt2 <- build.angio.data.time(sdm, age=tt[2], trait="seedmass", min.size=20)
sdm.tt3 <- build.angio.data.time(sdm, age=tt[3], trait="seedmass", min.size=20)
sdm.tt4 <- build.angio.data.time(sdm, age=tt[4], trait="seedmass", min.size=20)
sdm.tt5 <- build.angio.data.time(sdm, age=tt[5], trait="seedmass", min.size=20) 


lfn <- get.leafn.data()
lfn.fam <- build.angio.data.clade(lfn, rank="family", trait="leafn", min.size=20)
lfn.ord <- build.angio.data.clade(lfn, rank="order", trait="leafn", min.size=20)
lfn.tt1 <- build.angio.data.time(lfn, age=tt[1], trait="leafn", min.size=20)
lfn.tt2 <- build.angio.data.time(lfn, age=tt[2], trait="leafn", min.size=20)
lfn.tt3 <- build.angio.data.time(lfn, age=tt[3], trait="leafn", min.size=20)
lfn.tt4 <- build.angio.data.time(lfn, age=tt[4], trait="leafn", min.size=20)
lfn.tt5 <- build.angio.data.time(lfn, age=tt[5], trait="leafn", min.size=20) 



all.trait.data <- c(sla.fam, sla.ord, sla.tt1, sla.tt2, sla.tt3, sla.tt4, sla.tt5,
                         sdm.fam, sdm.ord, sdm.tt1, sdm.tt2, sdm.tt3, sdm.tt4, sdm.tt5,
                         lfn.fam, lfn.ord, lfn.tt1, lfn.tt2, lfn.tt3, lfn.tt4, lfn.tt5)

## convert the form of the data for downstream use
trait.data.to.vector <- function(x){
    x$states <- as.matrix(x$states)[,1]
    x
}

all.trait.data <- lapply(all.trait.data, function(x) trait.data.to.vector(x))

## write to rds
saveRDS(all.trait.data, file="data/angio-trait-data-all.rds")





