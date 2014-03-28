## functions for building rds file for all data and all traits
require(geiger, quietly=TRUE)

source("R/paths.R")

build.data <- function(dataset) {
  dataset <- match.arg(dataset, c("sla", "seedMass", "leafN"))

  t <- get.tree()
  raw <- read.csv(sprintf("output/species-mean-%s.csv", dataset),
                  stringsAsFactors=FALSE)
  t <- extract.clade(t, node="Angiospermae")

  ## NOTE: base 10 log
  dat <- structure(log10(raw$mean), names=raw$gs)

  ## Drop species from tree not in data, and v.v.
  ## This step is the slow point.
  phy <- geiger:::.drop.tip(phy=t, tip=setdiff(t$tip.label, names(dat)))
  dat <- dat[phy$tip.label]

  ## TODO: recompute SLA from data?
  ##   mean(log10(raw$sd[log10(raw$sd) > 0.0001]), na.rm=TRUE)
  se <- c(sla=0.1039405, seedMass=0.1551108, leafN=0.07626127)[[dataset]]

  list(phy=phy, states=dat, SE=se)
}

treedata.q <- function(...)
  geiger::treedata(..., warnings=FALSE)

## functions for extracting subtrees
## some of this code was written by Jon M. Eastman
get.node.heights <- function(t){
  edg <- as.matrix(arbutus:::edge.height(t))
  edg <- max(edg[,"start"]) - edg
  dimnames(edg) <- NULL
  # This is the height of the start and end of each node (so the first
  # row is tip 1, row Ntip(tree)+1 is the root, etc, so put it into
  # edge matrix order:
  edg[t$edge[,2],]
}

extract.all.nodes<-function(tree,node.height.min,node.height.max){
  #note that this returns a matrix of the same dimensions as tree$edge 
  #containing the height above the root of each node in edge
  tree$node.label<-rep("",length(tree$node.label))
  a<-get.node.heights(tree)
  interesting.nodes<-tree$edge[which(a[,2]>node.height.min&a[,2]<node.height.max)] #these are edges which pass through a certain age
  non.terminal.interesting.nodes<-interesting.nodes[interesting.nodes>length(tree$node.label)+2]
  return(non.terminal.interesting.nodes)
}

where.to.cut<-function(tree, age) {
  tree$node.label <- rep("", tree$Nnode)
  a <- get.node.heights(tree)
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

# TODO(RGF): Why is this 'temp.tree'?  Clean this up.
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
        td <- lapply(trees, function(x) treedata.q(phy=x, data=data))
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
        td <- lapply(trees, function(x) treedata.q(phy=x, data=data))
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

    td <- lapply(tr, function(x) treedata.q(phy=x, data=data))

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


logspace.f <- function(x, f, ...) {
  exp(f(log(x), ...))
}

geometric.mean <- function(x, ...) {
  logspace.f(x, mean, ...)
}

geometric.sd <- function(x, ...) {
  logspace.f(x, sd, ...)
}
