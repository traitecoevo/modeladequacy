
library(phytools)
#Now look at sub-trees
#old approach, now defunct:
#a<-branching.times(temp.tree) 

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
treedata.taxon <- function(phy, data, rank="family", min.size=20){
    
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

