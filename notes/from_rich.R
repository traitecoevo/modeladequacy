library(arbutus)
library(geiger)
source("from_rich_rewrite.R")

as.unit.tree <- function(x, ...)
  UseMethod("as.unit.tree")

## Compared with the previous version, we don't check the class of the
## tree, because it will only get here if it is a tree.
as.unit.tree.phylo <- function(x, data, ...) {
  ## check tree and data to make sure they match	
  td <- treedata(phy=x, data=data)
  phy <- td$phy
  data <- td$data
  
  ## calculate pics
  pics <- pic(data, phy, var.contrasts=TRUE)
  
  ## append all the object together
  unit.tree <- list(phy=phy, data=data, pics=pics)
  
  ## change the class of the unit.tree
  class(unit.tree) <- "unit.tree"
  
  unit.tree
}

## Geiger fits are of class c("gfit", "list"), so this will dispatch
## on "gfit".  Because we are using things in geiger that are
## essentially undocumented this could be a bit dangerous, but given
## that we have two of geigers developers (such that that term
## applies) perhaps we can set this up so that we can rely on it a bit
## better.  Think about this before geiger gets resubmitted to CRAN.
as.unit.tree.gfit <- function(x, ...) {
  obj <- parseFitContinuous(x)
  ## NOTE: We can use the .phylo method here now that we have our
  ## ducks in a row.  We could also append any extra information about
  ## the model fit if we wanted, but there doesn't seem to be much
  ## gained at this point.
  as.unit.tree(rescale.fitContinuous(obj), obj$data)
}

## Little utility function -- may come in useful later
is.unit.tree <- function(x)
  inherits(x, "unit.tree")

## So, now we have:
data(geospiza)
tree <- geospiza$phy
data <- geospiza$dat[,1]
fit.bm <- fitContinuous(phy=tree, dat=data, model="BM")

unit.tree <- as.unit.tree(fit.bm)
