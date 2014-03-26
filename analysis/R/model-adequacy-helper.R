## Assorted functions for helping analysis of adequacy results

## Load in all packages used in model-adequacy.R
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)
library(ape)
library(geiger)

## Functions for making a dataframe with results from only the best model
get.model.suffix <- function(x){
    st <- strsplit(x, split=".", fixed=TRUE)
    sf <- st[[1]][length(st[[1]])]
    sf
}


rm.model.suffix <- function(x){
    st <- strsplit(x, split=".", fixed=TRUE)[[1]]
    if ("ml" %in% st | "mcmc" %in% st)
        st <- st[-(length(st) - 1)]

    strm <- st[-length(st)]
    nofx <- paste(strm, collapse=".")
    nofx      
}


rm.model.suffix.colnames <- function(x){
    colnames(x) <- sapply(colnames(x), function(y) rm.model.suffix(y))
    x
}



prune.dataset.best.ml <- function(x){
    ## colnames to ignore
    ign <- c("taxa", "rank", "trait", "size", "age")
    tmp <- colnames(x)[-which(colnames(x) %in% ign)]
    tmp.df <- x[,tmp]

    ## get model for each column
    suf <- sapply(tmp, function(y) get.model.suffix(y))

    ## find best model
    aic <- c("aic.bm", "aic.ou", "aic.eb")
    mm <- sapply(seq_len(nrow(x)), function(y)
                 return(aic[which(x[y,aic] == min(x[y,aic]))]))
    model <- sapply(mm, function(y) get.model.suffix(y))

    ## build new dataframe inlcuding only the best model
    dd <- lapply(seq_len(nrow(tmp.df)), function(y)
                 return(tmp.df[y,which(suf == model[y])]))

    ## strip column names of last element
    dd <- lapply(dd, function(y) rm.model.suffix.colnames(y))

    df <- do.call(rbind, dd)

    df <- cbind.data.frame(x[,ign], df)
    rownames(df) <- NULL
    df
}



prune.dataset.best.bayes <- function(x){
    ## colnames to ignore
    ign <- c("taxa", "rank", "trait", "size", "age")
    tmp <- colnames(x)[-which(colnames(x) %in% ign)]
    tmp.df <- x[,tmp]

    ## get model for each column
    suf <- sapply(tmp, function(y) get.model.suffix(y))

    ## find best model
    dic <- c("dic.bm", "dic.ou", "dic.eb")
    mm <- sapply(seq_len(nrow(x)), function(y)
                 return(dic[which(x[y,dic] == min(x[y,dic]))]))
    model <- sapply(mm, function(y) get.model.suffix(y))

    ## build new dataframe inlcuding only the best model
    dd <- lapply(seq_len(nrow(tmp.df)), function(y)
                 return(tmp.df[y,which(suf == model[y])]))

    ## strip column names of last element
    dd <- lapply(dd, function(y) rm.model.suffix.colnames(y))

    df <- do.call(rbind, dd)

    df <- cbind.data.frame(x[,ign], df)
    rownames(df) <- NULL
    df
}



## Build table for plotting AIC support versus model adequacy
## Find difference in AIC between best model and BM
## Excluding datasets where BM is best supported model
build.table.adequacy.aic <- function(x){
    ## colnames to ignore
    ign <- c("taxa", "rank", "trait", "size", "age")
    tmp <- colnames(x)[-which(colnames(x) %in% ign)]
    tmp.df <- x[,tmp]

    ## get model for each column
    suf <- sapply(tmp, function(y) get.model.suffix(y))

    ## find best model
    aic <- c("aic.bm", "aic.ou", "aic.eb")
    mm <- sapply(seq_len(nrow(x)), function(y)
                 return(aic[which(x[y,aic] == min(x[y,aic]))]))
    model <- sapply(mm, function(y) get.model.suffix(y))

    ## difference in AIC between best supported model and BM
    diff.bm <- sapply(seq_len(nrow(tmp.df)), function(x)
                      return(tmp.df[x,"aic.bm"] - tmp.df[x,names(model)[x]]))

    ## build new dataframe inlcuding only the best model
    dd <- lapply(seq_len(nrow(tmp.df)), function(y)
                 return(tmp.df[y,which(suf == model[y])]))

    ## strip column names of last element
    dd <- lapply(dd, function(y) rm.model.suffix.colnames(y))

    df <- do.call(rbind, dd)

    df <- cbind.data.frame(x[,ign], df, diff.bm)
    rownames(df) <- NULL

    ## drop cases where bm is the best
    df[-which(df$diff.bm == 0),]
   
}




## Build table for plotting DIC support versus model adequacy
## Find difference in DIC between best model and BM
## Excluding datasets where BM is best supported model
build.table.adequacy.dic <- function(x){
    ## colnames to ignore
    ign <- c("taxa", "rank", "trait", "size", "age")
    tmp <- colnames(x)[-which(colnames(x) %in% ign)]
    tmp.df <- x[,tmp]

    ## get model for each column
    suf <- sapply(tmp, function(y) get.model.suffix(y))

    ## find best model
    dic <- c("dic.bm", "dic.ou", "dic.eb")
    mm <- sapply(seq_len(nrow(x)), function(y)
                 return(dic[which(x[y,dic] == min(x[y,dic]))]))
    model <- sapply(mm, function(y) get.model.suffix(y))

    ## difference in AIC between best supported model and BM
    diff.bm <- sapply(seq_len(nrow(tmp.df)), function(x)
                      return(tmp.df[x,"dic.bm"] - tmp.df[x,names(model)[x]]))

    ## build new dataframe inlcuding only the best model
    dd <- lapply(seq_len(nrow(tmp.df)), function(y)
                 return(tmp.df[y,which(suf == model[y])]))

    ## strip column names of last element
    dd <- lapply(dd, function(y) rm.model.suffix.colnames(y))

    df <- do.call(rbind, dd)

    df <- cbind.data.frame(x[,ign], df, diff.bm)
    rownames(df) <- NULL

    ## drop cases where bm is the best
    df[-which(df$diff.bm == 0),]
   
}










## function for counting the number of p-values below a threshold
count.pvalues <- function(dat, threshold){
    dd <- dat[,c("m.pic", "v.pic", "s.var", "s.anc", "s.hgt", "d.ks")]
    bt <- vector()
    for (i in seq_len(nrow(dd))){
        tmp <- length(which(dd[i,] <= threshold))
        bt <- c(bt,tmp)
    }
    bt
}



## little functions for renaming traits and ranks
rename.traits <- function(x){
    y <- x
    
    if (x == "seedmass")
        y <- "SeedMass"

    if (x == "leafn")
        y <- "LeafN"

    y
}

cap.ranks <- function(x){
    if (x == "family")
        y <- "Family"

    if (x == "timeslice")
        y <- "Timeslice"

    if (x == "order")
        y <- "Order"

    y
}


## Evaluate expression 'expr' that produces a figure as a side effect,
## saving the result in a pdf file.
to.pdf <- function(filename, width, height, expr,
                   ..., pointsize=12, verbose=TRUE) {
  if (verbose)
    cat(sprintf("Creating %s\n", filename))
  pdf(filename, width=width, height=height, pointsize=pointsize, ...)
  on.exit(dev.off())
  res <- eval.parent(substitute(expr))
  # Workaround for ggplot figures:
  if (inherits(res, "gg")) {
    print(res)
  }
}
