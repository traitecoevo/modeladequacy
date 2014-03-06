## Assorted functions for helping analysis of adequacy results

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
                 return(aic[which(x[y,dic] == min(x[y,dic]))]))
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
