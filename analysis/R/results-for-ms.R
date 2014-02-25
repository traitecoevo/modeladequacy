## read in ml data
ml <- read.csv("output/ml-results.csv", row.names=1, as.is=TRUE)


## number of clades examined
nrow(ml)


## pick out best model and create dataframe
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

ml.bst <- prune.dataset.best.ml(ml)


## Two clade example

## meliaceae p values
mel <- ml.bst[which(ml.bst$taxa == "Meliaceae"), ]
mel.pvalue <- mel[c("m.pic", "v.pic", "s.var", "s.anc", "s.hgt", "d.ks")]
mel.pvalue

## fagaceae p values
fg <- ml.bst[which(ml.bst$taxa == "Fagaceae"),]
fg <- fg[which(fg$trait == "seedmass"),]
fg.pvalue <- fg[c("m.pic", "v.pic", "s.var", "s.anc", "s.hgt", "d.ks")]
fg.pvalue




## Empirical results

## AIC support for various models
tmp <- c("aicw.bm", "aicw.ou", "aicw.eb")
dd <- ml[,tmp]
colnames(dd) <- c("BM", "OU", "EB")
dd <- dd[order(dd[,"OU"], dd[,"BM"], decreasing = TRUE),]

best.model <- vector()
for (i in 1:nrow(dd)){
    mx <- colnames(dd)[which(dd[i,] == max(dd[i,]))]
    best.model <- c(best.model, mx)
}

dd <- cbind(dd, best.model)

## how many clades best supported by OU
length(which(dd$best.model == "OU"))

## how many clades OU support == 100%
length(which(dd$OU == 1))

## how many clades OU support > 75%
length(which(dd$OU > 0.75))

## how many clades EB support > 75%
length(which(dd$EB > 0.75))

## include size
dd <- cbind(dd, ml[,"size"])
colnames(dd)[5] <- "size"

## large trees only
dd.lg <- dd[which(dd$size >= 100),]

## how many large trees supported by single model with 90% or more
aic90 <- sapply(seq_len(nrow(dd.lg)), function(x)
                return(max(dd.lg[x,c("BM", "OU", "EB")]) > 0.9))

length(which(aic90))

## how many large trees support OU with >90%
ou.aic90 <- dd.lg[which(aic90),]
length(which(ou.aic90$OU > 0.9))




## pvalues
ml.bst.sla <- ml.bst[which(ml.bst$trait == "SLA"),]
ml.bst.sdm <- ml.bst[which(ml.bst$trait == "seedmass"),]
ml.bst.lfn <- ml.bst[which(ml.bst$trait == "leafn"),]

## fxn for counting p-values below some threshold
count.pvalues <- function(dat, threshold){
    dd <- dat[,c("m.pic", "v.pic", "s.var", "s.anc", "s.hgt", "d.ks")]
    bt <- vector()
    for (i in seq_len(nrow(dd))){
        tmp <- length(which(dd[i,] <= threshold))
        bt <- c(bt,tmp)
    }
    bt
}


## number of clades with SLA data
nrow(ml.bst.sla)
## sla p values
sla.pval <- count.pvalues(ml.bst.sla, 0.05)
## number rejected by at least one
length(which(sla.pval >= 1))
## number rejected by at least two
length(which(sla.pval >= 2))
## number rejected by at least three
length(which(sla.pval >= 3))



## number of clades with seedmass data
nrow(ml.bst.sdm)
## sla p values
sdm.pval <- count.pvalues(ml.bst.sdm, 0.05)
## number rejected by at least one
length(which(sdm.pval >= 1))
## number rejected by at least two
length(which(sdm.pval >= 2))
## number rejected by at least three
length(which(sdm.pval >= 3))



## number of clades with leaf nitrogen data
nrow(ml.bst.lfn)
## sla p values
lfn.pval <- count.pvalues(ml.bst.lfn, 0.05)
## number rejected by at least one
length(which(lfn.pval >= 1))
## number rejected by at least two
length(which(lfn.pval >= 2))
## number rejected by at least three
length(which(lfn.pval >= 3))



## Overall
bst.pval <- count.pvalues(ml.bst, 0.05)
## number which did not violate any of the summary stats
length(which(bst.pval == 0))

## rejected by m.pic
length(which(ml.bst$m.pic <= 0.05))
## rejected by v.pic
length(which(ml.bst$v.pic <= 0.05))
## rejected by s.var
length(which(ml.bst$s.var <= 0.05))
## rejected by s.anc
length(which(ml.bst$s.anc <= 0.05))
## rejected by s.hgt
length(which(ml.bst$s.hgt <= 0.05))
## rejected by d.ks
length(which(ml.bst$d.ks <= 0.05))



## frequency by which obs mean less than simu mean
## get clades where model rejected by m.pic
mpic <- ml.bst[which(ml.bst$m.pic <= 0.05),]
## how many times is the mean less than simulated
length(which(mpic$mean.diag < 500))





## OVERALL DATA NUMBERS FOR METHODS SECTION
## import big tree
tree <- read.tree(file="data/tempo_scrubbed_CONSTRAINT_rooted.dated.tre")

## extract angiosperms
t <- extract.clade(tree, node="Angiospermae")
## number of total angiosperms in tree
Ntip(t)

## age of tree
max(branching.times(t))

## read in SLA data
sla <- read.csv(file="output/species-mean-sla.csv", header=TRUE, row.names=1)

## number of species for which we have sla
nrow(sla)

tmp.sla <-  t$tip.label[!t$tip.label %in% rownames(sla)]

## overlap between tree and sla
Ntip(t) - length(tmp.sla)


## read in seedmass data
sdm <- read.csv(file="output/species-mean-seedMass.csv", row.names=1, header=TRUE)

## number of species for which we have seedmass data
nrow(sdm)

tmp.sdm <- t$tip.label[!t$tip.label %in% rownames(sdm)]

## overlap between tree and seedmass
Ntip(t) - length(tmp.sdm)


## read in leafn data
lfn <- read.csv("output/species-mean-leafN.csv", row.names=1, as.is=TRUE)

## number of species for which we have leafn data
nrow(lfn)

tmp.lfn <- t$tip.label[!t$tip.label %in% rownames(lfn)]

## overlap between tree and leafn
Ntip(t) - length(tmp.lfn)

   
