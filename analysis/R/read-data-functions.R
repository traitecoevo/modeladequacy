## These functions are designed to read and process the data

require(geiger)

## get the data for the SLA
get.sla.data <- function(tree.file, species.mean.file){
    t <- read.tree(file=tree.file)
    sla.raw <- read.csv(file=species.mean.file, header=TRUE)

    ## natural log
    sla <- log10(sla.raw[,"x"])
    names(sla) <- sla.raw[,"X"]

    ## drop extra tips
    tmp <- t$tip.label[!t$tip.label %in% names(sla)]
    phy <- geiger:::.drop.tip(phy=t, tip = tmp)

    ## return tree and data
    list(phy=phy, states=sla, SE=0.1024202)
}
