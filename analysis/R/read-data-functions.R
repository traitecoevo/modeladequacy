## These functions are designed to read and process the data

require(geiger)


## get the data for the SLA
get.sla.data <- function(){
    t <- read.tree(file=file.path(getwd(), "data", "vascular_plant_phylogeny.tre"))
    sla.raw <- read.csv(file=file.path(getwd(), "output", "species-mean-sla.csv"), header=TRUE)

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


get.seedmass.data <- function(){
    t <- read.tree(file=file.path(getwd(), "data", "vascular_plant_phylogeny.tre"))
    sm.raw <- read.csv(file=file.path(getwd(), "output", "species-mean-seedMass.csv"), header=TRUE)

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

    ## drop one of any tips which has same value as sister species
 #   cher <- geiger:::cherries(phy)

  #  chck <- vector()
   # for (i in 1:nrow(cher)){
    #    tmp1 <- sm[cher[i,1]]
     #   tmp2 <- sm[cher[i,2]]
      #  same <- tmp1 == tmp2
       # chck <- c(chck, same)
    #}

    # phy <- geiger:::.drop.tip(phy=phy, tip=names(which(chck)))
    sm <- sm[phy$tip.label]
        
    

    ## return tree and data
    list(phy=phy, states=sm, SE=0.1551108)
}





get.leafn.data <- function(){
    t <- read.tree("data/vascular_plant_phylogeny.tre")
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



get.sla.v.leafn.data <-function(){
    sla <- get.sla.data()
    ln <- get.leafn.data()

    tmp <- sla$phy$tip.label[!(sla$phy$tip.label %in% ln$phy$tip.label)]
    phy <- geiger:::.drop.tip(phy=sla$phy, tip=tmp)

    states.sla <- sla$states[phy$tip.label] 
    states.ln <- ln$states[phy$tip.label]
    states <- cbind.data.frame(states.sla, states.ln)
    colnames(states) <- c("sla", "leafN")

    SE <- c(sla$SE, ln$SE)
    names(SE) <- c("sla", "leafN")

    list(phy=phy, states=states, SE=SE)
}



    
    
    



    
