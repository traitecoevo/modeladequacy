## wrapper functions


phy.model.check <- function(x, nsim=1000, stats=NULL, ...){

    ## create unit tree
    unit.tree <- as.unit.tree(x, ...)

    ## compute summary statistics on observed data
    ss.obs <- summ.stats(unit.tree, stats=stats)

    ## simulate n datasets
    sims <- sim.char.unit(unit.tree, nsim=nsim)

    ## compute summary statistics on simulated data
    ss.sim <- summ.stats(sims, stats=stats)

    ## compare summary statistics between observed and simulated
    p.values <- compare.summ.stats(ss.obs, ss.sim)

    p.values
}
