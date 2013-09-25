## wrapper functions


phyModelCheck <- function(x, nsim=1000, stats=NULL, ...){

    ## create unit tree
    unit.tree <- as.unit.tree(x, ...)

    ## compute summary statistics on observed data
    ss.obs <- summStat(unit.tree, stats=stats)

    ## simulate n datasets
    sims <- sim.charUnit(unit.tree, nsim=nsim)

    ## compute summary statistics on simulated data
    ss.sim <- summStat(sims, stats=stats)

    ## compare summary statistics between observed and simulated
    p.values <- compare.summStat(ss.obs, ss.sim)

    p.values
}
