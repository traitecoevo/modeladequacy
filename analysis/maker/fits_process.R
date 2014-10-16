combine_fits <- function(fits) {
  do.call(rbind, lapply(fits, process))
}

process <- function(res) {
  type <- res$BM$info$type
  col_ic <- if (type == "ml") "aic" else "dic"

  process1 <- function(x) {
    model <- tolower(x$info$model)
    # P values for the various statistics:
    pv <- pvalue_arbutus(x$ma)
    names(pv) <- sprintf("%s.%s.%s", names(pv), type, model)

    # Mean distance
    #
    # NOTE: This fails sometimes (computationally singular) but I
    # think I'm out of date with arbutus.  For now, just setting this
    # to use try(), and setting to NaN on error so we can track it
    # down later.
    mv <- try(mean(mahalanobis_arbutus(x$ma)), silent=TRUE)
    if (inherits(mv, "try-error")) {
      mv <- NaN
    }
    names(mv) <- sprintf("mv.%s.%s", type, model)

    # Mean diagnostic
    #
    # This is not strictly necessary but was used to test whether
    # m.pic was consistently under or overestimated This is the only
    # summary statistic which we expected to be systematically biased
    # in one direction
    md <- sum(x$ma$sim$m.sig < x$ma$obs$m.sig)
    names(md) <- sprintf("mean.diag.%s", model)

    as.list(c(pv, mv, md))
  }
  
  info <- res$BM$info[c("taxa", "rank", "trait", "size", "age")]

  ic <- sapply(res, function(x) x$info[[col_ic]])
  icw <- ic_weights(ic)
  names(ic)  <- sprintf("%s.%s",  col_ic, tolower(names(ic)))
  names(icw) <- sprintf("%sw.%s", col_ic, tolower(names(icw)))

  stats <- unlist_with_names(lapply(res, process1))
  
  as.data.frame(c(info, ic, icw, stats), stringsAsFactors=FALSE)
}

## Compute weights for AIC or DIC
## takes a named vector of AIC/DIC values
ic_weights <- function(x){
  d_x <- x - min(x)
  tmp <- exp(-0.5 * d_x)
  tmp / sum(tmp)
}

## Compute deviance information criterion from mcmcsamples
dic_mcmcsamples <- function(x, burnin=0){
  if (!inherits(x, "mcmcsamples"))
    stop("this function is only designed for diversitree's mcmcsamples object")

  lik <- attr(x, "func")
  p <- coef(x, burnin=burnin)
  dev <- -2 * apply(p, 1, lik)
  
  ## estimate effective number of parameters
  dbar <- mean(dev)

  ## deviance of posterior means:
  post_means <- colMeans(p)
  ## evaluate deviance at the mean posterior estimate
  dhat <- -2 * lik(post_means)

  pd <- dbar - dhat

  ## calculate dic
  dic <- dbar + pd
  unname(dic)
}
