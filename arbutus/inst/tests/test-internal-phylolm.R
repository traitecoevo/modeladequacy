source("helper-arbutus.R")

context("phylolm (internals)")

model.type <- arbutus:::model.type
model.data <- arbutus:::model.data
model.pars <- arbutus:::model.pars
model.info <- arbutus:::model.info

model.type.phylolm <- arbutus:::model.type.phylolm
model.data.phylolm <- arbutus:::model.data.phylolm
model.pars.phylolm <- arbutus:::model.pars.phylolm
model.info.phylolm <- arbutus:::model.info.phylolm

set.seed(1)
phy <- tree.bd(pars=c(1,0), max.taxa=100)
tx <- sim.character(phy, 1)
ty <- sim.character(phy, 1) + 3
data <- data.frame(x=tx, y=ty, row.names=names(tx))

## fit all of the available models
fit.lm.bm <- phylolm(y~x, data=data, phy=phy, model="BM")
fit.lm.or <- phylolm(y~x, data=data, phy=phy, model="OUrandomRoot")
fit.lm.of <- phylolm(y~x, data=data, phy=phy, model="OUfixedRoot")
fit.lm.la <- phylolm(y~x, data=data, phy=phy, model="lambda")
fit.lm.ka <- phylolm(y~x, data=data, phy=phy, model="kappa")
fit.lm.de <- phylolm(y~x, data=data, phy=phy, model="delta")
fit.lm.eb <- phylolm(y~x, data=data, phy=phy, model="EB")
## ignore trend for now
## fit.lm.tr <- phylolm(y~x, data=data, phy=phy, model="trend")

test_that("Model types are correct", {
    expect_that(model.type(fit.lm.bm), is_identical_to("BM"))
    expect_that(model.type(fit.lm.or), is_identical_to("OU"))
    expect_that(model.type(fit.lm.of), is_identical_to("OU"))
    expect_that(model.type(fit.lm.la), is_identical_to("lambda"))
    expect_that(model.type(fit.lm.ka), is_identical_to("kappa"))
    expect_that(model.type(fit.lm.de), is_identical_to("delta"))
    expect_that(model.type(fit.lm.eb), is_identical_to("EB"))
})

## note: correct data here are the residuals
test_that("Models return their source data", {
    drop.attributes <- function(x){
        attr(x, "node.state") <- NULL
        y <- as.numeric(x)
        names(y) <- rownames(x)
        y
    }

    cmp.bm <- list(phy=phy, data=drop.attributes(resid(fit.lm.bm)))
    cmp.or <- list(phy=phy, data=drop.attributes(resid(fit.lm.or)))
    cmp.of <- list(phy=phy, data=drop.attributes(resid(fit.lm.of)))
    cmp.la <- list(phy=phy, data=drop.attributes(resid(fit.lm.la)))
    cmp.ka <- list(phy=phy, data=drop.attributes(resid(fit.lm.ka)))
    cmp.de <- list(phy=phy, data=drop.attributes(resid(fit.lm.de)))
    cmp.eb <- list(phy=phy, data=drop.attributes(resid(fit.lm.eb)))
                   
    ## processing phylolm objects requires that phylogeny be provided
    ## not included in fitted object
    ## check to make sure this throws error
    expect_that(model.data(fit.lm.bm), throws_error())
    expect_that(model.data(fit.lm.or), throws_error())
    expect_that(model.data(fit.lm.of), throws_error())
    expect_that(model.data(fit.lm.la), throws_error())
    expect_that(model.data(fit.lm.ka), throws_error())
    expect_that(model.data(fit.lm.de), throws_error())
    expect_that(model.data(fit.lm.eb), throws_error())

    expect_that(model.data(fit.lm.bm, phy=phy), equals(cmp.bm))
    expect_that(model.data(fit.lm.or, phy=phy), equals(cmp.or))
    expect_that(model.data(fit.lm.of, phy=phy), equals(cmp.of))
    expect_that(model.data(fit.lm.la, phy=phy), equals(cmp.la))
    expect_that(model.data(fit.lm.ka, phy=phy), equals(cmp.ka))
    expect_that(model.data(fit.lm.de, phy=phy), equals(cmp.de))
    expect_that(model.data(fit.lm.eb, phy=phy), equals(cmp.eb))
})


## check that parameters are correct
## phylolm actuall returns sigsq so no need to extract it
test_that("Parameters correctly returned", {
    ## check sigsq
    s2.bm <- fit.lm.bm$sigma2
    s2.or <- fit.lm.or$sigma2
    s2.of <- fit.lm.of$sigma2
    s2.la <- fit.lm.la$sigma2
    s2.ka <- fit.lm.ka$sigma2
    s2.de <- fit.lm.de$sigma2
    s2.eb <- fit.lm.eb$sigma2

    expect_that(model.pars(fit.lm.bm)$sigsq, equals(s2.bm))
    expect_that(model.pars(fit.lm.or)$sigsq, equals(s2.or))
    expect_that(model.pars(fit.lm.of)$sigsq, equals(s2.of))
    expect_that(model.pars(fit.lm.la)$sigsq, equals(s2.la))
    expect_that(model.pars(fit.lm.ka)$sigsq, equals(s2.ka))
    expect_that(model.pars(fit.lm.de)$sigsq, equals(s2.de))
    expect_that(model.pars(fit.lm.eb)$sigsq, equals(s2.eb))

    ralpha <- fit.lm.or$optpar
    falpha <- fit.lm.of$optpar
    lambda <- fit.lm.la$optpar
    kappa  <- fit.lm.ka$optpar
    delta  <- fit.lm.de$optpar
    a      <- fit.lm.eb$optpar
    
    expect_that(model.pars(fit.lm.or)$alpha, equals(ralpha))
    expect_that(model.pars(fit.lm.of)$alpha, equals(falpha))
    expect_that(model.pars(fit.lm.la)$lambda, equals(lambda))
    expect_that(model.pars(fit.lm.ka)$kappa, equals(kappa))
    expect_that(model.pars(fit.lm.de)$delta, equals(delta))
    expect_that(model.pars(fit.lm.eb)$a, equals(a))

})


test_that("Coefficient names are as expected", {
    pars.bm <- c("sigsq", "z0", "SE")
    pars.ou <- c("sigsq", "alpha", "z0", "SE")
    pars.la <- c("sigsq", "lambda","z0", "SE")
    pars.ka <- c("sigsq", "kappa",  "z0", "SE")
    pars.de <- c("sigsq", "delta", "z0", "SE")
    pars.eb <- c("sigsq", "a",  "z0", "SE")

    expect_that(names(model.pars(fit.lm.bm)), is_identical_to(pars.bm))
    expect_that(names(model.pars(fit.lm.or)), is_identical_to(pars.ou))
    expect_that(names(model.pars(fit.lm.of)), is_identical_to(pars.ou))
    expect_that(names(model.pars(fit.lm.la)), is_identical_to(pars.la))
    expect_that(names(model.pars(fit.lm.ka)), is_identical_to(pars.ka))
    expect_that(names(model.pars(fit.lm.de)), is_identical_to(pars.de))
    expect_that(names(model.pars(fit.lm.eb)), is_identical_to(pars.eb))
})

test_that("Overall object looks legit",{
    obj.bm <- model.info(fit.lm.bm, phy=phy)
    obj.or <- model.info(fit.lm.or, phy=phy)
    obj.of <- model.info(fit.lm.of, phy=phy)
    obj.la <- model.info(fit.lm.la, phy=phy)
    obj.ka <- model.info(fit.lm.ka, phy=phy)
    obj.de <- model.info(fit.lm.de, phy=phy)
    obj.eb <- model.info(fit.lm.eb, phy=phy)

    obj.names <- c("data", "pars", "type")
    expect_that(names(obj.bm), is_identical_to(obj.names))
    expect_that(names(obj.or), is_identical_to(obj.names))
    expect_that(names(obj.of), is_identical_to(obj.names))
    expect_that(names(obj.la), is_identical_to(obj.names))
    expect_that(names(obj.ka), is_identical_to(obj.names))
    expect_that(names(obj.de), is_identical_to(obj.names))
    expect_that(names(obj.eb), is_identical_to(obj.names))

    expect_that(class(obj.bm), is_identical_to("fitC"))
    expect_that(class(obj.or), is_identical_to("fitC"))
    expect_that(class(obj.of), is_identical_to("fitC"))
    expect_that(class(obj.la), is_identical_to("fitC"))
    expect_that(class(obj.ka), is_identical_to("fitC"))
    expect_that(class(obj.de), is_identical_to("fitC"))
    expect_that(class(obj.eb), is_identical_to("fitC"))
})
    
