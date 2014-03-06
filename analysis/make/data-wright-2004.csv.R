#!/usr/bin/env Rscript
library(methods)
library(xlsx, quietly=TRUE)
library(digest, quietly=TRUE)

file.csv <- "data/wright-2004.csv"

## There are severl strategies for reading in an excel file, but this
## one works quite well.
d <- read.xlsx2("data/wright-2004.xls", sheetIndex=1, startRow=11,
                stringsAsFactors=FALSE, check.names=FALSE)

## A change here may cause a false positive warning.  I think changes
## in xlsx, among other things, can cause minor changes to the file
## and that causes the hash to change.
hash.obj.sha <- "8817c8f1adf28d41e7b138238e2ac11b0f540d2a"
if (digest(d, algo="sha1") != hash.obj.sha)
  warning("Imported data did not match expected hash")

## Do some name translations:
tr <- c("Code"="Code",
        "Dataset"="Dataset",
        "BIOME"="Biome",
        "Species"="Species",
        "GF"="GrowthForm",
        "Decid/E'green"="Deciduous",
        "Needle/Broad lf"="Needle",
        "C3C4"="C3",
        "N2-fixer"="N2fixer",
        "log LL"="LogLeafLifespan",
        "log LMA"="LogLMA",
        "log Nmass"="Log.N.mass",
        "log Narea"="Log.N.area",
        "log Pmass"="Log.P.mass",
        "log Parea"="Log.P.area",
        "log Amass"="Log.A.mass",
        "log Aarea"="Log.A.area",
        "log Gs"="Log.Gs",
        "log Rdmass"="Log.Rd.mass",
        "log Rdarea"="Log.Rd.area",
        "Ca - Ci"="CaCi")
names(d)[match(names(tr), names(d))] <- tr

## Be aware of these issues with species names -- probably we should
## clean these out at some point.
## grep(" .+ ", d$Species, value=TRUE)

category.to.logical <- function(x, trueval) {
  x[x == ""] <- NA
  x == trueval
}

## Drop blank columns
d <- d[names(d) != " "]

## Data tweaking:
d[["Code"]] <- as.integer(d[["Code"]])
d[["CaCi"]] <- as.numeric(d[["CaCi"]])

d[["Deciduous"]] <- category.to.logical(d[["Deciduous"]], "D")
d[["Needle"]]    <- category.to.logical(d[["Needle"]],    "N")
d[["C3"]]        <- category.to.logical(d[["C3"]],        "C3")
d[["N2fixer"]]   <- category.to.logical(d[["N2fixer"]],   "Y")

## The logged values are log10(), not log(), so create non-logged
## versions to remove doubt.  My (RGF) preference would be to drop the
## logged versions entirely and just log transform when needed.
re <- "^Log\\."
i.log <- grep(re, names(d))
d[i.log] <- lapply(d[i.log], as.numeric)
d.unlogged <- as.data.frame(10^d[i.log])
names(d.unlogged) <- sub(re, "", names(d.unlogged))
d <- cbind(d, d.unlogged)

write.csv(d, file.csv, row.names=FALSE)
