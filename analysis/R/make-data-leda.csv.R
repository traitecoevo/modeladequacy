#!/usr/bin/env Rscript

## An embedded ';' in three references causes some trouble with
## importing.  This substitutes the ';' for a ':'.
tmp <- readLines("data/leda.txt")[-(1:4)]
## This is how I identified the problem (the perl=TRUE is required
## because of some non-ASCII characters).
##   tmp.split <- strsplit(tmp, ";", perl=TRUE)
##   n <- sapply(tmp.split, length)
##   i <- which(n != 23)
##   tmp.split[i]
## Replace some troublesome lines:
tmp2 <- sub("daylight; lowered", "daylight: lowered", tmp,
            perl=TRUE)
f <- tempfile()
writeLines(tmp2, f)
d <- read.table(f, sep=";",
                stringsAsFactors=FALSE, check.names=FALSE,
                comment.char="", header=TRUE)

tr <- c("SBS name"="SBS.name",
        "general method"="GeneralMethod",
        "single value [mm^2/mg]"="SingleValue",
        "sample size"="SampleSize",
        "valid"="Valid",
        "leaf specific method"="LeafSpecificMethod",
        "reference"="Reference",
        "SBS number"="SBS.number",
        "original reference"="OriginalReference",
        "country"="Country",
        "UTM zone"="UTM.zone",
        "UTM easting"="UTM.easting",
        "UTM northing"="UTM.northing",
        "mean SLA [mm^2/mg]"="SLA.mean",
        "maximum SLA [mm^2/mg]"="SLA.max",
        "minimum SLA [mm^2/mg]"="SLA.min",
        "number of replicates"="Replicates",
        "standard deviation"="StandardDeviation",
        "standard error"="StandardError",
        "balance error [mg]"="BalanceError",
        "collection date"="CollectionDate",
        "general comment"="Comment",
        "plant stage"="PlantStage")
names(d)[match(names(tr), names(d))] <- tr
write.csv(d, "data/leda.csv", row.names=FALSE)
