#!/usr/bin/env Rscript

# Quite possible that this would be better to do with packrat:
#   http://rstudio.github.io/packrat/
# but this should do for now.

pkgs.cran <- c(ape="3.0-11",
               dplyr="0.1.2.0.99",
               digest="0.6.4",
               diversitree="0.9-7",
               geiger="2.0",
               ggplot2="0.9.3.1",
               xlsx="0.5.5",
               XML="3.98-1.1",
               grid="3.1.0",
               gridExtra="0.9.1",
               reshape2="1.2.2",
               knitr="1.5")
pkgs.github <- c("mwpennell/arbutus"="1.0",
                 "richfitz/sowsear"="0.1-1")

# There's a bunch of effort here to avoid downloading things
# needlessly.  If a simpler version would simply do:
#
#   install.packages(c(names(pkgs.cran), "devtools"))
#   library(devtools)
#   for (pkg in names(pkgs.github)) {
#     install_github(pkg)
#   }
#
# which would be less fragile, but involves downloading more things
# and modifying the users computer more than ideal.

pkgs <- installed.packages()

# CRAN:
to.install <- setdiff(names(pkgs.cran), rownames(pkgs))

installed <- intersect(names(pkgs.cran), rownames(pkgs))
to.upgrade <-
  installed[numeric_version(pkgs[installed, "Version"]) <
            numeric_version(pkgs.cran[installed])]

msg.install <- (if (length(to.install) == 0) character(0) else
                paste("installing:", paste(to.install, collapse=", ")))
msg.upgrade <- (if (length(to.upgrade) == 0) character(0) else
                paste("upgrading:", paste(to.upgrade, collapse=", ")))
msg <- c(msg.install, msg.upgrade)
if (length(msg) > 0) {
  message(paste(msg, collapse="\n"))
  install.packages(c(to.install, to.upgrade))
}

# GitHub:
pkgs.github.short <- sub("^.+/", "", names(pkgs.github))

tr <- structure(names(pkgs.github), names=pkgs.github.short)

to.install <- setdiff(pkgs.github.short, rownames(pkgs))
installed <- intersect(pkgs.github.short, rownames(pkgs))

to.upgrade <-
  installed[numeric_version(pkgs[installed, "Version"]) <
            numeric_version(pkgs.github[tr[installed]])]

msg.install <- (if (length(to.install) == 0) character(0) else
                paste("installing:", paste(to.install, collapse=", ")))
msg.upgrade <- (if (length(to.upgrade) == 0) character(0) else
                paste("upgrading:", paste(to.upgrade, collapse=", ")))
msg <- c(msg.install, msg.upgrade)
if (length(msg) > 0) {
  message(paste(msg, collapse="\n"))
  install.packages("devtools")
  library(devtools)
  for (pkg in tr[c(to.install, to.upgrade)]) {
    install_github(pkg)
  }
}
