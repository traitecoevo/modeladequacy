#!/usr/bin/env Rscript
library(dockertest)
build_remake("clean")
build_remake("all")
