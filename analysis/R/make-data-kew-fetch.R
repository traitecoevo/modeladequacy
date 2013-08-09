#!/usr/bin/env Rscript

## This is a proof of concept.  Probably need to be a bit careful in
## running this because you don't want to overwhelm the server and get
## in trouble!

source("R/kew.R")
kew.fetch.clades()
