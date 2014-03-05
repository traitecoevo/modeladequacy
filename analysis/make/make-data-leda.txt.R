#!/usr/bin/env Rscript
url <-
  "http://www.leda-traitbase.org/LEDAportal/objects/Data_files/SLA%20und%20geo.txt"
file <- "data/leda.txt"
hash.file.md5 <- "6f46c9f2beab30cfc45f7d0186400f69"

download.file(url, file)

if (tools::md5sum(file) != hash.file.md5)
  warning("Downloaded file did not match expected hash",
          immediate.=TRUE)
