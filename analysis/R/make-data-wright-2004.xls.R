#!/usr/bin/env Rscript
url <-
  "http://www.nature.com/nature/journal/v428/n6985/extref/nature02403-s2.xls"
file.xls <- "data/wright-2004.xls"
hash.file.md5 <- "5f0a94c5dda1ae5c6a39ae8c4e2e681c"
download.file(url, file.xls)
if (tools::md5sum(file.xls) != hash.file.md5)
  warning("Downloaded xls file did not match expected hash",
          immediate.=TRUE)
