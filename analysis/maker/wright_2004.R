download_wright_2004 <- function(destination_filename) {
  url <-
    "http://www.nature.com/nature/journal/v428/n6985/extref/nature02403-s2.xls"
  download.file(url, destination_filename)
}

process_wright_2004 <- function(filename) {
  ## There are several strategies for reading in an excel file, but
  ## this one works quite well.
  ##
  ## I'm keeping the xlsx dependency here rather than within the
  ## makerfile because it loads really slowly.
  library(methods) # Serious, WTF Rscript?
  library(xlsx, quietly=TRUE) # Yay, using xlsx::read won't work...
  d <- read.xlsx2(filename, sheetIndex=1, startRow=11,
                  stringsAsFactors=FALSE, check.names=FALSE)

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

  ## Drop blank columns
  d <- d[names(d) != " "]

  ## Data tweaking:
  d[["Code"]] <- as.integer(d[["Code"]])
  d[["CaCi"]] <- as.numeric(d[["CaCi"]])

  d[["Deciduous"]] <- category_to_logical(d[["Deciduous"]], "D")
  d[["Needle"]]    <- category_to_logical(d[["Needle"]],    "N")
  d[["C3"]]        <- category_to_logical(d[["C3"]],        "C3")
  d[["N2fixer"]]   <- category_to_logical(d[["N2fixer"]],   "Y")

  ## The logged values are log10(), not log(), so create non-logged
  ## versions to remove doubt.  My (RGF) preference would be to drop the
  ## logged versions entirely and just log transform when needed.
  re <- "^Log\\."
  i_log <- grep(re, names(d))
  d[i_log] <- lapply(d[i_log], as.numeric)
  d_unlogged <- as.data.frame(10^d[i_log])
  names(d_unlogged) <- sub(re, "", names(d_unlogged))
  d <- cbind(d, d_unlogged)

  d
}
