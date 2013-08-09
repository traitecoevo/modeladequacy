kew.url <- function(family=NULL, clade=NULL) {
  if (is.null(family) && is.null(clade))
    stop("One of family or clade must be given")
  if (!is.null(family) && !is.null(clade))
    stop("Only one of family or clade must be given")
  if (is.null(family))
    family <- ""
  if (is.null(clade))
    clade <- ""
  url.fmt <-
    "http://data.kew.org/sid/SidServlet?Clade=%s&Order=&Family=%s&APG=off&Genus=&Species=&StorBehav=0&WtFlag=on"
  sub(" ", "%20", sprintf(url.fmt, clade, family))
}

kew.fetch <- function(filename, family=NULL, clade=NULL) {
  already.existed <- file.exists(filename)
  if (already.existed) {
    message(sprintf("Skipping %s -- already exists", filename))
  } else {
    dat <- readLines(kew.url(family, clade))
    writeLines(dat, filename)
  }
  invisible(!already.existed)
}

kew.clades <- function()
  c("BASAL DICOTS", "MONOCOTS", "COMMELINOIDS", "EUDICOTS",
    "CORE EUDICOTS", "ROSIDS", "EUROSIDS I", "EUROSIDS II",
    "ASTERIDS", "EUASTERIDS I", "EUASTERIDS II", "UNCERTAIN")
kew.clades.tr <- function(x)
  tolower(sub(" ", "_", x))

kew.fetch.clades <- function(timeout=60) {
  for (clade in kew.clades()) {
    message(sprintf("Fetching %s", clade))
    dest <- sprintf("data/kew_apg_%s.html", kew.clades.tr(clade))
    did.download <- kew.fetch(dest, clade=clade)
    if (did.download && timeout > 0) {
      message("...resting")
      Sys.sleep(timeout)
    }
  }
}

kew.html.to.csv <- function(filename) {
  write.csv(kew.load.html(filename),
            paste0(tools::file_path_sans_ext(filename), ".csv"),
            row.names=FALSE)
}

kew.load.html <- function(filename) {
  require(XML)
  html <- htmlParse(filename)

  ## Get the data container:
  container <- xpathSApply(html, "//*[@id='sid']")[[1]]

  ## Next bit is done manually.
  kids <- xmlChildren(container)
  if (sum(names(kids) == "b") != 1)
    stop("Unexpected input (expected one 'b'))")
  i <- match("b", names(kids))
  header.1 <- container[[i]]   # number of records
  header.2 <- container[[i+1]] # key to style
  content  <- container[[i+2]]

  ## Parse the header
  info <- kew.parse.header.key(header.2)
  key <- names(info)[info == "Mean 1000 Seed Weight"]
  ret <- kew.parse.content(content, key)
  ret <- do.call(rbind, ret)
  data.frame(species=unname(ret[,"species"]),
             value=as.numeric(sub("g$", "", ret[,"value"])),
             stringsAsFactors=FALSE)
}

kew.parse.header.key <- function(header) {
  tmp <- xmlChildren(header)
  if (!all(names(tmp) == "span"))
    stop("Unknown type in header")
  attrs <- unname(sapply(tmp, function(x)
                         xmlAttrs(x)[["style"]]))
  txt <- sub(", *", "", unname(sapply(tmp, xmlValue)))
  names(txt) <- attrs
  txt
}

kew.parse.content <- function(content, info) {
  kids <- xmlChildren(content)
  idx <- unname(which(names(kids) == "a"))
  grp <- rep(seq_along(idx), diff(c(idx, length(kids)+1)))
  lapply(split(kids, grp), kew.process.species, info)
}

kew.process.species <- function(x, key) {
  vals <- x[names(x) == "span"]
  type <- unname(sapply(vals, function(x)
                        xmlAttrs(x)[["style"]]))
  value <- strip(xmlValue(vals[[which(type == key)]]))
  c(species=xmlValue(x[["a"]]), value=value)
}

strip <- function(x)
  gsub("(^ +| +$)", "", x)
