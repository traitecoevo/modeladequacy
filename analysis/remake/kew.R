## This is the main file that still does things with organising
## downloading only if not there already and things like that.  Ideas
## for fixing this are in my backlog.  But the Kew set is a mess and
## doing this should not be really needed anyway.
##
## Note that that the kew_path() bits are also hard coded here:
## Ideally they'd be set in the remakefile.  But such is life.
kew_build <- function() {
  clades <- kew_clades()
  kew_fetch_clades(clades)
  dat <- lapply(clades, kew_load_clade)
  id <- rep(kew_clades_tr(clades), sapply(dat, nrow))
  combined <- cbind(clade=I(id), do.call(rbind, dat))
  rownames(combined) <- NULL
  combined
}

kew_clades <- function() {
  c("BASAL DICOTS", "MONOCOTS", "COMMELINOIDS", "EUDICOTS",
    "CORE EUDICOTS", "ROSIDS", "EUROSIDS I", "EUROSIDS II",
    "ASTERIDS", "EUASTERIDS I", "EUASTERIDS II", "UNCERTAIN")
}

kew_clades_tr <- function(x) {
  tolower(sub(" ", "_", x))
}

kew_load_clade <- function(clade, verbose=TRUE) {
  filename <- kew_csv_name(clade)
  if (!file.exists(filename)) {
    message("Reading ", clade)
    dat <- kew_load_html(kew_html_name(clade))
    write.csv(dat, filename, row.names=FALSE)
  }
  read.csv(filename, stringsAsFactors=FALSE)  
}

kew_url <- function(clade=NULL) {
  family <- ""
  url_fmt <-
    "http://data.kew.org/sid/SidServlet?Clade=%s&Order=&Family=%s&APG=off&Genus=&Species=&StorBehav=0&WtFlag=on"
  sub(" ", "%20", sprintf(url_fmt, clade, family))
}

kew_path <- function() {
  "data/kew"
}
kew_html_name <- function(clade) {
  file.path(kew_path(),sprintf("kew_apg_%s.html", kew_clades_tr(clade)))
}
kew_csv_name <- function(clade) {
  file.path(kew_path(), sprintf("kew_apg_%s.csv", kew_clades_tr(clade)))
}

## Simple throttling:
kew_fetch_clades <- function(clades, timeout=20, verbose=TRUE) {
  for (clade in clades) {
    did_download <- kew_fetch_clade(clade, verbose)
    if (did_download && timeout > 0) {
      if (verbose) {
        message("...resting")
      }
      Sys.sleep(timeout)
    }
  }
}

kew_fetch_clade <- function(clade, verbose=TRUE) {
  filename <- kew_html_name(clade)
  already_exists <- file.exists(filename)
  if (already_exists) {
  } else {
    if (verbose) {
      message(sprintf("Fetching %s", clade))
    }
    dir.create(dirname(filename), FALSE, TRUE)
    download.file(kew_url(clade), filename, quiet=!verbose)
  }
  invisible(!already_exists)
}

## Below here is parsing code.  It uses XML and is slow.
kew_load_html <- function(filename) {
  html <- htmlParse(filename)

  ## Get the data container:
  container <- xpathSApply(html, "//*[@id='sid']")[[1]]

  ## Next bit is done manually.
  kids <- xmlChildren(container)
  if (sum(names(kids) == "b") != 1)
    stop("Unexpected input (expected one 'b'))")
  i <- match("b", names(kids))
  header_1 <- container[[i]]   # number of records
  header_2 <- container[[i+1]] # key to style
  content  <- container[[i+2]]

  ## Parse the header
  info <- kew_parse_header_key(header_2)
  key <- names(info)[info == "Mean 1000 Seed Weight"]
  ret <- kew_parse_content(content, key)
  ret <- do.call(rbind, ret)
  data.frame(species=unname(ret[,"species"]),
             value=as.numeric(sub("g$", "", ret[,"value"])),
             stringsAsFactors=FALSE)
}

kew_parse_header_key <- function(header) {
  tmp <- xmlChildren(header)
  if (!all(names(tmp) == "span"))
    stop("Unknown type in header")
  attrs <- unname(sapply(tmp, function(x)
                         xmlAttrs(x)[["style"]]))
  txt <- sub(", *", "", unname(sapply(tmp, xmlValue)))
  names(txt) <- attrs
  txt
}

kew_parse_content <- function(content, info) {
  kids <- xmlChildren(content)
  idx <- unname(which(names(kids) == "a"))
  grp <- rep(seq_along(idx), diff(c(idx, length(kids)+1)))
  lapply(split(kids, grp), kew_process_species, info)
}

kew_process_species <- function(x, key) {
  vals <- x[names(x) == "span"]
  type <- unname(sapply(vals, function(x)
                        xmlAttrs(x)[["style"]]))
  value <- strstrip(xmlValue(vals[[which(type == key)]]))
  c(species=xmlValue(x[["a"]]), value=value)
}
