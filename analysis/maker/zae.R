zae_prefix <- function() {
  'http://datadryad.org/bitstream/handle/10255/'
}

download_zae_trees <- function(destination_filename) {
  url <- paste0(zae_prefix(),
                'dryad.55548/PhylogeneticResources.zip?sequence=1')
  download.file(url, destination_filename)
}

download_zae_genera <- function(destination_filename) {
  url <- paste0(zae_prefix(),
                'dryad.55304/Spermatophyta_Genera.csv?sequence=2')
  download.file(url, destination_filename)
}

unpack_tree <- function(resources_zip) {
  tmp <- tempdir()
  ## unzip() will not throw an error if a the file does not exist in
  ## the archive, so this is here to force an error to happen.  This
  ## occurs if the zip file is corrupted for example.
  oo <- options(warn=2)
  on.exit(options(oo))
  unzip(resources_zip,
        'PhylogeneticResources/Vascular_Plants_rooted.dated.tre',
        junkpaths=TRUE, exdir=tmp)
  read.tree(file.path(tmp, "Vascular_Plants_rooted.dated.tre"))
}

unpack_genus_order_lookup <- function(filename) {
  lookup <- read.csv(filename, stringsAsFactors=FALSE)
  names(lookup)[1] <- "genus"
  lookup <- lookup[c("genus", "family", "order")]
  ## Two genera need families set:
  lookup$family[lookup$genus == "Peltanthera"] <- "Gesneriaceae"
  lookup$family[lookup$genus == "Brachynema" ] <- "Olacaceae"
  lookup
}
