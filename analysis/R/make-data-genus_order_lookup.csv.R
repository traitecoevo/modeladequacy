#!/usr/bin/env Rscript

lookup <- read.csv(file.path("data/zae/Spermatophyta_Genera.csv"),
                   stringsAsFactors=FALSE)
names(lookup)[1] <- "genus"
lookup <- lookup[c("genus", "family", "order")]

# Two genera need families set:
lookup$family[lookup$genus == "Peltanthera"] <- "Gesneriaceae"
lookup$family[lookup$genus == "Brachynema" ] <- "Olacaceae"

write.csv(lookup, file.path("data/genus_order_lookup.csv"),
          row.names=FALSE)
