all: data-preprocess

data-fetch: data/wright-2004.xls data/leda.txt kew-fetch zae-fetch
data-preprocess: data/wright-2004.csv data/leda.csv data/kew.csv \
	data/genus_order_lookup.csv data/vascular_plant_phylogeny.tre

# Little helpers:
kew-fetch: data/kew_apg_basal_dicots.html
zae-fetch: data/zae/PhylogeneticResources.zip data/zae/Spermatophyta_Genera.csv

# Data from Wright et al.:
data/wright-2004.xls:
	Rscript R/make-data-wright-2004.xls.R
data/wright-2004.csv: data/wright-2004.xls R/make-data-wright-2004.csv.R
	Rscript R/make-data-wright-2004.csv.R

# Data from the LEDA database:
data/leda.txt:
	Rscript R/make-data-leda.txt.R
data/leda.csv: data/leda.txt R/make-data-leda.csv.R
	Rscript R/make-data-leda.csv.R

# Data from Kew:
data/kew_apg_basal_dicots.html:
	Rscript R/make-data-kew-fetch.R
data/kew.csv: data/kew_apg_basal_dicots.html
	Rscript R/make-data-kew.csv.R

# Data from Zanne et al.:
data/zae/Spermatophyta_Genera.csv:
	Rscript R/make-data-zae-Spermatophyta_Genera.csv.R
data/zae/PhylogeneticResources.zip:
	Rscript R/make-data-zae-PhylogeneticResources.zip.R
data/genus_order_lookup.csv: data/zae/Spermatophyta_Genera.csv
	Rscript R/make-data-genus_order_lookup.csv.R
data/vascular_plant_phylogeny.tre: data/zae/PhylogeneticResources.zip
	Rscript R/make-data-vascular_plant_phylogeny.tre.R

.PHONY: all data-fetch data-preprocess kew-fetch zae-fetch
