all: data-preprocess

data-fetch: data/wright-2004.xls data/leda.txt kew-fetch zae-fetch
data-preprocess: data/wright-2004.csv data/leda.csv data/kew.csv \
	data/genus_order_lookup.csv data/vascular_plant_phylogeny.tre
data-species-mean: output/species-mean-leafN.csv \
	output/species-mean-sla.csv output/species-mean-seedMass.csv
data-with-trees: output/data-leafN.rds output/data-sla.rds \
	output/data-seedMass.rds

# Little helpers:
kew-fetch: data/kew/kew_apg_uncertain.html
zae-fetch: data/zae/PhylogeneticResources.zip data/zae/Spermatophyta_Genera.csv

# Data from Wright et al.:
data/wright-2004.xls:
	Rscript make/make-data-wright-2004.xls.R
data/wright-2004.csv: data/wright-2004.xls make/make-data-wright-2004.csv.R
	Rscript make/make-data-wright-2004.csv.R

# Data from the LEDA database:
data/leda.txt:
	Rscript make/make-data-leda.txt.R
data/leda.csv: data/leda.txt make/make-data-leda.csv.R
	Rscript make/make-data-leda.csv.R

# Data from Kew:
data/kew/kew_apg_uncertain.html:
	Rscript make/make-data-kew-fetch.R
data/kew.csv: data/kew/kew_apg_uncertain.html
	Rscript make/make-data-kew.csv.R

# Data from Zanne et al.:
data/zae/Spermatophyta_Genera.csv:
	Rscript make/make-data-zae-Spermatophyta_Genera.csv.R
data/zae/PhylogeneticResources.zip:
	Rscript make/make-data-zae-PhylogeneticResources.zip.R
data/genus_order_lookup.csv: data/zae/Spermatophyta_Genera.csv
	Rscript make/make-data-genus_order_lookup.csv.R
data/vascular_plant_phylogeny.tre: data/zae/PhylogeneticResources.zip
	Rscript make/make-data-vascular_plant_phylogeny.tre.R

# Data processing
output/vascular_plant_phylogeny.rds: data/vascular_plant_phylogeny.tre
	Rscript make/make-output-vascular_plant_phylogeny.rds.R
## We'll refer to this a bunch, so shortcut:
TREE = output/vascular_plant_phylogeny.rds

output/species-mean-leafN.csv: ${TREE} data/wright-2004.csv make/make-output-species-mean-leafN.csv.R
	Rscript make/make-output-species-mean-leafN.csv.R
output/species-mean-sla.csv: ${TREE} data/wright-2004.csv data/leda.csv make/make-output-species-mean-sla.csv.R
	Rscript make/make-output-species-mean-sla.csv.R
output/species-mean-seedMass.csv: ${TREE} data/kew.csv make/make-output-species-mean-seedMass.csv.R
	Rscript make/make-output-species-mean-seedMass.csv.R

output/data-leafN.rds: output/species-mean-leafN.csv R/read-data-functions.R
	Rscript make/make-output-data-leafN.rds.R
output/data-sla.rds: output/species-mean-sla.csv R/read-data-functions.R
	Rscript make/make-output-data-sla.rds.R
output/data-seedMass.rds: output/species-mean-seedMass.csv R/read-data-functions.R
	Rscript make/make-output-data-seedMass.rds.R

.PHONY: all data-fetch data-preprocess data-species-mean kew-fetch zae-fetch
