all: data

data-fetch: data/wright-2004.xls data/leda.txt kew-fetch zae-fetch
data-preprocess: data/wright-2004.csv data/leda.csv data/kew.csv \
	data/genus_order_lookup.csv data/vascular_plant_phylogeny.tre
data-species-mean: output/species-mean-leafN.csv \
	output/species-mean-sla.csv output/species-mean-seedMass.csv

DATA_PHY = output/data-leafN.rds output/data-sla.rds output/data-seedMass.rds
data-with-trees: ${DATA_PHY}
data: output/angio-data

run: output/results-ml.csv output/results-mcmc.csv

# Little helpers:
kew-fetch: data/kew/kew_apg_uncertain.html
zae-fetch: data/zae/PhylogeneticResources.zip data/zae/Spermatophyta_Genera.csv

# Data from Wright et al.:
data/wright-2004.xls:
	Rscript make/data-wright-2004.xls.R
data/wright-2004.csv: data/wright-2004.xls make/data-wright-2004.csv.R
	Rscript make/data-wright-2004.csv.R

# Data from the LEDA database:
data/leda.txt:
	Rscript make/data-leda.txt.R
data/leda.csv: data/leda.txt make/data-leda.csv.R
	Rscript make/data-leda.csv.R

# Data from Kew:
data/kew/kew_apg_uncertain.html:
	Rscript make/data-kew-fetch.R
data/kew.csv: data/kew/kew_apg_uncertain.html
	Rscript make/data-kew.csv.R

# Data from Zanne et al.:
data/zae/Spermatophyta_Genera.csv:
	Rscript make/data-zae-Spermatophyta_Genera.csv.R
data/zae/PhylogeneticResources.zip:
	Rscript make/data-zae-PhylogeneticResources.zip.R
data/genus_order_lookup.csv: data/zae/Spermatophyta_Genera.csv
	Rscript make/data-genus_order_lookup.csv.R
data/vascular_plant_phylogeny.tre: data/zae/PhylogeneticResources.zip
	Rscript make/data-vascular_plant_phylogeny.tre.R

# Data processing
output/vascular_plant_phylogeny.rds: data/vascular_plant_phylogeny.tre
	Rscript make/output-vascular_plant_phylogeny.rds.R
## We'll refer to this a bunch, so shortcut:
TREE = output/vascular_plant_phylogeny.rds

output/synonyms.rds: ${TREE} data/spermatophyta_synonyms_PLANTLIST.csv
	Rscript make/output-synonyms.rds.R
output/corrections.rds: ${TREE} data/names-tr.txt
	Rscript make/output-corrections.rds.R
SYN  = output/synonyms.rds output/corrections.rds

output/species-mean-leafN.csv: ${TREE} ${SYN} data/wright-2004.csv make/output-species-mean-leafN.csv.R
	Rscript make/output-species-mean-leafN.csv.R
output/species-mean-sla.csv: ${TREE} ${SYN} data/wright-2004.csv data/leda.csv make/output-species-mean-sla.csv.R
	Rscript make/output-species-mean-sla.csv.R
output/species-mean-seedMass.csv: ${TREE} ${SYN} data/kew.csv make/output-species-mean-seedMass.csv.R
	Rscript make/output-species-mean-seedMass.csv.R

output/data-leafN.rds: output/species-mean-leafN.csv R/read-data-functions.R
	Rscript make/output-data-leafN.rds.R
output/data-sla.rds: output/species-mean-sla.csv R/read-data-functions.R
	Rscript make/output-data-sla.rds.R
output/data-seedMass.rds: output/species-mean-seedMass.csv R/read-data-functions.R
	Rscript make/output-data-seedMass.rds.R

output/angio-data: make/output-angio-data.R ${DATA_PHY}
	Rscript $<

output/results-ml.csv: model-adequacy-ml.R output/angio-data
	Rscript $<

output/results-bayes.csv: model-adequacy-bayes.R output/angio-data
	Rscript $<

# Save on some farting about with data:
DOWNLOADED_DATA = 		\
	data/kew/*.html 	\
	data/leda.txt		\
	data/wright-2004.xls 	\
	data/zae/PhylogeneticResources.zip data/zae/Spermatophyta_Genera.csv
DOWNLOADED_DATA_SAVE = .downloaded_data.tar.gz

downloaded-data-delete:
	rm -f ${DOWNLOADED_DATA}
downloaded-data-save:
	tar -zcvf .downloaded_data.tar.gz ${DOWNLOADED_DATA}
downloaded-data-unpack:
	tar -zxvf .downloaded_data.tar.gz
downloaded-data-bulk-fetch:
	curl -o .downloaded_data.tar.gz http://www.zoology.ubc.ca/~fitzjohn/files/modeladequacy_data.tar.gz

clean-data: clean-processed
	rm -rf data/kew data/kew.csv
	rm -f data/leda.txt data/leda.csv
	rm -f data/vascular_plant_phylogeny.tre
	rm -f data/genus_order_lookup.csv
	rm -f data/wright-2004.xls data/wright-2004.csv
	rm -rf data/zae

## More to add here, perhaps
clean-processed:
	rm -rf output/angio-data
	rm -f output/data-*rds
	rm -f output/species-mean-*csv
	rm -f ${TREE}
	rm -f ${SYN}

clean-outputs:
	rm -rf output/figs
	rm -rf output/results-ml
	rm -rf output/results-bayes

.PHONY: all data data-fetch data-preprocess data-species-mean \
	kew-fetch zae-fetch \
	downloaded-data-delete downloaded-data-save \
	downloaded-data-unpack downloaded-data-bulk-fetch \
	clean-data clean-processed
