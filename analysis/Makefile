all: data/wright-2004.csv

data-fetch: data/wright-2004.xls
data-preprocess: data/wright-2004.csv

data/wright-2004.xls:
	Rscript R/make-data-wright-2004.xls.R
data/wright-2004.csv: data/wright-2004.xls R/make-data-wright-2004.csv.R
	Rscript R/make-data-wright-2004.csv.R

.PHONY: all data-fetch data-preprocess
