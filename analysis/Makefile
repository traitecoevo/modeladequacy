all: data-preprocess

data-fetch: data/wright-2004.xls data/leda.txt
data-preprocess: data/wright-2004.csv data/leda.csv

data/wright-2004.xls:
	Rscript R/make-data-wright-2004.xls.R
data/wright-2004.csv: data/wright-2004.xls R/make-data-wright-2004.csv.R
	Rscript R/make-data-wright-2004.csv.R

data/leda.txt:
	Rscript R/make-data-leda.txt.R
data/leda.csv: data/leda.txt R/make-data-leda.csv.R
	Rscript R/make-data-leda.csv.R

.PHONY: all data-fetch data-preprocess
