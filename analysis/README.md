# Model adequacy and the macroevolution of angiosperm functional traits: Analysis

We developed this repository and structured our analysis with the goal of making our methods transparent and our results completely reproducible. 

To repeat our analyses, first ensure that [make](https://www.gnu.org/software/make) is installed on your computer. The analyses also rely on the following R packages:

* arbutus
* ape
* dplyr (>= 0.1.2)
* digest
* diversitree
* geiger
* ggplot2
* xlsx
* XML
* grid
* gridExtra
* reshape2
* knitr
* sowsear

All packages except [arbutus](https://github.com/mwpennell/arbutus) and [sowsear](https://github.com/richfitz/sowsear) are currently avaialbe on CRAN. To install `arbutus` and `sowsear`, the easiest way to install is to use [devtools](https://github.com/hadley/devtools). Install `devtools` then type
```
library(devtools)
install_github("mwpennell/arbutus")
install_github("richfitz/sowsear")
```

Once everything is installed, clone the repository
```
git clone https://github.com/richfitz/modeladequacy.git
```
move into the `/analysis` folder and type

```
make
```
This will download all of the data, process the data, run all analyses and output a document `model-adequacy.html` presenting all results used in the paper and the figures for the main text and supplementary materials (in `output/figs`). Note that this will take several hours to run (see below for how to utilize parallel processing to speed up the analyses). 



The steps done by `make` can be run separately:

### Download the required data sets and do some initial data cleaning:

```
make data-preprocess
```

Note that one of the data sources (Kew) is continuously updated without a versioning control system and therefore the downloaded data may be different from what is presented in the paper. To use the same dataset, pull down a cached data set. This also has the added benefit of not hammering Kew's servers.

```
make downloaded-data-bulk-fetch
make downloaded-data-unpack
make data-preprocess
```


### Clean up data sets taxonomically and generate simple data sets by species (in the `output` directory).

```
make data
```

Currently:
* `leafN-process.R`: `data/wright-2004.csv` ->  `output/species-mean-leafN.csv`
* `seedMass-process.R`: `data/kew.csv` -> `output/species-mean-seedMass.csv`
* `sla-process.R:trait`: {`data/wright-2004`, `data/leda.csv`} -> `output/species-mean-sla.csv`

(that list out of date)


### Run the analyses

```
make fits
```

This runs the code in `model-adequacy-ml.R` and `model-adequacy-bayes.R`, which goes through and fits models using ML or MCMC (respectively) and assesses adequacy using arbutus' `phy.model.check` function. If multiple processes are available, this can be sped up by changing the line `options(mc.cores=2)` in `model-adequacy-ml.R` and/or `model-adequacy-bayes.R` 

Individual fits are stored in `data/results-ml` and `data/results-bayes` and then summarised together in `data/results-ml.csv` and `data/results-bayes.csv`.


### Process the analysis

```
make analysis
```

This will run all the code in `model-analysis.R`, creating figures (in `output/figs`).  It also creates the file `model-adequacy.html` which shows the process of running all the code.  This makes use of knitr.

**Note**: For now, the dependencies for the fits and analysis section aren't loaded automatically (to save time during development).  So to run everything in one fell swoop:

```
make downloaded-data-bulk-fetch downloaded-data-unpack
make data fits analysis
```




## Overview of contents:

### The `data` directory.

Most of these files are downloaded from elsewhere, but a few are local:

* `errors.csv`: Data value error fixing (see `R/data-process-taxonomic.R`)
* `names-tr.txt`: Species name translations to fix spelling errors, etc (see `R/data-process-taxonomic.R`)
* `spermatophyta_synonyms_PLANTLIST.csv` list of synonyms derived from the Plant List.  This will change to use the woodiness data shortly.

Data downloaded from other sources::

* `kew.csv`, `kew/`: Seed weight data from Kew; http://data.kew.org/sid
* `leda.csv` (and `leda.txt`): SLA data from LEDA traitbase; Kleyer et al. 2008, Journal of Ecology 96:1266-1274.
* `wright-2004.csv` (and `wright-2004.xls`): SLA and leaf nitrogen data from Wright et al. 2004, Nature 428:821-827.

### The `R` directory.

Functions for handling the data, which are all sourced by the R files in `make/`:

* `data-process-angio.R`
* `data-process-kew.R`
* `data-process-taxonomic.R`

Functions that are called by `analysis/model-adequacy-ml.R` and `analysis/model-adequacy-bayes.R` to assist with the model fitting and analysis of model adequacy:

* `model-adequacy-fit.R`

Functions for assisting with plotting and the analysis of final results:

* `model-adequacy-analysis.R`

The file `paths.R` contains paths and helpers to load partly processed data.


### Main analysis scripts

There are two primary scripts for fitting models and assessing model adequacy::

* `model-adequacy-ml.R`
* `model-adequacy-bayes.R`

These two scripts fit the 3 evolutionary models (using ml or mcmc, respectively) and assess model adequacy. Both of these use multiple processes and the analyses will take some time.

All of the analysis of the results produced from the above functions are handled with:

* `model-adequacy.R`

This file contains scripts which produce all the figures and numbers that appear in the manuscript (with the exception of those that are not done in R).






