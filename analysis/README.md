# Model adequacy and the macroevolution of angiosperm functional traits: Analysis

If you have time on your hands and `make` is installed, you only have to run `make` and everything will run.

In slightly more detail:

1. Download the required data sets and do some initial data cleaning:

```
make data-preprocess
```

Or, to grab a cached data set (avoiding hammering Kew)

```
make downloaded-data-bulk-fetch
make downloaded-data-unpack
make data-preprocess
```


2. Clean up data sets taxonomically and generate simple data sets by species (in the `output` directory).

```
make data
```

Currently:
* `leafN-process.R`: `data/wright-2004.csv` ->  `output/species-mean-leafN.csv`
* `seedMass-process.R`: `data/kew.csv` -> `output/species-mean-seedMass.csv`
* `sla-process.R:trait`: {`data/wright-2004`, `data/leda.csv`} -> `output/species-mean-sla.csv`

(that list out of date)

3. Run the analyses

```
make fits
```

This runs the code in `model-adequacy-ml.R` and `model-adequacy-bayes.R`, which goes through and fits models using ML or MCMC (respectively) and assesses adequacy using arbutus' `phy.model.check` function.

Individual fits are stored in `data/results-ml` and `data/results-bayes` and then summarised together in `data/results-ml.csv` and `data/results-bayes.csv`.

4. Process the analysis

```
make analysis
```

This will run all the code in `model-analysis.R`, creating figures (in `output/figs`).  It also creates the file `model-adequacy.html` which shows the process of running all the code.  This makes use of knitr.

**Note**: For now, the dependencies for the fits and analysis section aren't done automatically (to save time during development).  So to run everything in one fell swoop:

```
make downloaded-data-bulk-fetch downloaded-data-unpack
make data fits analysis
```

# Other information:

Required packages

* XML
* arbutus
* ape
* dplyr (>= 0.1.2)
* digest
* diversitree
* geiger
* ggplot2
* xlsx

Running the analysis.  Because of the length of time this takes, there are two steps.  First, fit models and run the model adequacy assessment on each data set (there are 337 of these, and each takes a few minutes or more).  This is done in the files `model-adequacy-ml.R` and `model-adequacy-bayes.R`.

# Contents:

## The `data` directory.

Most of these files are downloaded from elsewhere, but a few are local:

* `errors.csv`: Data value error fixing (see `R/data-process-taxonomic.R`)
* `names-tr.txt`: Species name translations to fix spelling errors, etc (see `R/data-process-taxonomic.R`)
* `spermatophyta_synonyms_PLANTLIST.csv` list of synonyms derived from the Plant List.  This will change to use the woodiness data shortly.

Data downloaded from other sources::

* `kew.csv`, `kew/`: Seed weight data from Kew; http://data.kew.org/sid
* `leda.csv` (and `leda.txt`): SLA data from LEDA traitbase; Kleyer et al. 2008, Journal of Ecology 96:1266-1274.
* `wright-2004.csv` (and `wright-2004.xls`): SLA and leaf nitrogen data from Wright et al. 2004, Nature 428:821-827.

## The `R` directory.

Functions for handling the data, which are all sourced by the R files in `make/`:

* `data-process-angio.R`
* `data-process-kew.R`
* `data-process-taxonomic.R`

Functions that are called by `analysis/model-adequacy-ml.R` and `analysis/model-adequacy-bayes.R` to assist with the model fitting and analysis of model adequacy:

* `model-adequacy-fit.R`

Functions for assisting with plotting and the analysis of final results:

* `model-adequacy-analysis.R`

The file `paths.R` contains paths and helpers to load partly processed data.

## Main analysis scripts

There are two primary scripts for fitting models and assessing model adequacy::

* `model-ad-angio-ml.R`
* `model-ad-angio-bayes.R`

These two scripts fit the 3 evolutionary models (using ml or mcmc, respectively) and assess model adequacy. Both of these use multiple processes and the analyses will take some time.

All of the analysis of the results produced from the above functions are handled with:

* `model-adequacy.R`

This file contains scripts which produce all the figures and numbers that appear in the manuscript (with the exception of those that are not done in R).






