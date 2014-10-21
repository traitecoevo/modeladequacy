# Model adequacy and the macroevolution of angiosperm functional traits: Analysis

We developed this repository and structured our analysis with the goal of making our methods transparent and our results completely reproducible.

## Dependencies:

First, install [maker](https://github.com/richfitz/maker) from github using devtools (maker is not yet on CRAN as it is under active development).

```r
devtools::install_github("richfitz/maker")
```

We depend on quite a few packages for the analysis.  The full list is:

  - XML
  - ape
  - arbutus
  - diversitree
  - dplyr
  - geiger
  - ggplot2
  - grid
  - gridExtra
  - parallel
  - reshape2
  - xslx
  - sowsear

All packages except [arbutus](https://github.com/mwpennell/arbutus) and [sowsear](https://github.com/richfitz/sowsear) are currently avaialbe on CRAN.

To install `arbutus` and `sowsear`, use devtools:

```r
devtools::install_github("mwpennell/arbutus")
devtools::install_github("richfitz/sowsear")
```

The remaining packages can then by installed via CRAN.

Eventually you should be able to run:

```r
maker::make("deps")
```

to install all dependencies, but this currently only works for CRAN dependencies and will skip xlsx.

## Running the analysis

From within the `analysis/` directory run `maker::make()`.  This will take several hours, depending on your computer.

Or, in several steps:

1. Construct a "maker" object that we'll interact with:

```r
m <- maker$new()
```

2. Download all the data and build into data sets

```r
m$make("data")
```

3. Split the data into a large number of subsets and run the analysis on each (both ML and Bayesian)

```r
m$make("fits")
```

This takes a while!  If you have a powerful computer, you will find that setting the `mc.cores` option helps (e.g. `options(mc.cores=8)` to use 8 cores).

4. Generate the report and figures

*(still under development as we move to the new system)*


**NOTE** all the instructions below are out of date as we transition to using [maker](https://github.com/richfitz/maker) to control the workflow.  They may contain useful information though.

---

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
make downloaded-data-cache-fetch
make downloaded-data-cache-unpack
make data-preprocess
```


### Clean up data sets taxonomically and generate simple data sets by species (in the `output` directory).

```
make data
```

On a fresh install this will run the following R scripts:

* `make/data-vascular_plant_phylogeny.tre.R`
* `make/output-vascular_plant_phylogeny.rds.R`
* `make/output-synonyms.rds.R`
* `make/output-corrections.rds.R`
* `make/data-wright-2004.csv.R`
* `make/output-species-mean-leafN.csv.R`
* `make/output-data-leafN.rds.R`
* `make/data-leda.csv.R`
* `make/output-species-mean-sla.csv.R`
* `make/output-data-sla.rds.R`
* `make/data-kew.csv.R`
* `make/output-species-mean-seedMass.csv.R`
* `make/output-data-seedMass.rds.R`
* `make/output-angio-data.R`

These all build proceessed data objects (for example `output-species-mean-seedMass.csv.R` builds the file `output/species-mean-seedMass.csv`).

Running all these scripts should take a couple of minutes (probably 3-10 minutes).

### Run the analyses

```
make fits
```

This runs the code in `model-adequacy-ml.R` and `model-adequacy-bayes.R`, which goes through and fits models using ML or MCMC (respectively) and assesses adequacy using arbutus' `phy.model.check` function. If multiple processes are available, this can be sped up by changing the line `options(mc.cores=2)` in `model-adequacy-ml.R` and/or `model-adequacy-bayes.R`.

This can be run in parallel: by default we use two cores, but this can be modified by changing the line:

```
options(mc.cores=2)
```

in both files.

Individual fits are stored in `output/results-ml` and `output/results-bayes` directories.  If the fitting process is interrupted, rerunning `make fits` will skip over fits that have completed.  Note that this will generate about 1GB of files.

It takes about 5 CPU hours to run the ML analyses, and 13 CPU hours for the Bayesian analysis (so about 9 hours total over two processors).

All the fits are then summarised together in the csv files `output/results-ml.csv` (for ML) and `output/results-bayes.csv` (for Bayesian).

### Process the analysis

```
make analysis
```

This will run all the code in `model-analysis.R`, creating figures (in `output/figs`).  It also creates the file `model-adequacy.html` which shows the process of running all the code.  This makes use of knitr. The files summarizing the final results for the paper `output/results-ml.csv` and `output/results-bayes.csv` are included in the repository so `make analysis` can be used to generate the output files without re-running the analyses. This may be useful for researchers who want to know exactly where the numbers in the paper come from without having to rerun the analyses from scratch (which will take several hours). Also note that owing to the stochastic nature of the procedure, re-running the analyses will result in slightly different numbers from those reported in the paper.

**Note**: For now, the dependencies for the fits and analysis section aren't loaded automatically (to save time during development).  So to run everything in one fell swoop:

```
make downloaded-data-cache-fetch downloaded-data-cache-unpack
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






