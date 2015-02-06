# Model adequacy and the macroevolution of angiosperm functional traits: Analysis

We developed this repository and structured our analysis with the goal of making our methods transparent and our results completely reproducible.  Be warned though that doing this takes some time.

## Dependencies:

First, install [remake](https://github.com/richfitz/remake) from github using devtools (remake is not yet on CRAN as it is under active development).

```r
devtools::install_github("richfitz/remake")
```

We depend on quite a few packages for the analysis.  To install any missing packages, run:

```r
remake::install_missing_packages()
```

This step also requires `devtools` to install [arbutus](https://github.com/mwpennell/arbutus) and [sowsear](https://github.com/richfitz/sowsear) from github (see [remake_sources.yml](remake_sources.yml))

The full list of required packages is in [the `packages:` section of remake.yml](analysis/remake.yml), plus the `xlsx` package to build `wright_2004`.

## Running the analysis

From within the `analysis/` directory run `remake::make()`.  This will take several hours, depending on your computer, and print out some information to what it is doing.

Or, in several steps:


1. Download all the data and build into data sets

```r
remake::make("data")
```

Note that one of the data sources (Kew) is continuously updated without a versioning control system and therefore the downloaded data may be different from what is presented in the paper.  We will organise a cached verison of this data.

3. Split the data into a large number of subsets and run the analysis on each (both ML and Bayesian)

```r
remake::make("fits")
```

This step takes a while (about 5 CPU hours to run the ML analyses, and 13 CPU hours for the Bayesian analysis: so about 9 hours total over two processors).
If you have a powerful computer, you will find that setting the global `mc.cores` option before running this helps (e.g. `options(mc.cores=8)` to use 8 cores).  It will also produce close to 1GB of output, which would be a problem if you are very tight for disk space.

4. Generate the report and figures

```
remake::make()
```

This will build an html report (`model-adequacy.html`) that includes all numbers referenced in the paper, plus all the figures (pdf versions in `output/figs`, plus summary tables (`output/results-ml.csv` and `output/results-bayes.csv`).  Note that owing to the stochastic nature of the procedure, re-running the analyses will result in slightly different numbers from those reported in the paper.

To explore further, run:
```
e <- remake:::make_dependencies(remake::remake(), "model-adequacy.md")
remake::remake_attach(e)
```

and then any of the code in `model-adequacy.R` will be runnable.

## The `data` directory.

Most of the data used in this analysis are downloaded from elsewhere, but a few are local:

* `errors.csv`: Data value error fixing
* `names-tr.txt`: Species name translations to fix spelling errors, etc
* `spermatophyta_synonyms_PLANTLIST.csv` list of synonyms derived from the Plant List, by Jon Eastman.

Data downloaded from other sources:

* `kew.csv`, `kew/`: Seed weight data from Kew; http://data.kew.org/sid
* `leda.csv` (and `leda.txt`): SLA data from LEDA traitbase; Kleyer et al. 2008, Journal of Ecology 96:1266-1274.
* `wright_2004.csv` (and `wright_2004.xls`): SLA and leaf nitrogen data from Wright et al. 2004, Nature 428:821-827.
