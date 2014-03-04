# Common phylogenetic trait models are widely inadequate: Analysis

If you have time on your hands and `make` is installed, you only have to run `make` and everything will run.

Required packages (in addition to arbutus)

* LaplacesDemon (currently not on CRAN?)
* XML
* ape
* cwhmisc
* digest
* diversitree
* geiger
* ggplot2
* phytools
* scales
* xlsx

# Contents:

## The `data` directory.

Most of these files are downloaded from elsewhere, but a few are local:

* `errors.csv`: Data value error fixing (see `R/load-scrubbing-tools.R`)
* `names-tr.txt`: Species name translations to fix spelling errors, etc (see `R/load-scrubbing-tools.R`)
* `spermatophyta_synonyms_PLANTLIST.csv` list of synonyms derived from the Plant List.  This will change to use the woodiness data shortly.

Data downloaded from other sources::

* `kew.csv`, `kew/`: Seed weight data from kew (CITATION)
* `leda.csv` (and `leda.txt`): SLA data (CITATION)
* `wright-2004.csv` (and `wright-2004.xls`): SLA data (CITATION)

---

I'd (RGF) suggest that the final directory structure could look like
this:

* analysis/ -- directory with all the analysis mess in it
  - data/   -- downloaded data sets
  - R/      -- scripts unrelated to the package
  -         -- plus other assorted files to make it run
* pkg/      -- the package, but directory correctly named
  - R/      -- } all the usual
  - man/    -- } package stuff
* ms/       -- manuscript, as currently set up

So, as we start accumulating things that are associated primarily with
the analyis (such as downloading data, cleaning up, etc), they could
go in here.

I'm structuring this with a Makefile; typing "make" should create
everything so far.  That will get more interesting as the analysis
gets larger.

I'll list interesting make targets as we get them.  So far we have:

* `all`: (the default target) -- do everything
* `data-fetch`: Fetches all remote data sets.  Currently just
  `data/wright-2004.xls`.
* `data-preprocess`: Turn all fetched data sets into useful csv data
  sets.

The convention I'm using for file naming is that a script called
`R/make-foo.R` will make a file called something like `foo`.  There
will probably end up being some additional scripts that are used by
several of these files.
