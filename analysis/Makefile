all: data

output/angio-data: make/output-angio-data.R ${DATA_PHY}
	Rscript $<

# NOTE: the targets here are potentially screwy, and things may not
# always get rebuilt at the appropriate times.  Because of the manual
# caching that goes on, we might be best having these be phony
# targets.
output/results-ml.csv: model-adequacy-ml.R output/angio-data
	Rscript $<
output/results-bayes.csv: model-adequacy-bayes.R output/angio-data
	Rscript $<

# Here's the main analysis:
model-adequacy.Rmd: model-adequacy.R
	Rscript -e "library(sowsear); sowsear('$<', 'Rmd')"
model-adequacy.md: model-adequacy.Rmd
	Rscript -e "library(knitr); knit('$<')"
model-adequacy.html: model-adequacy.md
	Rscript -e "library(markdown);\
	 opts <- setdiff(markdownHTMLOptions(TRUE), 'base64_images');\
	 markdownToHTML('$<', '$@', options=opts)"
