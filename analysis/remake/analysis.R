guess_ic <- function(fits) {
  if ("aic.bm" %in% names(fits)) "aic" else "dic"
}

prune_dataset_best <- function(fits) {
  ## colnames to ignore
  ign <- c("taxa", "rank", "trait", "size", "age")
  keep <- setdiff(names(fits), ign)
  tmp_df <- fits[, keep]

  ## get model for each column
  suf <- get_model_suffix(names(tmp_df))

  ## find best model
  ic <- guess_ic(fits)
  ic_cols <- sprintf("%s.%s", ic, c("bm", "ou", "eb"))
  mm <- ic_cols[apply(fits[ic_cols], 1, which.min)]
  model <- get_model_suffix(mm)

  ## build new dataframe inlcuding only the best model, and strip the
  ## model suffix from this.
  dd <- lapply(seq_len(nrow(tmp_df)), function(i)
               rm_model_suffix_colnames(tmp_df[i, suf == model[i]]))
  df <- do.call(rbind, dd)

  df <- cbind.data.frame(fits[, ign], df)
  rownames(df) <- NULL
  df
}


## Functions for making a dataframe with results from only the best model
get_model_suffix <- function(x) {
  sub(".*\\.([^.]+)$", "\\1", x)
}

rm_model_suffix <- function(x) {
  sub("\\.(ml\\.|bayes\\.)?([^.]+)$", "", x, perl=TRUE)
}

rm_model_suffix_colnames <- function(x){
  colnames(x) <- rm_model_suffix(colnames(x))
  x
}

pvalue_names_arbutus <- function() {
  # "m.pic", "v.pic", "s.var", "s.anc", "s.hgt", "d.ks"  
  c("m.sig", "c.var", "s.var", "s.asr", "s.hgt", "d.cdf")
}

count_pvalues <- function(dat, alpha) {
  dd <- as.matrix(dat[, pvalue_names_arbutus()])
  rowSums(dd <= alpha)
}

## Functions for making a dataframe with results from only the best model
get_model_suffix <- function(x) {
  sub(".*\\.([^.]+)$", "\\1", x)
}

rm_model_suffix <- function(x) {
  sub("\\.(ml\\.|bayes\\.)?([^.]+)$", "", x, perl=TRUE)
}

rm_model_suffix_colnames <- function(x){
  colnames(x) <- rm_model_suffix(colnames(x))
  x
}

build_ic <- function(fits) {
  models <- c("bm", "ou", "eb")
  type <- if ("aic.bm" %in% names(fits)) "aic" else "dic"
  col_names <- paste0(type, "w.", models)
  ic <- fits[col_names]
  colnames(ic) <- toupper(models)
  ic <- ic[order(ic[,"OU"], ic[,"BM"], decreasing = TRUE),]
  ic
}

## Build table for plotting AIC support versus model adequacy
## Find difference in AIC between best model and BM
## Excluding datasets where BM is best supported model
build_table_adequacy_ic <- function(fits) {
  fits_best <- prune_dataset_best(fits)

  ic <- guess_ic(fits)
  ic_bm <- paste0(ic, ".bm")

  ## difference in AIC between best supported model and BM
  diff_bm <- fits[[ic_bm]] - fits_best[[ic]]

  ret <- cbind(fits_best, diff.bm=diff_bm)
  rownames(ret) <- NULL
  ret[ret$diff.bm > 0,]
}

prepare_df_for_ggplot <- function(df) {
  ## Capitalize ranks
  df$rank <- capitalise_first(df$rank)
  ## Rename and reorder trait
  df$trait <- rename_traits(df$trait, as_factor=TRUE)
  df
}

extract_example_fits <- function(d) {
  clades <- c("Meliaceae", "Fagaceae")
  taxa <- sapply(d, function(x) x$BM$info$taxa)
  i <- match(clades, taxa)
  ret <- d[i]
  names(ret) <- clades
  ret
}
