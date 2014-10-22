prune_dataset_best <- function(fits) {
  ## colnames to ignore
  ign <- c("taxa", "rank", "trait", "size", "age")
  keep <- setdiff(names(fits), ign)
  tmp_df <- fits[, keep]

  ## get model for each column
  suf <- get_model_suffix(names(tmp_df))

  ## find best model
  ic <- if ("aic.bm" %in% names(fits)) "aic" else "dic"
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
