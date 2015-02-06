## Creates the split data sets
make_data_subsets <- function(data) {
  ret <- c(make_data_subsets_times(data),
           make_data_subsets_clades("order", data),
           make_data_subsets_clades("family", data))
  names(ret) <- NULL
  ret
}

treedata_subset <- function(tr, data, rank, taxon) {
  x_log <- treedata_q(tr, data$states_log)
  states_raw <- x_log$data
  states_raw[] <- data$states_raw[rownames(x_log$data)]
  list(phy=x_log$phy,
       states_log=x_log$data, states_raw=states_raw,
       se_log=data$se_log, se_raw=data$se_raw,
       trait=data$trait, rank=rank, taxa=taxon)
}
