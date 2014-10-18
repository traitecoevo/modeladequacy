make_data_subsets_clades <- function(rank, data) {
  min_size <- 20L
  treedata_rank(rank, data, min_size)
}

treedata_rank <- function(rank, data, min_size) {
  tr <- rank_slice_tree(rank, data$phy, min_size)
  lapply(names(tr), function(t)
         treedata_subset(tr[[t]], data, rank, t))
}

rank_slice_tree <- function(rank, tree, min_size) {
  rank <- match.arg(rank, c("family", "order"))
  pattern <- if (rank == "family") "[A-z]+ceae" else "[A-z]+ales"

  ## get all family/order nodes in tree, and extract associated subtrees:
  tax <- tree$node.label[grep(pattern, tree$node.label, perl=TRUE)]
  tree_sub <- lapply(tax, extract.clade, phy=tree)
  names(tree_sub) <- tax

  ## Only keep those that are big enough:
  ## get number of tips
  tree_sub[sapply(tree_sub, Ntip) >= min_size]
}
