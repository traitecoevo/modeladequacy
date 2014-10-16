make_data_families <- function(data) {
  make_data_clades("family", data)
}

make_data_orders <- function(data) {
  make_data_clades("order", data)
}

make_data_clades <- function(rank, data) {
  min_size <- 20L

  tree <- data$phy
  states <- data$states
  
  td <- treedata_taxon(rank, tree, states, min_size)
  trait <- data$trait
  se <- data$se

  for (taxa in names(td)) {
    td[[taxa]] <- add_metadata(td[[taxa]], rank, taxa, trait, se)
  }

  td
}

treedata_taxon <- function(rank, tree, data, min_size) {
  rank <- match.arg(rank, c("family", "order"))
  pattern <- if (rank == "family") "[A-z]+ceae" else "[A-z]+ales"

  ## get all family/order nodes in tree
  tax <- tree$node.label[grep(pattern, tree$node.label, perl=TRUE)]
  ## extract subtrees
  trees <- lapply(tax, function(x) extract.clade(tree, node=x))
  ## use treedata to match to family level
  td <- lapply(trees, treedata_q, data)
  names(td) <- tax
  ## get number of tips
  tips <- lapply(td, function(x) Ntip(x$phy))
  ## extract those which meet threshold
  td[tips >= min_size]
}
