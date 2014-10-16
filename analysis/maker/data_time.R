make_data_times <- function(data) {
  ## These are constant across all data sets:
  min_size <- 20L
  times <- c(0.2687, 50.2697, 100.2697, 150.2697, 200.2697)

  trait <- data$trait
  se <- data$se

  td <- lapply(times, treedata_time, data$phy, data$states, min_size)

  ## Collapse all the different time spaces together:
  td <- unlist(td, FALSE)
  for (i in seq_along(td)) {
    td[[i]] <- add_metadata(td[[i]], "timeslice", "random", trait, se)
  }

  td
}

treedata_time <- function(age, tree, data, min_size) {
  tr <- time_slice_tree(age, tree, min_size)
  lapply(tr, treedata_q, data)
}

## functions for extracting subtrees
## some of this code was written by Jon M. Eastman
time_slice_tree <- function(age, tree, min_size) {
  edge_matrix_id <- where_to_cut(age, tree)
  extract_sub_trees(tree, edge_matrix_id, min_size)
}

## Locate nonterminal edges which pass through a certain age.
where_to_cut <- function(age, tree) {
  a <- get_node_heights(tree)
  pos <- tree$edge[a[,1] < age & a[,2] > age, 2]
  pos[pos > Ntip(tree)]
}

get_node_heights <- function(t){
  edg <- as.matrix(arbutus:::edge.height(t))
  edg <- max(edg[,"start"]) - edg
  dimnames(edg) <- NULL
  # This is the height of the start and end of each node (so the first
  # row is tip 1, row Ntip(tree)+1 is the root, etc, so put it into
  # edge matrix order:
  edg[t$edge[,2],]
}

extract_sub_trees <- function(tree, edge_matrix_id, min_size, plot=FALSE) {
  node_label <- find_node_label(tree, edge_matrix_id, plot)
  extract_sub_trees_orig(tree, min_size, node_label)
}

## TODO: This seems to be used in conjunction with
## extract_sub_trees(), so I wrapped the pair up below.
find_node_label <- function(tree, edge_matrix_id, plot.tree=TRUE) {
  new_node_label <- rep("", length(tree$node.label))
  node_id <- edge_matrix_id - length(tree$tip.label)
  new_node_label[node_id] <- seq_along(node_id)
  tree$node.label <- new_node_label
  if (plot.tree) {
    plot(tree, show.tip.label=FALSE, show.node.label=TRUE)
  }
  new_node_label
}

# Extracts subclades with a species richness criterion
# Args:
#   tree: the big tree.
#   min_size: the cut off for which sub-clade to return
#   poss_nodes: the possible nodes to check
# Returns:
#   list of sub-clades that have enough species
extract_sub_trees_orig <- function(tree, min_size, poss_nodes) {
  tree_list <- list()
  tree$node.label <- poss_nodes
  pn <- poss_nodes[poss_nodes != ""]
  tree_list <- lapply(pn, extract.clade, phy=tree)
  sr <- sapply(tree_list, Ntip)
  if (!any(sr > min_size)) {
    stop("no clades within the specified age range have enough species")
  }
  tree_list[sr > min_size]
}
