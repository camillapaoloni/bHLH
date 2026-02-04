library(ape)
library(ggtree)
project_root <- Sys.getenv("BHLH_PROJECT_ROOT", unset = ".")
p <- function(...) file.path(project_root, ...)


# 1. Read the Newick tree
tree <- read.tree(p("data", "raw", "Tree23sp.newick"))

# 2. Order species according to species_order_rev
# rotateConstr orders the tree to satisfy the given order
tree_ordered <- rotateConstr(tree, species_order_rev)

# 3. Function to find the MRCA of a clade and assign the clade name to that node
add_clade_labels <- function(tree, clades) {
  # Create a vector of NA for all internal nodes
  node_labels <- rep(NA, tree$Nnode)
  
  # Internal nodes in ape are indexed from (Ntip + 1) to (Ntip + Nnode)
  internal_nodes <- (Ntip(tree) + 1):(Ntip(tree) + tree$Nnode)
  
  # For each clade
  for (clade_name in names(clades)) {
    species <- clades[[clade_name]]
    
    # MRCA (most recent common ancestor) of this clade
    mrca_node <- getMRCA(tree, species)
    
    if (!is.null(mrca_node)) {
      # Index in the node_labels vector
      idx <- mrca_node - Ntip(tree)
      node_labels[idx] <- clade_name
    }
  }
  
  # Assign clade names to internal nodes (NA if no name)
  tree$node.label <- node_labels
  return(tree)
}

# 4. Apply the function
tree_labeled <- add_clade_labels(tree_ordered, clades)

# 5. Plot with ggtree (with labeled nodes)
p <- ggtree(tree_labeled) + 
  geom_tiplab() + 
  geom_text2(aes(subset = !isTip, label = node.label), hjust = -0.2, vjust = -0.5, fontface = "bold")

print(p)
