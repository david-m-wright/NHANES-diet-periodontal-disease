# Helper functions for treelet transform analysis

require(haven)
require(tidyverse)
#require(sca)
require(treelet)
require(ggdendro)

## Functions ##

# Create a dendrogram from a TT tree
# Args: object output from TT tree build with Build_JTree or Run_JTree
# Value: object of class ("dendro")
# Note that only a full height tree can be plotted
ConvertTTDendro <- function(tt_tree){
  
  # Maximum height of the tree
  maxlev <- length(tt_tree$TreeCovs)
  
  # Matrix of node merges
  zpd <- tt_tree$Zpos
  # Unique node values in the matrix
  node_values <- unique(as.numeric(zpd))
  
  # Convert matrix of node merges to format used for 'hclust' trees
  # Mark original variables as negative 
  # Set repeated node values to equal the clustering stage at which they occured (tree height)
  node_addr <- vector(mode = "list", length = length(node_values))
  
  for(i in 1:length(node_addr)){
    node_addr[[i]] <- which(zpd == node_values[i], arr.ind = T)  
  }
  
  for(j in 1:length(node_addr)){
    if(nrow(node_addr[[j]]) > 1){
      zpd[node_addr[[j]][order(node_addr[[j]][, "row"]),]] <- c(node_values[j] * -1, sort(node_addr[[j]][, "row"])[1:nrow(node_addr[[j]])-1])
      
    } else {
      zpd[node_addr[[j]]] <- node_values[j] * -1
    }
  }
  
  # Create hclust object
  ttd <- list()
  ttd$merge <- zpd
  ttd$height <- 1:maxlev
  ttd$order <- 1:ncol(tt_tree$TreeCovs[[maxlev]])
  ttd$labels <- rownames(tt_tree$TreeCovs[[maxlev]])
  class(ttd) <- "hclust"
  
  #plot(ttd)
  # Convert to dendrogram
  as.dendrogram(ttd)
  
}


# Function to create summary table of variances for treelet components
# Args: tt_tree = treelet transform tree produced by Build_JTree() or Run_JTree()
#       c_matrix = correlation or covariance matrix on which the tree was computed
#       cut_level = height at which to cut the tree and extract the basis functions, defaults to max
#       components = number of components to display, defaults to the maximum
# Value: a tibble with the components and associated variances
# Note that the ordering of the TCs with identical variance may differ between the standard
# and adjusted variance calculations (Corrected Sum of Variance from 'sca' package)
TTVariances <- function(tt_tree, c_matrix, cut_level = NULL, components = NULL){
  
  # Maximum tree height
  if(is.null(cut_level)){
    cut_level <- length(tt_tree$TreeCovs)
  } 
  
  # Variance associated with each Treelet Component at the selected cut level
  tc_v <- diag(tt_tree$TreeCovs[[cut_level]]) 
  # Decreasing order of variances
  tc_order <- order(tc_v, decreasing = T)
  
  # Summary table of variances
  tc_var_summary <- enframe(tc_v[tc_order]) %>% 
    transmute(Component = 1:nrow(.),
              Variance = value,
              Proportion = Variance/sum(Variance),
              Cumulative = cumsum(Proportion))
              # Variances adjusted for correlations between TCs
              #scc = sccrit(c_matrix, tt_tree$basis[[cut_level]], criterion = "csv", sortP = TRUE),
              #adj_proportion = c(scc[1], diff(scc))) %>% 
    #select(-scc)
  
  
  # Extract the treelet components (TCs)
  # Each column in the basis matrix is a TC.
  tcomps <- tt_tree$basis[[cut_level]][, tc_order]
  rownames(tcomps) <- rownames(tt_tree$TreeCovs[[cut_level]])
  colnames(tcomps) <- paste0("TC", 1:ncol(tcomps))
  
  # Components to report
  if(!is.null(components)){
    tcomps <- tcomps[,1:components]
  }
  
  list(tc = tcomps,
       tcv = tc_var_summary)
}


# Function to cross validate a treelet transform
# Args: dataset = dataset on which to perform the cross validation
#       m = number of components to extract
#       nfolds = number of folds to use
# Value: vector of average cross validation scores across cut levels
# Note that this constructs the tree using the correlation matrix and 
# so only scaled variables should be used
CrossValidateTT <- function(dataset, m, nfolds = 10){
  dataset <- as.matrix(dataset)
  # Define the folds
  cv_folds <- sample(ntile(1:nrow(dataset), nfolds))
  # Maximum tree height
  maxlev <- ncol(dataset) -1
  # Initialise matrix for cross validation scores at each level
  cv_scores <- matrix(NA, nrow = maxlev, ncol = nfolds, dimnames = list(paste0("level", 1:maxlev), NULL))
  
  # Build tree for all folds except the focal one
  for(i in 1:nfolds){
    tt_data <- dataset[cv_folds != i, ]
    tt_tree <- Run_JTree(cor(tt_data), maxlev = maxlev, whichsave = 1:maxlev)
    
    # Extract the m highest variance components from each cut level
    comp_var <- vector("list", length = maxlev)
    for(j in 1:maxlev){
      m_comps <- tt_tree$basis[[j]][, order(diag(tt_tree$TreeCovs[[j]]), decreasing = T)[1:m]]
      # Calculate sum of variances based on these components in focal fold
      comp_var[[j]] <- sum(diag(cov(dataset[cv_folds == i, ] %*% m_comps)))
    }
    
    cv_scores[,i] <- as.numeric(comp_var)
  }
  apply(cv_scores, 1, mean)
}
