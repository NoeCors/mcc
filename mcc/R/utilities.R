# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Select Edge Direction Based on Node Density
#'
#' For each edge connecting two nodes with different densities, selects the
#' direction from lower density to higher density node.
#'
#' @param edge_list Data frame with columns V1, V2, w (first node, second node, weight)
#' @param node_density Numeric vector of density values for each node
#'
#' @return A list with:
#' \describe{
#'   \item{keep}{Data frame of edges retained with assigned direction.}
#'   \item{not_keep}{Data frame of edges discarded (opposite directions).}
#' }
#' @keywords internal
#'
assign_edge_directions <- function(edge_list, node_density){

  v1 <- edge_list[,1]
  v2 <- edge_list[,2]

  # Classify edge relationships
  equal <- node_density[v1] == node_density[v2]
  smaller <- node_density[v1] < node_density[v2]
  greater <- !equal & !smaller

  edges_to_keep <- rbind(
    edge_list[equal, c(1,2,3)],
    edge_list[equal,c(4,5,6)],
    edge_list[smaller, c(1,2,3)],
    edge_list[greater, c(4,5,6)]
  )

  edges_not_keep <- rbind(
    edge_list[equal,c(4,5,6)],
    edge_list[equal,c(1,2,3)],
    edge_list[smaller, c(4,5,6)],
    edge_list[greater, c(1,2,3)]
  )

  return(list(
    keep = as.data.frame(edges_to_keep),
    not_keep = as.data.frame(edges_not_keep)
  ))
}

#' Convert Individual-Level Data to Frequency Format
#'
#' Aggregates individual observations into frequency counts for each unique
#' combination of variable levels
#'
#' @param data Data frame with categorical variables
#' @return A data frame containing all unique combinations of variable values,
#' along with a frequency column \code{y} giving the count for each combination.
#'
#' @examples
#' df <- data.frame(A = c(1,1,2,2), B = c("x","x","y","y"))
#' convert_to_frequency_format(df)
#'
#' @importFrom dplyr %>% n select slice everything ungroup mutate filter across group_by summarise
#' @importFrom stats glm na.omit poisson
#' @importFrom igraph graph_from_data_frame components
#'
#' @export
#'
convert_to_frequency_format <- function(data) {
  data %>%
    mutate(y = 1) %>%
    group_by(across(everything())) %>%
    summarise(y = n(), .groups = 'drop') %>%
    as.data.frame()
}


#' Perform Mode-Shift Algorithm
#'
#' Iteratively reassigns each node to its most connected higher-density neighbor
#'  until convergence. The final assignments identify modal regions in the density.
#'
#' @param edge_list A data frame with V1, V2 columns defining node connections
#' @param n_nodes Integer specifying the total number of nodes in the graph
#'
#' @return Vector of cluster assignments.
#' @keywords internal
perform_mode_shift <- function(edge_list, n_nodes) {

  path_history <- list()
  converged <- FALSE
  iteration <- 1

  # Initialize: each node points to itself
  path_history[[iteration]] <- 1:n_nodes

  # Iteratively shift to higher-density neighbors
  while (!converged) {
    iteration <- iteration + 1

    # Shift each node to its target (V2) in edge list
    path_history[[iteration]] <- edge_list$V2[path_history[[iteration - 1]]]

    # Check convergence: no change from last iteration
    converged <- all(path_history[[iteration]] == path_history[[iteration - 1]])

    # Also check for cycles (return to earlier state)
    if (iteration > 2 && !converged) {
      converged <- any(unlist(lapply(path_history[-iteration], function(x) all(path_history[[iteration]] == x))))
    }
  }

  return(path_history[[iteration]])
}


#' Perform Level-Set Clustering
#'
#' Identifies connected components at different density threshold levels,
#' creating a hierarchical clustering structure.
#'
#' @param edge_list A data frame with V1, V2, w columns representing the graph
#' edges and their associated weights (densities)
#' @param freq_counts Numeric vector of frequency counts corresponding to each node.
#' @param n_cells Total number of nodes (categorical combinations).
#' @param n_total Total number of observations in the dataset.
#'
#' @return Vector of cluster assignments.
#' @keywords internal
perform_level_set_clustering <- function(edge_list, freq_counts, n_cells, n_total) {

  thresholds <- sort(unique(as.numeric(edge_list$w)), decreasing = TRUE)
  n_thresholds <- length(thresholds)

  # Initialize output structure
  numcon <- list()
  numcon$id <- matrix(-1, n_cells, n_thresholds + 2)
  numcon$p <- numeric(n_thresholds + 2)
  numcon$nc <- numeric(n_thresholds + 2)

  # Level-set algorithm: threshold at each unique weight value
  for (i in 1:n_thresholds) {

    # Extract edges above current threshold
    subgraph_edges <- edge_list[which(edge_list$w >= thresholds[i]),]
    subgraph <- graph_from_data_frame(d = subgraph_edges, directed = TRUE)

    # Identify connected components
    temp <- components(subgraph)

    numcon$id[as.numeric(names(temp$membership)), i + 1] <- temp$membership
    numcon$p[i + 1] <- sum(freq_counts[as.matrix(subgraph_edges[,1:2])])/n_total
    numcon$nc[i + 1] <- temp$no
  }

  # Add boundary conditions
  numcon$nc[1] <- 0
  numcon$id[,ncol(numcon$id)] <- 1
  numcon$nc[length(numcon$nc)] <- 1

  # Building the tree
  tree <- con2tree(numcon)

  return(tree$g)
}

#' Convert Level-Set Results to Level-Set Tree
#'
#' Constructs a level-set tree structure from level-set clustering results.
#' This is a complex recursive algorithm that builds the clustering hierarchy.
#'
#' @param object Output from perform_level_set_clustering()
#'
#' @return A list with:
#' \describe{
#'   \item{g}{Vector of final group memberships for each node.}
#'   \item{tree}{A \code{dendrogram} object describing the hierarchy.}
#'   \item{bad}{Logical indicating numerical instability in height estimation.}
#'   \item{noc}{Number of clusters detected.}
#' }
#' @keywords internal

con2tree <- function(object, f = NULL){
  # object is the output of num.con()
  # f density estimate
  ow <- options("warn")
  nc <- object$nc
  p <- object$p
  index <- which(diff(nc) != 0) # posizione punti di salto
  K <- length(index)  					# numero di punti di salto
  ps <- p[index]      					# frazione di dati inclusi ai vari punti di salto
  lista <- list()
  if (K == 1){
    gruppi <- as.vector(object$id[,ncol(object$id)])
    M <- 1 } else {
      for(j in 1:K) lista[[j]] <- as.vector(object$id[,index[j]])  # list of connected sets at the jumps
      gruppi <- rep(0, length(lista[[1]]))
      M <- 0
      insiemi <- list()
      # step
      k <- 1
      allocated <- (gruppi > 0)
      insiemi[[k]] <- setdiff(unique(gruppi), 0)
      # loop k = 2,...,K
      while(k < K) {
        k <- k + 1
        sets  <- lista[[k]]                #elenco connessi al salto k
        insieme <- list()
        for(m in 1:max(sets)) {
          set <- which(sets == m)
          # which objects are in connected component m?
          new <- setdiff(set, which(allocated))
          # which objects in connected component m are not allocated yet?
          if(length(new) > 0){
            # are there new objects in the connected component?
            g.new <- unique(gruppi[intersect(set, which(allocated))])
            # which is the label group of objects in connected component m and has been already allocated?
            if(length(g.new) == 0)  gruppi[set] <- M <- M + 1
            if(length(g.new) == 1)  gruppi[set] <- g.new
            allocated <- (gruppi > 0)
          }
          gg <- sort(setdiff(unique(gruppi[set]), 0))
          if(length(gg) > 0) insieme[[length(insieme) + 1]] <- gg
        }
        insiemi[[k]] <- insieme
      }
      g <- gruppi
      # if(!missing(f)){
      #   u <- unique(gruppi[rev(order(f))])
      #   g <- rep(0, length(f))
      #   u0 <- u[u > 0]
      #   for(i in 1:max(gruppi)) g[gruppi == u0[i]] <- i
      #   gruppi <- g
      # }
    }
  salti <- diff(nc)
  salta.giu <- rev(which(diff(nc)[index] < 0))
  altezza <- numeric(M)
  m <- 0
  salti.su <- salti[salti > 0]
  options(warn = -1)
  while(m < M){
    m <- m + 1
    r <- min(which(cumsum(salti.su) >= m))
    altezza[m] <- p[salti > 0][r]
  }
  bad.grid <- any(is.na(altezza))
  sotto.albero <- function(tree, k, set){
    insieme <- insiemi[[salta.giu[k]]]
    r <- 0
    branch <- list()
    for(item0 in insieme){
      item <- intersect(set, unlist(item0))
      if(length(item) == 1){
        r <- r + 1
        u <- item
        attr(u, "members") <- 1
        attr(u, "height") <- altezza[item]
        attr(u, "label") <- paste(as.character(item), " ", sep = "")
        attr(u, "leaf") <- TRUE
        branch[[r]] <- u
      }
      if(length(item) > 1) {
        r <- r + 1
        u <- sotto.albero(list(), k + 1, item)
        attr(u, "members") <- length(unlist(item))
        attr(u, "height") <- max(ps[salta.giu[k + 1]])
        attr(u, "label") <- paste("{", paste(item, collapse = ","), "}", sep = "")
        attr(u, "leaf") <- FALSE
        branch[[r]] <- u
        # browser()
      }
    }
    branch
  }
  if(M > 1) {
    tree <- sotto.albero(list(), 1, 1:M)
    attr(tree, "members") <- M
    attr(tree, "height") <-  max(ps[salta.giu[1]])
    attr(tree, "label") <- paste("{", paste(1:M, collapse = ","),"}", sep = "")
    noc <- M } else {
      tree <- list()
      tree[[1]]<-1
      attr(tree[[1]], "members") <- 1
      attr(tree[[1]], "height") <- 0
      attr(tree[[1]], "label") <- "1"
      attr(tree[[1]], "leaf") <- TRUE
      attr(tree, "members") <- 1
      attr(tree, "height") <- 1
      attr(tree, "label") <- paste("{", paste(1, collapse = ","),"}", sep = "")
      noc <- 1
    }
  tree <- list(tree)
  attr(tree, "members") <- M
  attr(tree, "height") <- 1
  attr(tree, "class") <- "dendrogram"
  options(warn = ow$warn)
  invisible(list(g = gruppi, tree = tree, bad = bad.grid, noc = noc))
}


#' Get Individual Cluster Labels from Frequency Data
#'
#' Helper function that maps cluster assignments obtained at the frequency level
#' back to individual observations by matching their categorical combinations.
#'
#' @param individual_data Original data frame with individual observations.
#' @param frequency_data Data frame containing unique categorical combinations.
#' @param cl_freq Vector of cluster labels corresponding to the rows in \code{frequency_data}.
#'
#' @return A vector of cluster labels aligned with the rows of \code{individual_data}.
#'
#' @export
get_clusters <- function(individual_data, frequency_data, cl_freq) {

  # Preliminary checks
  if (ncol(individual_data) != ncol(frequency_data)) {
    stop("Individual data and Frequency data must have the same number of columns")
  }

  if (length(cl_freq) != nrow(frequency_data)) {
    stop("length(Cluster column in Frequency data) must be equal to number of rows of Frequency data")
  }

  # Convert rows to strings for matching
  ind_strings <- apply(individual_data, 1, paste, collapse = "_")
  freq_strings <- apply(frequency_data, 1, paste, collapse = "_")

  # Find matching positions
  positions <- match(ind_strings, freq_strings)

  # Return corresponding cluster labels
  return(cl_freq[positions])
}


