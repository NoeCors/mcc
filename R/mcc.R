#' Modal Clustering for Categorical Data
#'
#' Performs modal clustering (mode-shift and level-set) on categorical data
#' as introduced by Corsini and Menardi (2025).
#'
#' @param data A data frame containing only categorical variables.
#' @param min_overlap Minimum number of shared categories required to connect two cells (default = 1).
#' @param rounding Number of decimal places used when checking for independence (default = 2).
#' @param independence_tol Tolerance level for determining independence between cells (default = 0.0001).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{original_data}{Original data frame.}
#'   \item{individual_data_with_cl}{Individual-level data with assigned cluster labels.}
#'   \item{frequency_data_with_cl}{Frequency-level data with cluster assignments.}
#'   \item{edge_list}{Data frame of edges with weights used to build the graph.}
#'   \item{graph}{igraph object representing the connectivity among categorical combinations.}
#'   \item{clustering}{Final clustering labels (if mode-shift and level-set results coincide).}
#' }
#'
#' @details
#' This function constructs a weighted directed graph whose nodes represent unique
#' categorical combinations and whose edges connect statistically dependent combinations.
#' The direction of each edge is determined by relative density. Clustering is then performed
#' using both mode-shift and level-set procedures. If both yield the same partition, that
#' is returned as the final clustering.
#'
#' @references
#' Corsini, N. and Menardi, G. (2025). *Modal Clustering for Categorical Data.*
#'
#' @examples
#' \dontrun{
#' data <- data.frame(A = sample(letters[1:3], 100, TRUE),
#'                    B = sample(LETTERS[1:2], 100, TRUE))
#' res <- mcc(data)
#' }
#'
#' @importFrom dplyr %>% select everything slice ungroup filter n mutate across group_by summarise
#' @importFrom igraph graph_from_data_frame components
#' @importFrom stats glm na.omit poisson
#' @importFrom mclust adjustedRandIndex
#'
#' @export

mcc <- function(data, min_overlap = 1, rounding = 2, independence_tol = 0.0001){

  # Ensure factors
  data <- data %>% mutate(across(everything(), as.factor))

  # Convert data to frequency format
  freq_data <- convert_to_frequency_format(data)

  n_vars <- ncol(data)
  n_cells <- nrow(freq_data)
  n_total <- sum(freq_data$y)

  # Compute observed probabilities
  freq_data$observed <- freq_data$y/n_total

  # Fit independence model
  independence_model <- glm(y ~ ., family = poisson, data = freq_data[, 1:(n_vars + 1)])
  expected <- independence_model$fitted.values

  # Check if data shows independence
  independent_cells <- which(
    round(freq_data$y/expected, rounding) >= (1 - independence_tol) &
      round(freq_data$y/expected, rounding) <= (1 + independence_tol))
  dependent_cells <- setdiff(1:n_cells, independent_cells)

  # Compute the normalized mutual information contribution for each cell
  mi_contribution <- freq_data$observed*log(freq_data$observed/(expected/n_total))
  mi_normalized <- (mi_contribution - min(mi_contribution, na.rm = TRUE)) /
    (max(mi_contribution, na.rm = TRUE) - min(mi_contribution, na.rm = TRUE))

  # Initialize edge list with self-loops
  edge_list <- data.frame(V1 = 1:n_cells, V2 = 1:n_cells, w = mi_normalized)

  # Build edges between dependent cells
  if (length(dependent_cells) > 1) {

    edge_candidates_i <- matrix(NA, nrow = 1, ncol = 3)
    edge_candidates_j <- matrix(NA, nrow = 1, ncol = 3)

    for (idx in 1:(length(dependent_cells) - 1)) {

      i <- dependent_cells[idx]

      for (j in dependent_cells[(idx + 1):length(dependent_cells)]) {

        # Count shared and different categories
        matching_categories <- which(freq_data[i, 1:n_vars] == freq_data[j, 1:n_vars])
        nomatching_categories <- setdiff(1:n_vars, matching_categories)

        if(length(matching_categories) >= min_overlap & length(nomatching_categories) >= 1) {
          weight_i <- mi_normalized[i]
          weight_j <- mi_normalized[j]
          } else {
            weight_i <- NA
            weight_j <- NA
            }

        edge_candidates_i <- rbind(edge_candidates_i, c(i, j, weight_i))
        edge_candidates_j <- rbind(edge_candidates_j, c(j, i, weight_j))
      }
    }

    # Clean up candidate edges
    rownames(edge_candidates_i) <- rownames(edge_candidates_j) <- NULL
    colnames(edge_candidates_i) <- colnames(edge_candidates_j) <- c("V1", "V2", "w")
    edge_candidates <- cbind(
      as.data.frame(edge_candidates_i[-1, ]),
      as.data.frame(edge_candidates_j[-1, ])
    )
    edge_candidates <- na.omit(edge_candidates)

    # Select edge directions based on density
    edge_directions <- assign_edge_directions(edge_candidates, mi_normalized)

    # Combine edges to keep and to discard
    combined_edges <- cbind(edge_directions$keep, edge_directions$not_keep)
    colnames(combined_edges) <- c("V1", "V2", "w1", "V1_NO", "V2_NO", "w2")

    # Select edge with minimum weight (more conservative)
    combined_edges <- cbind(combined_edges, min_weight = pmin(combined_edges$w1, combined_edges$w2))

    # Build maximal spanning tree: keep maximum weight per source node
    maximal_spanning_tree <- combined_edges %>%
      group_by(V1) %>%
      slice(which(min_weight == max(min_weight))) %>%
      select(-c(V1_NO, V2_NO)) %>%
      as.data.frame()
    final_edges <- maximal_spanning_tree %>%
      group_by(V1) %>%
      slice(which(w2 == max(w2)))  %>%
      select(-c(w2, w1)) %>%
      as.data.frame()
    colnames(final_edges) <- c("V1", "V2", "w")

    # Handle ties by preferring edges with more shared categories
    edge_ties <- final_edges %>%
      group_by(V1) %>%
      mutate(n_ties = n()) %>%
      ungroup()

    if(any(edge_ties$n_ties > 1)){
      message("Tie detected in edge selection. Resolving by choosing link among cells with maximum number of categories in common")
      final_edges <- edge_ties %>%
        mutate(n_shared = rowSums(freq_data[edge_ties$V1, 1:n_vars] == freq_data[edge_ties$V2, 1:n_vars])) %>%
        group_by(V1) %>%
        filter(n_shared == max(n_shared)) %>%
        select(-c(n_ties, n_shared)) %>%
        as.data.frame()
    } else final_edges <- edge_ties %>% select(-n_ties) %>% as.data.frame()

    # Creating the final edge list
    edge_list <- as.data.frame(rbind(edge_list[setdiff(edge_list[,1], final_edges[,1]),], final_edges))
    edge_list <- edge_list[order(edge_list[,1]),]
    rownames(edge_list) <- NULL

    # Perform level-set clustering
    level_set_clusters <- perform_level_set_clustering(edge_list, freq_data$y, n_cells, n_total)

    # Perform mean-shift clustering
    mode_shift_clusters <- perform_mode_shift(edge_list, n_cells)

    # Obtain the final graph
    final_graph <- graph_from_data_frame(edge_list, directed = TRUE)

    } else {
      message("Data shows near-perfect independence: each cell forms its own cluster")
      edge_list <- data.frame(V1 = 1:n_cells, V2 = 1:n_cells, w = mi_normalized)
      mode_shift_clusters <- 1:n_cells
      level_set_clusters <- 1:n_cells
      final_graph <- graph_from_data_frame(edge_list, directed = TRUE)
    }

  freq_data$mode_shift_cl <- mode_shift_clusters
  freq_data$level_set_cl <- level_set_clusters

  # Check that the mode-shift clustering and the level-set clustering are equal
  if(adjustedRandIndex(freq_data$mode_shift_cl, freq_data$level_set_cl) != 1){

    message("The mean-shift and level-set clustering results might not be equal")

    individual_mode_shift_cl <- get_clusters(data, freq_data[,1:n_vars], mode_shift_clusters)
    individual_level_set_cl <- get_clusters(data, freq_data[,1:n_vars], level_set_clusters)

    freq_data <- freq_data[, setdiff(names(freq_data), c("observed"))]

    individual_data <- data
    individual_data$mode_shift_cl <- individual_mode_shift_cl
    individual_data$level_set_cl <- individual_level_set_cl

    output <- list(
      original_data = data,
      individual_data_with_cl = individual_data,
      frequency_data_with_cl = freq_data,
      edge_list = edge_list,
      mode_shift_clusters = freq_data$mode_shift_cl,
      level_set_clusters = freq_data$level_set_cl,
      graph = final_graph
    )

  } else {

    # Obtain a final clustering results
    final_clustering <- rep(NA, length(mode_shift_clusters))
    unique_clusters <- sort(unique(mode_shift_clusters))
    for (i in seq_along(unique_clusters)) {
      final_clustering[mode_shift_clusters == unique_clusters[i]] <- i
    }

    individual_cl <- get_clusters(data, freq_data[,1:n_vars], final_clustering)

    freq_data <- freq_data[, setdiff(names(freq_data), c("observed", "mode_shift_cl", "level_set_cl"))]
    freq_data$final_cl <- final_clustering

    individual_data <- data
    individual_data$cl <- individual_cl

    output <- list(
      original_data = data,
      individual_data_with_cl = individual_data,
      frequency_data_with_cl = freq_data,
      edge_list = edge_list,
      clustering = final_clustering,
      graph = final_graph
    )
  }
  return(output)
}






















