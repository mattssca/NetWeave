#' Recalculate Network Metrics on a Filtered Node Set
#'
#' After filtering \code{node_metrics} (e.g., with
#' \code{filter_network_nodes()}), the original centrality values reflect the
#' full pre-filter network. This function rebuilds the \code{igraph} object
#' from the surviving nodes and recomputes all structural metrics.
#'
#' @param node_metrics A data frame of network node metrics, already filtered
#'   to the genes of interest. Must contain a \code{name} column with HGNC
#'   gene symbols.
#' @param edge_data A data frame of network edges as returned in the
#'   \code{edge_data} element of \code{create_and_expand_network()}. Must
#'   contain \code{from} and \code{to} columns.
#' @param verbose Logical. If \code{TRUE} (default), progress messages
#'   including node and edge counts before and after filtering are printed.
#'
#' @return The \code{node_metrics} data frame with the following columns
#'   replaced with values recalculated on the filtered subnetwork:
#'   \code{degree}, \code{betweenness}, \code{closeness}, \code{eigenvector},
#'   \code{hub_score}, and \code{community}.
#'
#' @details
#' Edges are filtered to retain only those where both endpoints exist in
#' \code{node_metrics$name}. The \code{igraph} graph is then rebuilt from
#' this reduced edge list and all metrics recomputed. Genes not connected
#' to any other node in the filtered set (isolates) will have zero degree,
#' betweenness, and closeness. The \code{origin} column and all annotation
#' columns are preserved unchanged.
#'
#' @seealso \code{\link{create_and_expand_network}}, \code{\link{filter_network_nodes}}
#'
#' @import igraph
#'
#' @export
recalculate_network_metrics <- function(node_metrics = NULL,
                                        edge_data = NULL,
                                        verbose = TRUE){
  
  # checks
  if(is.null(node_metrics)){
    stop("User must provide an incoming object with network node metrics...")
  }
  
  if(!"name" %in% colnames(node_metrics)){
    stop("node_metrics must contain a 'name' column...")
  }
  
  if(is.null(edge_data)){
    stop("User must provide edge_data from create_and_expand_network()...")
  }
  
  if(!all(c("from", "to") %in% colnames(edge_data))){
    stop("edge_data must contain 'from' and 'to' columns...")
  }
  
  required_metric_cols <- c("degree", "betweenness", "closeness",
                            "eigenvector", "hub_score", "community")
  missing_cols <- required_metric_cols[!required_metric_cols %in% colnames(node_metrics)]
  if(length(missing_cols) > 0){
    stop(sprintf(
      "node_metrics is missing expected metric column(s): %s",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  retained_genes <- node_metrics$name
  
  if(verbose){
    message("Recalculating network metrics on filtered subnetwork...")
    message(sprintf("  -> Retained nodes: %d", length(retained_genes)))
    message(sprintf("  -> Original edges: %d", nrow(edge_data)))
  }
  
  # filter edges to only those where both endpoints are retained
  filtered_edges <- edge_data[
    edge_data$from %in% retained_genes & edge_data$to %in% retained_genes, 
  ]
  
  if(nrow(filtered_edges) == 0){
    stop("No edges remain after filtering to retained nodes. Cannot recalculate metrics...")
  }
  
  if(verbose){
    message(sprintf("  -> Filtered edges: %d", nrow(filtered_edges)))
  }
  
  # rebuild graph from filtered edges
  g <- graph_from_data_frame(filtered_edges, directed = FALSE)
  
  # ensure all retained genes are present as vertices, even if isolated
  missing_nodes <- setdiff(retained_genes, V(g)$name)
  if(length(missing_nodes) > 0){
    g <- g + igraph::vertices(missing_nodes)
    if(verbose){
      message(sprintf("  -> Added %d isolated nodes (no edges in filtered network)",
                      length(missing_nodes)))
    }
  }
  
  if(verbose){
    message(sprintf("  -> Recalculating metrics for %d nodes, %d edges",
                    length(V(g)), length(E(g))))
  }
  
  # recompute all metrics
  hits_result  <- hits_scores(g, scale = TRUE)
  eigen_scores <- eigen_centrality(g)$vector
  
  new_metrics <- data.frame(
    name        = V(g)$name,
    degree      = igraph::degree(g),
    betweenness = betweenness(g),
    closeness   = closeness(g),
    eigenvector = eigen_scores,
    hub_score   = hits_result$hub,
    community   = as.numeric(membership(cluster_louvain(g))),
    stringsAsFactors = FALSE
  )
  
  # drop old metric columns from node_metrics and join fresh values
  node_metrics <- node_metrics[, !colnames(node_metrics) %in% required_metric_cols]
  node_metrics <- merge(node_metrics, new_metrics, by = "name", all.x = TRUE)
  
  # restore original column order: name, origin, metrics, then everything else
  core_cols    <- c("name", "origin", required_metric_cols)
  core_cols    <- core_cols[core_cols %in% colnames(node_metrics)]
  other_cols   <- setdiff(colnames(node_metrics), core_cols)
  node_metrics <- node_metrics[, c(core_cols, other_cols)]
  
  # restore original row order
  node_metrics <- node_metrics[match(retained_genes, node_metrics$name), ]
  rownames(node_metrics) <- NULL
  
  if(verbose){
    message(sprintf("  -> Done. %d nodes with recalculated metrics.", nrow(node_metrics)))
  }
  
  return(node_metrics)
  
}
