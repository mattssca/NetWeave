#' Write Network Files for Cytoscape Import
#'
#' Exports network node attributes and edges in formats ready for direct
#' import into Cytoscape. Outputs a node attribute table, an edge table,
#' and/or a GraphML file containing the full annotated network.
#'
#' @param node_metrics A data frame of network node metrics with all
#'   annotations appended by the pipeline. Must contain a \code{name} column
#'   with HGNC gene symbols.
#' @param edge_data A data frame of network edges as returned in the
#'   \code{edge_data} element of \code{create_and_expand_network()}. Must
#'   contain \code{from} and \code{to} columns.
#' @param out_dir Character string. Directory in which to write output files.
#'   Will be created if it does not exist. Defaults to \code{"."}.
#' @param file_prefix Character string. Prefix applied to all output filenames.
#'   Defaults to \code{"network"}.
#' @param format Character string controlling which files are written. One of:
#'   \describe{
#'     \item{\code{"both"}}{(Default) Writes node table, edge table, and
#'       GraphML.}
#'     \item{\code{"tables"}}{Writes \code{<prefix>_nodes.tsv} and
#'       \code{<prefix>_edges.tsv} only. Use this for Cytoscape's
#'       File → Import → Network / Table from File workflow.}
#'     \item{\code{"graphml"}}{Writes \code{<prefix>.graphml} only. All node
#'       attributes are embedded as vertex properties. Cytoscape opens this
#'       directly via File → Import → Network from File.}
#'   }
#' @param verbose Logical. If \code{TRUE} (default), progress messages and
#'   the paths of written files are printed.
#'
#' @return A named list of the file paths written:
#'   \describe{
#'     \item{node_table}{Path to the node TSV, or \code{NULL} if not written.}
#'     \item{edge_table}{Path to the edge TSV, or \code{NULL} if not written.}
#'     \item{graphml}{Path to the GraphML file, or \code{NULL} if not written.}
#'   }
#'
#' @details
#' \strong{Node table} (\code{<prefix>_nodes.tsv}): All columns from
#' \code{node_metrics}, with \code{name} as the key column. Import into
#' Cytoscape via File → Import → Table from File, setting the key column
#' to \code{name}.
#'
#' \strong{Edge table} (\code{<prefix>_edges.tsv}): The \code{from} and
#' \code{to} columns from \code{edge_data}, filtered to edges where both
#' endpoints are present in \code{node_metrics$name}. Import as a network
#' via File → Import → Network from File.
#'
#' \strong{GraphML} (\code{<prefix>.graphml}): A single file containing
#' both the network topology and all node attributes. Logical columns are
#' coerced to integer (\code{0}/\code{1}) for broad GraphML parser
#' compatibility. Open directly in Cytoscape via File → Import → Network
#' from File.
#'
#' If \code{node_metrics} contains duplicate \code{name} values (e.g. from
#' multiple HPA cell type rows), only the first occurrence of each gene is
#' retained and a warning is issued.
#'
#' @import igraph
#' @importFrom utils write.table
#'
#' @export
create_cytoscape_objects <- function(node_metrics = NULL,
                                     edge_data = NULL,
                                     out_dir = ".",
                                     file_prefix = "network",
                                     format = c("both", "tables", "graphml"),
                                     verbose = TRUE){
  
  format <- match.arg(format)
  
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
  
  # handle duplicate gene rows (e.g. multiple HPA cell types)
  if(anyDuplicated(node_metrics$name)){
    n_dupes <- sum(duplicated(node_metrics$name))
    warning(sprintf(
      "%d duplicate name(s) found in node_metrics — keeping first occurrence per gene. Consider summarising HPA cell type rows before export.",
      n_dupes
    ))
    node_metrics <- node_metrics[!duplicated(node_metrics$name), ]
  }
  
  # create output directory if needed
  if(!dir.exists(out_dir)){
    dir.create(out_dir, recursive = TRUE)
    if(verbose) message(sprintf("  -> Created output directory: %s", out_dir))
  }
  
  retained_genes  <- node_metrics$name
  filtered_edges  <- edge_data[
    edge_data$from %in% retained_genes & edge_data$to %in% retained_genes, 
  ]
  
  if(verbose){
    message("Writing Cytoscape files...")
    message(sprintf("  -> Nodes: %d", length(retained_genes)))
    message(sprintf("  -> Edges: %d", nrow(filtered_edges)))
    message(sprintf("  -> Format: %s", format))
  }
  
  paths <- list(node_table = NULL, edge_table = NULL, graphml = NULL)
  
  write_tsv <- function(df, path){
    write.table(df, file = path, sep = "\t", quote = FALSE,
                row.names = FALSE, na = "")
  }
  
  # --- node table ---
  if(format %in% c("both", "tables")){
    node_path <- file.path(out_dir, paste0(file_prefix, "_nodes.tsv"))
    write_tsv(node_metrics, node_path)
    paths$node_table <- node_path
    if(verbose) message(sprintf("  -> Node table written: %s", node_path))
  }
  
  # --- edge table ---
  if(format %in% c("both", "tables")){
    edge_path <- file.path(out_dir, paste0(file_prefix, "_edges.tsv"))
    write_tsv(filtered_edges, edge_path)
    paths$edge_table <- edge_path
    if(verbose) message(sprintf("  -> Edge table written: %s", edge_path))
  }
  
  # --- graphml ---
  if(format %in% c("both", "graphml")){
    
    graphml_path <- file.path(out_dir, paste0(file_prefix, ".graphml"))
    
    # rebuild graph from filtered edges; add isolated nodes
    g            <- graph_from_data_frame(filtered_edges, directed = FALSE)
    missing_nodes <- setdiff(retained_genes, V(g)$name)
    if(length(missing_nodes) > 0){
      g <- g + igraph::vertices(missing_nodes)
    }
    
    # attach node attributes from node_metrics
    attr_cols <- setdiff(colnames(node_metrics), "name")
    
    for(col in attr_cols){
      vals <- node_metrics[[col]]
      
      # coerce logical to integer for GraphML compatibility
      if(is.logical(vals)){
        vals <- as.integer(vals)
      }
      
      # match to vertex order in graph
      vals_ordered <- vals[match(V(g)$name, node_metrics$name)]
      g <- set_vertex_attr(g, name = col, value = vals_ordered)
    }
    
    write_graph(g, file = graphml_path, format = "graphml")
    paths$graphml <- graphml_path
    if(verbose) message(sprintf("  -> GraphML written:   %s", graphml_path))
  }
  
  if(verbose){
    message("  -> Done.")
  }
  
  return(invisible(paths))
  
}
