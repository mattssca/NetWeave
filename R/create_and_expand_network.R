#' Create and Expand a Protein Interaction Network from Pathway Genes
#'
#' Builds a gene interaction network using STRING database protein-protein
#' interaction data. Starting from a set of pathway genes, the network is
#' expanded by adding the highest-degree STRING neighbours that are also
#' present in the supplied expression data. Node-level network metrics are
#' computed and returned alongside the edge list.
#'
#' @param expr_data A numeric matrix or data frame of expression values with
#'   genes as rows (rownames = HGNC gene symbols) and samples as columns.
#' @param pathway_genes A character vector of HGNC gene symbols defining the
#'   seed pathway. Must not be \code{NULL}.
#' @param seed_gene An optional single character string naming a specific gene
#'   of interest within the pathway. Used to label the \code{origin} field in
#'   node metrics as \code{"seed"}. If \code{NULL}, no gene is labelled as
#'   seed.
#' @param max_added_genes Integer. Maximum number of high-degree neighbour
#'   genes to add when expanding the network beyond the pathway genes.
#'   Defaults to \code{20}.
#' @param string_score_threshold Integer. Minimum combined STRING interaction
#'   score (0–1000) to retain an edge. Defaults to \code{400}.
#' @param verbose Logical. If \code{TRUE} (default), progress messages are
#'   printed at each step.
#' @param genes_blacklist An optional character vector of gene symbols to
#'   exclude from the neighbour expansion step. Pathway seed genes are not
#'   affected.
#' @param string_data_dir Character string path to a directory containing
#'   pre-downloaded STRING flat files. If the required files are present,
#'   they are used instead of downloading. Defaults to
#'   \code{"data/string_db"}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{edge_data}{A data frame of network edges (from/to STRING
#'       interaction pairs) as returned by \code{igraph::as_data_frame}.}
#'     \item{node_metrics}{A data frame with one row per node containing:
#'       \code{name} (HGNC symbol), \code{origin} (\code{"seed"},
#'       \code{"pathway"}, \code{"extended set"}, or \code{"unknown"}),
#'       \code{degree}, \code{betweenness}, \code{closeness},
#'       \code{eigenvector}, \code{hub_score}, and \code{community}
#'       (Louvain cluster membership).}
#'   }
#'
#' @details
#' The function proceeds in five steps:
#' \enumerate{
#'   \item Initialise the STRING database via \code{STRINGdb}, using cached
#'     local files if all required \code{.txt.gz} files are present in
#'     \code{string_data_dir}.
#'   \item Map pathway genes to STRING protein IDs.
#'   \item Retrieve STRING first-degree neighbours, convert them to HGNC
#'     symbols using the STRING protein info table (preferred_name field),
#'     filter to genes present in \code{expr_data}, optionally remove
#'     blacklisted genes, and rank by STRING interaction degree to select
#'     the top \code{max_added_genes}.
#'   \item Build an \code{igraph} undirected graph from STRING interactions
#'     among the expanded gene set, retaining only nodes present in
#'     \code{expr_data}.
#'   \item Compute per-node metrics: degree, betweenness, closeness,
#'     eigenvector centrality, HITS hub score, and Louvain community.
#' }
#'
#' Requires an internet connection only if STRING files are not cached locally.
#' Gene symbol mapping uses \code{string_db$get_proteins()} which reads the
#' cached STRING protein info file, guaranteeing full coverage and consistency
#' with the STRING IDs used throughout.
#'
#' @import STRINGdb igraph
#' @importFrom utils capture.output
#'
#' @export
create_and_expand_network <- function(expr_data,
                                      pathway_genes = NULL,
                                      seed_gene = NULL,
                                      max_added_genes = 20,
                                      string_score_threshold = 400,
                                      verbose = TRUE,
                                      genes_blacklist = NULL,
                                      string_data_dir = "data/string_db"){
  
  # checks
  if(is.null(pathway_genes) || length(pathway_genes) == 0){
    stop("User must provide a non-empty character vector of pathway genes...")
  }
  
  if(!is.null(seed_gene) && (!is.character(seed_gene) || length(seed_gene) != 1)){
    stop("seed_gene must be a single character string or NULL...")
  }
  
  if(!is.null(genes_blacklist) && !is.character(genes_blacklist)){
    stop("genes_blacklist must be a character vector or NULL...")
  }
  
  # get available genes in expression data
  available_genes <- rownames(expr_data)
  
  # filter pathway genes to those present in expression data
  available_pathway_genes <- intersect(pathway_genes, available_genes)
  
  if(length(available_pathway_genes) == 0){
    stop("None of the pathway genes are present in the expression data...")
  }
  
  if(verbose && length(available_pathway_genes) < length(pathway_genes)){
    message(sprintf("  -> Warning: %d pathway genes not found in expression data",
                    length(pathway_genes) - length(available_pathway_genes)))
  }
  
  # step 1: initialise STRING database
  if(verbose) message("Initializing STRING database...")
  
  required_files <- c(
    "9606.protein.aliases.v11.5.txt.gz",
    "9606.protein.info.v11.5.txt.gz",
    "9606.protein.links.v11.5.txt.gz"
  )
  
  use_cache <- dir.exists(string_data_dir) &&
    all(required_files %in% list.files(string_data_dir, pattern = "\\.txt\\.gz$"))
  
  if(use_cache){
    if(verbose) message("  -> Using cached STRING files")
    string_db <- STRINGdb$new(version = "11.5",
                              species = 9606,
                              score_threshold = string_score_threshold,
                              input_directory = string_data_dir)
  } else {
    if(verbose) message(ifelse(dir.exists(string_data_dir),
                               "  -> Cached files incomplete, downloading...",
                               "  -> No cached files found, downloading..."))
    string_db <- STRINGdb$new(version = "11.5",
                              species = 9606,
                              score_threshold = string_score_threshold,
                              input_directory = "")
  }
  
  # build a complete STRING protein_id -> HGNC symbol lookup from the
  # STRING protein info table — guaranteed full coverage, no external db needed
  if(verbose) message("  -> Building STRING ID -> symbol lookup...")
  all_proteins <- string_db$get_proteins()
  id_to_symbol <- setNames(all_proteins$preferred_name,
                           all_proteins$protein_external_id)
  
  # step 2: map pathway genes to STRING IDs
  if(verbose) message("Mapping pathway genes to STRING IDs...")
  
  mapped <- suppressWarnings(suppressMessages(
    string_db$map(data.frame(gene = available_pathway_genes), "gene", removeUnmappedRows = TRUE)
  ))
  signature_ids <- mapped$STRING_id
  
  if(verbose) message(sprintf("  -> Mapped %d genes", length(signature_ids)))
  
  # step 3: find and filter neighbours
  if(verbose) message("Finding neighbor genes...")
  
  neighbors <- tryCatch(string_db$get_neighbors(signature_ids),
                        error = function(e) character(0))
  
  if(length(neighbors) == 0 || all(is.na(neighbors))){
    stop("No neighbors found for pathway genes in STRING...")
  }
  
  neighbors <- as.character(neighbors[!is.na(neighbors) & neighbors != ""])
  if(length(neighbors) == 0) stop("No valid STRING neighbor IDs...")
  
  # convert neighbour STRING IDs to gene symbols using the pre-built lookup
  if(verbose) message("  -> Mapping neighbor IDs to gene symbols...")
  neighbor_genes <- unique(id_to_symbol[neighbors])
  neighbor_genes <- neighbor_genes[!is.na(neighbor_genes) & neighbor_genes != ""]
  neighbor_genes <- setdiff(neighbor_genes, available_pathway_genes)
  neighbor_genes <- intersect(neighbor_genes, available_genes)
  
  # remove blacklisted genes if provided
  if(!is.null(genes_blacklist)){
    n_removed <- length(intersect(neighbor_genes, genes_blacklist))
    if(verbose) message(sprintf("  -> Removing %d blacklisted genes", n_removed))
    neighbor_genes <- setdiff(neighbor_genes, genes_blacklist)
  }
  
  if(verbose) message(sprintf("  -> Found %d neighbor genes in expression data",
                              length(neighbor_genes)))
  
  # rank neighbours by degree and select top N
  all_genes <- c(available_pathway_genes, neighbor_genes)
  
  if(verbose) message("Building interaction network (single query)...")
  
  all_mapped <- suppressWarnings(suppressMessages(
    string_db$map(data.frame(gene = all_genes), "gene", removeUnmappedRows = TRUE)
  ))
  all_ids <- all_mapped$STRING_id
  
  # cache this result - reused for the final graph
  network_edges_all <- string_db$get_interactions(all_ids)
  if(nrow(network_edges_all) == 0) stop("No network edges found...")
  
  degree_table   <- table(c(network_edges_all$from, network_edges_all$to))
  degree_symbols <- id_to_symbol[names(degree_table)]
  
  top_neighbors <- setdiff(degree_symbols[order(degree_table, decreasing = TRUE)],
                           available_pathway_genes)
  top_neighbors <- unique(top_neighbors[!is.na(top_neighbors) & top_neighbors != ""])
  top_neighbors <- intersect(top_neighbors, available_genes)
  top_neighbors <- head(top_neighbors, max_added_genes)
  
  expanded_genes <- unique(c(available_pathway_genes, top_neighbors))
  
  if(verbose){
    message(sprintf("  -> Selected top %d neighbors (total expanded: %d genes)",
                    length(top_neighbors), length(expanded_genes)))
  }
  
  # step 4: build network graph
  if(verbose) message("Building network graph...")
  
  expanded_ids <- all_mapped$STRING_id[all_mapped$gene %in% expanded_genes]
  
  network_edges <- network_edges_all[
    network_edges_all$from %in% expanded_ids &
      network_edges_all$to   %in% expanded_ids, ]
  
  if(nrow(network_edges) == 0) stop("No network edges found for expanded gene set...")
  
  g <- graph_from_data_frame(network_edges, directed = FALSE)
  
  # rename nodes using the pre-built STRING lookup
  V(g)$name <- id_to_symbol[V(g)$name]
  
  g <- induced_subgraph(g, which(V(g)$name %in% available_genes & !is.na(V(g)$name) & V(g)$name != ""))
  
  if(verbose){
    message(sprintf("  -> Network created: %d nodes, %d edges",
                    length(V(g)), length(E(g))))
  }
  
  # step 5: compute node metrics
  if(verbose) message("Computing network metrics...")
  
  gene_names <- V(g)$name
  
  origin <- dplyr::case_when(
    !is.null(seed_gene) & gene_names == seed_gene ~ "seed",
    gene_names %in% available_pathway_genes        ~ "pathway",
    gene_names %in% top_neighbors                  ~ "extended set",
    TRUE                                           ~ "unknown"
  )
  
  hits_result  <- hits_scores(g, scale = TRUE)
  eigen_scores <- eigen_centrality(g)$vector
  
  node_metrics <- data.frame(
    name        = gene_names,
    origin      = origin,
    degree      = igraph::degree(g),
    betweenness = betweenness(g),
    closeness   = closeness(g),
    eigenvector = eigen_scores,
    hub_score   = hits_result$hub,
    community   = as.numeric(membership(cluster_louvain(g)))
  )
  
  if(verbose){
    message(sprintf("  -> Computed metrics for %d nodes", nrow(node_metrics)))
    origin_counts <- table(node_metrics$origin)
    for(orig in names(origin_counts)){
      message(sprintf("    - %s: %d", orig, origin_counts[orig]))
    }
  }
  
  edge_data <- as_data_frame(g, what = "edges")
  
  return(list(
    edge_data    = edge_data,
    node_metrics = node_metrics
  ))
  
}