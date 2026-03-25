#' Draw a Network Graph Within R
#'
#' Renders the filtered and annotated network as a static \code{ggraph} plot
#' or an interactive \code{visNetwork} widget. Nodes can be coloured and sized
#' by any column in \code{node_metrics}; seed genes and community membership
#' are highlighted by default.
#'
#' @param node_metrics A data frame of network node metrics. Must contain a
#'   \code{name} column with HGNC gene symbols.
#' @param edge_data A data frame of network edges from
#'   \code{create_and_expand_network()}. Must contain \code{from} and \code{to}
#'   columns. Edges where either endpoint is absent from \code{node_metrics}
#'   are removed automatically.
#' @param layout Character string. igraph layout algorithm to use for the
#'   static plot. Common options: \code{"fr"} (Fruchterman-Reingold, default),
#'   \code{"kk"} (Kamada-Kawai), \code{"circle"}, \code{"stress"},
#'   \code{"drl"}. Ignored when \code{interactive = TRUE}.
#' @param color_by Character string. Name of the column in \code{node_metrics}
#'   to use for node fill colour. Defaults to \code{"community"}. Use
#'   \code{"origin"} to colour by pathway/seed/extended set membership.
#' @param size_by Character string. Name of a numeric column in
#'   \code{node_metrics} to scale node size by. Defaults to \code{"degree"}.
#'   Pass \code{NULL} for uniform node sizes.
#' @param label_nodes Character string controlling which nodes are labelled.
#'   One of:
#'   \describe{
#'     \item{\code{"seed"}}{(Default) Label only seed genes
#'       (\code{origin == "seed"}).}
#'     \item{\code{"all"}}{Label every node.}
#'     \item{\code{"none"}}{No labels.}
#'   }
#'   Requires an \code{origin} column for \code{"seed"}.
#' @param top_label_n Integer. When \code{label_nodes = "none"} is not set and
#'   there are many nodes, label only the top \code{top_label_n} by
#'   \code{size_by}. Set to \code{Inf} to label all. Defaults to \code{20}.
#' @param edge_alpha Numeric in \code{[0, 1]}. Transparency of edges.
#'   Defaults to \code{0.3}.
#' @param node_size_range Numeric vector of length 2. Minimum and maximum node
#'   sizes for scaling. Defaults to \code{c(3, 12)}.
#' @param title Character string. Plot title. Defaults to \code{NULL} (no
#'   title).
#' @param interactive Logical. If \code{TRUE}, renders an interactive
#'   \code{visNetwork} widget instead of a static \code{ggraph} plot. Requires
#'   the \code{visNetwork} package. Defaults to \code{FALSE}.
#' @param output_file Character string. Path to save the static plot as a PDF
#'   or PNG (extension determines format). Ignored when
#'   \code{interactive = TRUE}. Defaults to \code{NULL} (no save).
#' @param width Numeric. Width of saved file in inches. Defaults to \code{10}.
#' @param height Numeric. Height of saved file in inches. Defaults to
#'   \code{10}.
#' @param verbose Logical. If \code{TRUE} (default), progress messages are
#'   printed.
#'
#' @return When \code{interactive = FALSE}, a \code{ggplot} object (invisibly).
#'   When \code{interactive = TRUE}, a \code{visNetwork} object.
#'
#' @importFrom igraph graph_from_data_frame set_vertex_attr V
#' @importFrom ggraph ggraph geom_edge_link geom_node_point geom_node_label
#'   geom_node_text create_layout
#' @importFrom ggplot2 aes scale_size_continuous scale_fill_discrete
#'   scale_fill_gradient theme_void theme labs ggsave guides guide_legend
#'
#' @export
plot_network <- function(node_metrics = NULL,
                         edge_data = NULL,
                         layout = "fr",
                         color_by = "community",
                         size_by = "degree",
                         label_nodes = c("seed", "all", "none"),
                         top_label_n = 20,
                         edge_alpha = 0.3,
                         node_size_range = c(3, 12),
                         title = NULL,
                         interactive = FALSE,
                         output_file = NULL,
                         width = 10,
                         height = 10,
                         verbose = TRUE){
  
  label_nodes <- match.arg(label_nodes)
  
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
  
  if(!is.null(color_by) && !color_by %in% colnames(node_metrics)){
    stop(sprintf("color_by column '%s' not found in node_metrics...", color_by))
  }
  
  if(!is.null(size_by) && !size_by %in% colnames(node_metrics)){
    stop(sprintf("size_by column '%s' not found in node_metrics...", size_by))
  }
  
  if(label_nodes == "seed" && !"origin" %in% colnames(node_metrics)){
    warning("label_nodes = 'seed' requires an 'origin' column. Falling back to label_nodes = 'none'.")
    label_nodes <- "none"
  }
  
  # handle duplicate gene rows
  if(anyDuplicated(node_metrics$name)){
    warning("Duplicate names found in node_metrics — keeping first occurrence per gene.")
    node_metrics <- node_metrics[!duplicated(node_metrics$name), ]
  }
  
  # filter edges to retained nodes
  retained_genes <- node_metrics$name
  filtered_edges <- edge_data[
    edge_data$from %in% retained_genes & edge_data$to %in% retained_genes,
  ]
  
  if(nrow(filtered_edges) == 0){
    stop("No edges remain after filtering to retained nodes...")
  }
  
  if(verbose){
    message("Drawing network...")
    message(sprintf("  -> Nodes: %d, Edges: %d", length(retained_genes), nrow(filtered_edges)))
    message(sprintf("  -> Colour by: %s | Size by: %s", color_by, size_by %||% "uniform"))
  }
  
  # build igraph and attach all node attributes
  g             <- igraph::graph_from_data_frame(filtered_edges, directed = FALSE)
  missing_nodes <- setdiff(retained_genes, V(g)$name)
  if(length(missing_nodes) > 0){
    g <- g + igraph::vertices(missing_nodes)
  }
  
  attr_cols <- setdiff(colnames(node_metrics), "name")
  for(col in attr_cols){
    vals <- node_metrics[[col]][match(V(g)$name, node_metrics$name)]
    g    <- igraph::set_vertex_attr(g, name = col, value = vals)
  }
  
  # --- interactive plot ---
  if(interactive){
    
    if(!requireNamespace("visNetwork", quietly = TRUE)){
      stop("Package 'visNetwork' is required for interactive plots. Install with: install.packages('visNetwork')")
    }
    
    vis_nodes <- data.frame(
      id    = node_metrics$name,
      label = node_metrics$name,
      stringsAsFactors = FALSE
    )
    
    if(!is.null(size_by)){
      size_vals          <- node_metrics[[size_by]]
      size_vals          <- 10 + 30 * (size_vals - min(size_vals, na.rm = TRUE)) /
        (max(size_vals, na.rm = TRUE) - min(size_vals, na.rm = TRUE) + 0.001)
      vis_nodes$value    <- size_vals
    }
    
    if(!is.null(color_by)){
      color_col <- node_metrics[[color_by]]
      if(is.numeric(color_col)){
        # assign palette colours by group
        groups           <- as.character(color_col)
        vis_nodes$group  <- groups
      } else {
        vis_nodes$group  <- as.character(color_col)
      }
    }
    
    if("origin" %in% colnames(node_metrics)){
      vis_nodes$title <- paste0("<b>", node_metrics$name, "</b><br>Origin: ",
                                node_metrics$origin)
    } else {
      vis_nodes$title <- paste0("<b>", node_metrics$name, "</b>")
    }
    
    vis_edges <- data.frame(
      from = filtered_edges$from,
      to   = filtered_edges$to,
      stringsAsFactors = FALSE
    )
    
    vn <- visNetwork::visNetwork(vis_nodes, vis_edges, main = title) %>%
      visNetwork::visOptions(
        highlightNearest = list(enabled = TRUE, degree = 1, hover = TRUE),
        nodesIdSelection = TRUE
      ) %>%
      visNetwork::visPhysics(
        solver = "forceAtlas2Based",
        forceAtlas2Based = list(gravitationalConstant = -50)
      ) %>%
      visNetwork::visEdges(color = list(color = "#cccccc", highlight = "#e85d04")) %>%
      visNetwork::visLayout(randomSeed = 42)
    
    if(verbose) message("  -> Interactive visNetwork widget created.")
    return(vn)
  }
  
  # --- static ggraph plot ---
  if(!requireNamespace("ggraph", quietly = TRUE)){
    stop("Package 'ggraph' is required for static plots. Install with: install.packages('ggraph')")
  }
  
  if(!requireNamespace("ggplot2", quietly = TRUE)){
    stop("Package 'ggplot2' is required. Install with: install.packages('ggplot2')")
  }
  
  # determine which nodes to label
  node_df         <- as.data.frame(igraph::as_data_frame(g, what = "vertices"))
  node_df$name    <- rownames(node_df)  # ggraph exposes vertex names as 'name'
  
  if(label_nodes == "seed"){
    node_df$show_label <- !is.na(node_df$origin) & node_df$origin == "seed"
  } else if(label_nodes == "all"){
    if(!is.null(size_by) && is.finite(top_label_n)){
      size_vals          <- node_df[[size_by]]
      rank_cutoff        <- sort(size_vals, decreasing = TRUE)[min(top_label_n, length(size_vals))]
      node_df$show_label <- !is.na(size_vals) & size_vals >= rank_cutoff
    } else {
      node_df$show_label <- TRUE
    }
  } else {
    node_df$show_label <- FALSE
  }
  
  # detect whether color_by is continuous or categorical
  color_vals      <- igraph::vertex_attr(g, color_by)
  color_is_numeric <- is.numeric(color_vals)
  
  # coerce to factor only for categorical (non-continuous) numeric columns
  # like community IDs — identified as integer with few unique values
  if(is.integer(color_vals) && length(unique(color_vals)) <= 20){
    g <- igraph::set_vertex_attr(g, color_by, value = factor(color_vals))
    color_is_numeric <- FALSE
  }
  
  p <- ggraph::ggraph(g, layout = layout) +
    ggraph::geom_edge_link(alpha = edge_alpha, colour = "#aaaaaa") +
    {
      if(!is.null(size_by)){
        ggraph::geom_node_point(
          ggplot2::aes(fill = .data[[color_by]], size = .data[[size_by]]),
          shape = 21, colour = "white", stroke = 0.4
        )
      } else {
        ggraph::geom_node_point(
          ggplot2::aes(fill = .data[[color_by]]),
          shape = 21, colour = "white", stroke = 0.4, size = 5
        )
      }
    } +
    ggraph::geom_node_text(
      ggplot2::aes(label = ifelse(node_df$show_label, name, "")),
      repel = TRUE, size = 3, colour = "black"
    ) +
    {
      if(!is.null(size_by)){
        ggplot2::scale_size_continuous(
          range  = node_size_range,
          name   = size_by,
          guide  = ggplot2::guide_legend(override.aes = list(fill = "grey50"))
        )
      }
    } +
    {
      if(color_is_numeric){
        mid_val <- mean(color_vals, na.rm = TRUE)
        ggplot2::scale_fill_gradient2(
          name  = color_by,
          low   = "#3B82F6",
          mid   = "#F9FAFB",
          high  = "#E85D04",
          midpoint = mid_val,
          na.value = "grey80"
        )
      } else {
        ggplot2::scale_fill_discrete(name = color_by)
      }
    } +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(legend.position = "right") +
    ggplot2::labs(title = title)
  
  if(!is.null(output_file)){
    ggplot2::ggsave(output_file, plot = p, width = width, height = height)
    if(verbose) message(sprintf("  -> Plot saved to: %s", output_file))
  } else {
    print(p)
  }
  
  if(verbose) message("  -> Done.")
  
  invisible(p)
  
}

# internal null-coalescing helper (not exported)
`%||%` <- function(a, b) if(!is.null(a)) a else b