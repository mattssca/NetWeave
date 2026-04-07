#' Plot Subtype-Specific Hub Gene Scores as a Heatmap
#'
#' Visualises the subtype hub scores computed by
#' \code{calculate_subtype_hub_scores()} as a \code{ComplexHeatmap} heatmap,
#' with subtypes as rows and top-ranked genes as columns.
#'
#' @param node_metrics A data frame of network node metrics containing
#'   \code{<subtype>_hub_score} columns as produced by
#'   \code{calculate_subtype_hub_scores()}. Must also contain a \code{name}
#'   column. An \code{origin} column is optional but, if present, gene names
#'   are coloured by origin: \code{"seed"} = red, \code{"pathway"} = orange,
#'   \code{"extended set"} = teal.
#' @param top_n Integer. Number of genes to display. Which genes are selected
#'   is controlled by \code{select_by}. Defaults to \code{30}.
#' @param hub_direction Character string. One of \code{"up"} (default) to show
#'   the genes with the highest positive hub scores, or \code{"down"} to show
#'   genes with the most negative hub scores (downregulated high-centrality
#'   genes).
#' @param select_by Character string controlling which metric is used to pick
#'   the \code{top_n} genes. One of:
#'   \describe{
#'     \item{\code{"hub_score"}}{(Default) Select genes with the highest
#'       maximum hub score across subtypes — the most active overall hubs.}
#'     \item{\code{"delta"}}{Select genes with the largest difference between
#'       their top and second-highest absolute hub score — the most
#'       subtype-exclusive hubs. Requires \code{subtype_delta} to be present
#'       in \code{node_metrics} (computed by
#'       \code{calculate_subtype_hub_scores()}).}
#'   }
#' @param order_by Character string controlling column (gene) ordering after
#'   selection. One of:
#'   \describe{
#'     \item{\code{"hub_score"}}{(Default) Ordered by descending maximum hub
#'       score across subtypes.}
#'     \item{\code{"delta"}}{Ordered by descending subtype-specificity delta.
#'       Genes unique to one subtype appear leftmost. Ties broken by
#'       descending maximum hub score.}
#'     \item{\code{"alphabetical"}}{Alphabetical by gene name.}
#'     \item{\code{"none"}}{Hierarchical clustering.}
#'   }
#' @param color_low Character. Hex colour for the low end of the
#'   upregulated (\code{hub_direction = "up"}) colour scale.
#'   Defaults to \code{"#FFF5E6"}.
#' @param color_mid Character. Hex colour for the midpoint.
#'   Defaults to \code{"#E85D04"}.
#' @param color_high Character. Hex colour for the high end.
#'   Defaults to \code{"#370617"}.
#' @param title Character string. Heatmap column title. If \code{NULL}
#'   (default) a title is generated automatically.
#' @param fontsize_row Numeric. Font size for row (subtype) labels.
#'   Defaults to \code{10}.
#' @param fontsize_col Numeric. Font size for column (gene) labels.
#'   Defaults to \code{11}.
#' @param output_file Character string. Path to save the heatmap as a PDF.
#'   If \code{NULL} (default) the plot is drawn to the active device only.
#' @param width Numeric. Width of the saved PDF in inches. Defaults to \code{10}.
#' @param height Numeric. Height of the saved PDF in inches. Defaults to \code{8}.
#' @param verbose Logical. If \code{TRUE} (default), messages about the
#'   detected subtypes and selected genes are printed.
#'
#' @return A \code{Heatmap} object (invisibly). The heatmap is drawn to the
#'   active graphics device as a side effect.
#'
#' @seealso \code{\link{calculate_subtype_hub_scores}}
#'
#' @importFrom ComplexHeatmap Heatmap draw Legend
#' @importFrom circlize colorRamp2
#' @importFrom grid gpar
#' @importFrom tibble column_to_rownames
#'
#' @export
plot_subtype_hub_scores <- function(node_metrics = NULL,
                                    top_n = 30,
                                    hub_direction = c("up", "down"),
                                    select_by     = c("hub_score", "delta"),
                                    order_by      = c("hub_score", "delta", "alphabetical", "none"),
                                    color_low  = "#FFF5E6",
                                    color_mid  = "#E85D04",
                                    color_high = "#370617",
                                    title = NULL,
                                    fontsize_row = 10,
                                    fontsize_col = 11,
                                    output_file = NULL,
                                    width = 10,
                                    height = 8,
                                    verbose = TRUE){
  
  hub_direction <- match.arg(hub_direction)
  select_by     <- match.arg(select_by)
  order_by      <- match.arg(order_by)
  
  # checks
  if(is.null(node_metrics)){
    stop("User must provide an incoming object with network node metrics...")
  }
  
  if(!"name" %in% colnames(node_metrics)){
    stop("node_metrics must contain a 'name' column...")
  }
  
  if(!requireNamespace("ComplexHeatmap", quietly = TRUE)){
    stop("Package 'ComplexHeatmap' is required. Install with: BiocManager::install('ComplexHeatmap')")
  }
  
  if(!requireNamespace("circlize", quietly = TRUE)){
    stop("Package 'circlize' is required. Install with: install.packages('circlize')")
  }
  
  if(select_by == "delta" && !"subtype_delta" %in% colnames(node_metrics)){
    stop("select_by = 'delta' requires a 'subtype_delta' column. Run calculate_subtype_hub_scores() first...")
  }
  
  # auto-detect subtypes from <subtype>_hub_score columns
  hub_score_cols <- grep("_hub_score$", colnames(node_metrics), value = TRUE)
  subtypes       <- sub("_hub_score$", "", hub_score_cols)
  
  if(length(hub_score_cols) < 2){
    stop("At least two <subtype>_hub_score columns are required. Run calculate_subtype_hub_scores() first...")
  }
  
  if(verbose){
    message("Plotting subtype hub score heatmap...")
    message(sprintf("  -> Detected subtypes: %s", paste(subtypes, collapse = ", ")))
    message(sprintf("  -> Selecting top %d genes by: %s", top_n, select_by))
  }
  
  # always compute max/min hub score — used for selection and/or tie-breaking
  node_metrics$max_hub_score <- apply(node_metrics[, hub_score_cols, drop = FALSE], 1, max, na.rm = TRUE)
  node_metrics$min_hub_score <- apply(node_metrics[, hub_score_cols, drop = FALSE], 1, min, na.rm = TRUE)
  
  # select top_n genes
  if(select_by == "delta"){
    # select the most subtype-exclusive hubs regardless of absolute magnitude
    top_genes <- node_metrics[order(node_metrics$subtype_delta, decreasing = TRUE), ]
    # for "down", restrict candidate pool to genes with any negative hub score
    if(hub_direction == "down"){
      top_genes <- top_genes[top_genes$min_hub_score < 0, ]
    } else {
      top_genes <- top_genes[top_genes$max_hub_score > 0, ]
    }
    top_genes <- head(top_genes, top_n)
    
  } else {
    # default: select by highest magnitude hub score
    if(hub_direction == "up"){
      top_genes <- node_metrics[order(node_metrics$max_hub_score, decreasing = TRUE), ]
    } else {
      top_genes <- node_metrics[order(node_metrics$min_hub_score, decreasing = FALSE), ]
    }
    top_genes <- head(top_genes, top_n)
  }
  
  if(verbose){
    message(sprintf("  -> Selected %d genes", nrow(top_genes)))
  }
  
  # compute sort_delta on the selected genes for use in ordering
  score_mat <- top_genes[, hub_score_cols, drop = FALSE]
  
  if(hub_direction == "up"){
    top_genes$sort_delta <- apply(score_mat, 1, function(x) {
      s <- sort(x, decreasing = TRUE, na.last = TRUE)
      if(length(s) >= 2) s[1] - s[2] else s[1]
    })
  } else {
    top_genes$sort_delta <- apply(score_mat, 1, function(x) {
      s <- sort(x, decreasing = FALSE, na.last = TRUE)
      if(length(s) >= 2) s[2] - s[1] else abs(s[1])
    })
  }
  
  # build matrix — subtypes as rows, genes as columns
  hub_matrix           <- as.matrix(score_mat)
  rownames(hub_matrix) <- top_genes$name
  hub_matrix           <- t(hub_matrix)
  rownames(hub_matrix) <- subtypes
  
  # gene name colours by origin if column exists
  if("origin" %in% colnames(top_genes)){
    origin_colour_map <- c("seed" = "#CB0404", "pathway" = "#FF9F00", "extended set" = "#309898")
    gene_name_colours <- origin_colour_map[top_genes$origin]
    gene_name_colours[is.na(gene_name_colours)] <- "black"
    present_origins   <- intersect(names(origin_colour_map), unique(top_genes$origin))
    origin_legend     <- ComplexHeatmap::Legend(
      labels     = present_origins,
      legend_gp  = grid::gpar(fill = origin_colour_map[present_origins]),
      title      = "Gene Origin",
      type       = "points",
      pch        = NA_integer_,
      background = origin_colour_map[present_origins]
    )
  } else {
    gene_name_colours <- rep("black", ncol(hub_matrix))
    origin_legend     <- NULL
  }
  names(gene_name_colours) <- colnames(hub_matrix)
  
  # column ordering
  if(order_by == "alphabetical"){
    col_order         <- order(colnames(hub_matrix))
    cluster_cols      <- FALSE
    gene_name_colours <- gene_name_colours[colnames(hub_matrix)[col_order]]
    
  } else if(order_by == "hub_score"){
    if(hub_direction == "up"){
      col_order <- order(top_genes$max_hub_score, decreasing = TRUE)
    } else {
      col_order <- order(top_genes$min_hub_score, decreasing = FALSE)
    }
    cluster_cols      <- FALSE
    gene_name_colours <- gene_name_colours[colnames(hub_matrix)[col_order]]
    
  } else if(order_by == "delta"){
    if(hub_direction == "up"){
      col_order <- order(top_genes$sort_delta, top_genes$max_hub_score, decreasing = TRUE)
    } else {
      col_order <- order(top_genes$sort_delta, -top_genes$min_hub_score, decreasing = TRUE)
    }
    cluster_cols      <- FALSE
    gene_name_colours <- gene_name_colours[colnames(hub_matrix)[col_order]]
    if(verbose){
      message(sprintf("  -> Ordering by subtype delta (most unique gene: %s)",
                      colnames(hub_matrix)[col_order[1]]))
    }
    
  } else {
    col_order         <- NULL
    cluster_cols      <- TRUE
    gene_name_colours <- gene_name_colours[colnames(hub_matrix)]
  }
  
  # colour scale
  if(hub_direction == "up"){
    max_val  <- max(hub_matrix, na.rm = TRUE)
    col_fun  <- circlize::colorRamp2(
      c(0, max_val / 2, max_val),
      c(color_low, color_mid, color_high)
    )
    legend_param <- list(
      title     = "Hub Score",
      at        = c(0, max_val / 2, max_val),
      labels    = c("0", "Mid", "Max"),
      direction = "vertical"
    )
  } else {
    min_val  <- min(hub_matrix, na.rm = TRUE)
    mid_val  <- min_val / 2
    col_fun  <- circlize::colorRamp2(
      c(min_val, mid_val, 0),
      c("#213C51", "#6594B1", "#FFF5E6")
    )
    legend_param <- list(
      title     = "Hub Score",
      at        = c(0, mid_val, min_val),
      labels    = c("0", "Mid", "Min"),
      direction = "vertical"
    )
  }
  
  # auto title
  if(is.null(title)){
    title <- sprintf("Subtype-Specific Hub Gene Scores (%s, selected by %s)",
                     hub_direction, select_by)
  }
  
  # build heatmap
  ht <- ComplexHeatmap::Heatmap(
    hub_matrix,
    name                  = "Hub\nScore",
    col                   = col_fun,
    cluster_rows          = FALSE,
    cluster_columns       = cluster_cols,
    column_order          = col_order,
    show_row_names        = TRUE,
    show_column_names     = TRUE,
    rect_gp               = grid::gpar(col = "black", lwd = 2),
    row_names_side        = "left",
    row_names_gp          = grid::gpar(fontsize = fontsize_row),
    column_names_gp       = grid::gpar(fontsize = fontsize_col, col = gene_name_colours, fontface = "bold"),
    column_title          = title,
    heatmap_legend_param  = legend_param
  )
  
  # save to file if requested
  if(!is.null(output_file)){
    pdf(output_file, width = width, height = height)
    ComplexHeatmap::draw(ht, annotation_legend_list = if(!is.null(origin_legend)) list(origin_legend))
    dev.off()
    if(verbose) message(sprintf("  -> Heatmap saved to: %s", output_file))
  } else {
    ComplexHeatmap::draw(ht, annotation_legend_list = if(!is.null(origin_legend)) list(origin_legend))
  }
  
  invisible(ht)
  
}