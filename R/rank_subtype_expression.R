#' Rank Subtype Mean Expression Levels Across Network Genes
#'
#' For each gene in \code{node_metrics}, ranks the mean expression values across
#' subtypes (rank 1 = highest expression). Rank columns named \code{rank_<subtype>}
#' are appended to \code{node_metrics}. Requires that \code{add_subtype_expression()}
#' has already been run so that \code{mean_expr_<subtype>} columns are present.
#'
#' @param node_metrics A data frame of network node metrics containing
#'   \code{mean_expr_<subtype>} columns. Must also contain a \code{name} column.
#' @param subtypes An optional character vector of subtype labels to rank. If
#'   \code{NULL} (default), subtypes are detected automatically from columns
#'   matching \code{mean_expr_*}, excluding the \code{mean_expr_all_samples}
#'   column added by \code{annotate_low_expressed_genes()}.
#' @param ties_method Character string specifying how tied expression values are
#'   handled. Passed to \code{\link[base]{rank}}. Defaults to \code{"average"}.
#' @param verbose Logical. If \code{TRUE} (default), progress messages are printed.
#'
#' @return The \code{node_metrics} data frame with one additional \code{rank_<subtype>}
#'   column per subtype. Rank 1 indicates the highest mean expression for that gene.
#'   \code{NA} expression values are excluded from ranking (\code{na.last = "keep"}).
#'
#' @details
#' When \code{subtypes = NULL}, subtype labels are inferred by stripping the
#' \code{mean_expr_} prefix from all matching column names. The \code{mean_expr_all_samples}
#' column produced by \code{annotate_low_expressed_genes()} is explicitly excluded
#' from auto-detection. The order of subtypes follows the column order in
#' \code{node_metrics}. Ranking is performed per gene (row-wise). Genes where all
#' subtype expression values are \code{NA} will have \code{NA} for all rank columns.
#'
#' @export
rank_subtype_expression = function(node_metrics = NULL,
                                   subtypes = NULL,
                                   ties_method = "average",
                                   verbose = TRUE){
  
  # checks
  if(is.null(node_metrics)){
    stop("User must provide an incoming object with network node metrics...")
  }
  
  if(!"name" %in% colnames(node_metrics)){
    stop("node_metrics must contain a 'name' column...")
  }
  
  if(!ties_method %in% c("average", "first", "last", "random", "max", "min")){
    stop("ties_method must be one of: 'average', 'first', 'last', 'random', 'max', 'min'...")
  }
  
  # auto-detect subtypes from mean_expr_* columns if not provided
  if(is.null(subtypes)){
    detected_cols <- grep("^mean_expr_", colnames(node_metrics), value = TRUE)
    detected_cols <- detected_cols[detected_cols != "mean_expr_all_samples"]
    if(length(detected_cols) == 0){
      stop("No 'mean_expr_*' columns found in node_metrics — run add_subtype_expression() first, or provide subtypes explicitly...")
    }
    subtypes <- sub("^mean_expr_", "", detected_cols)
    if(verbose){
      message(sprintf("  -> Auto-detected %d subtypes: %s",
                      length(subtypes), paste(subtypes, collapse = ", ")))
    }
  } else {
    if(!is.character(subtypes) || length(subtypes) == 0){
      stop("subtypes must be a non-empty character vector...")
    }
    missing_cols <- paste0("mean_expr_", subtypes)[!paste0("mean_expr_", subtypes) %in% colnames(node_metrics)]
    if(length(missing_cols) > 0){
      stop(sprintf("Missing expected expression column(s) in node_metrics: %s",
                   paste(missing_cols, collapse = ", ")))
    }
  }
  
  if(verbose){
    message("Ranking subtype expression levels...")
  }
  
  # extract expression matrix for the relevant columns (genes x subtypes)
  subtype_expr_cols <- paste0("mean_expr_", subtypes)
  expr_matrix <- as.matrix(node_metrics[, subtype_expr_cols, drop = FALSE])
  
  # apply row-wise ranking (rank 1 = highest expression)
  rank_matrix <- t(apply(expr_matrix, 1, function(x){
    rank(-x, ties.method = ties_method, na.last = "keep")
  }))
  
  # assign rank columns
  rank_col_names <- paste0("rank_", subtypes)
  colnames(rank_matrix) <- rank_col_names
  node_metrics[, rank_col_names] <- rank_matrix
  
  if(verbose){
    message(sprintf("  -> Ranked expression levels for %d subtypes", length(subtypes)))
    message(sprintf("  -> Rank columns added: %s", paste(rank_col_names, collapse = ", ")))
  }
  
  return(node_metrics)
  
}