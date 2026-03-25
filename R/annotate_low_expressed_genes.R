#' Annotate Low Expressed Genes in Network Node Metrics
#'
#' Computes the mean expression across all samples for each network gene and
#' annotates genes as low-expressed based on a percentile or absolute threshold.
#' Two columns are appended to \code{node_metrics}: \code{mean_expr_all_samples}
#' and \code{is_low_expression}.
#'
#' @param node_metrics A data frame of network node metrics. Must contain a
#'   \code{name} column with gene identifiers matching the rownames of
#'   \code{expr_data}.
#' @param expr_data A numeric matrix or data frame of expression values with genes
#'   as rows (rownames = gene IDs) and samples as columns (colnames = sample IDs).
#' @param low_expr_threshold Numeric threshold controlling low-expression flagging.
#'   Values between 0 and 1 (exclusive) are interpreted as a percentile cutoff
#'   (e.g., \code{0.25} flags the bottom 25\% of network genes). Values >= 1 are
#'   treated as an absolute expression cutoff. A value of 0 disables filtering.
#'   Defaults to \code{0.25}.
#' @param verbose Logical. If \code{TRUE} (default), progress messages are printed.
#'
#' @return The \code{node_metrics} data frame with two additional columns:
#'   \describe{
#'     \item{mean_expr_all_samples}{Mean expression across all samples for each gene.
#'       \code{NA} for genes not found in \code{expr_data}.}
#'     \item{is_low_expression}{Logical flag indicating whether the gene is below
#'       the expression threshold. \code{NA} for genes not found in \code{expr_data}.}
#'   }
#'
#' @details
#' Gene matching is performed by intersecting \code{node_metrics$name} with
#' \code{rownames(expr_data)}. The percentile threshold is computed from the
#' distribution of mean expression values across matched network genes only.
#' Genes at exactly the cutoff value are included in the low-expression set.
#'
#' @importFrom stats rowMeans quantile
#'
#' @export
#' 
annotate_low_expressed_genes = function(node_metrics = NULL, 
                                        expr_data = NULL, 
                                        low_expr_threshold = 0.25,
                                        verbose = TRUE){
  
  # checks
  if(is.null(node_metrics)){
    stop("User must provide an incoming object with network node metrics...")
  }
  
  if(!"name" %in% colnames(node_metrics)){
    stop("node_metrics must contain a 'name' column...")
  }
  
  if(is.null(expr_data)){
    stop("User must provide expression data...")
  }
  
  if(!is.matrix(expr_data) && !is.data.frame(expr_data)){
    stop("expr_data must be a matrix or data.frame...")
  }
  
  if(!is.numeric(low_expr_threshold) || length(low_expr_threshold) != 1 || low_expr_threshold < 0){
    stop("low_expr_threshold must be a single non-negative numeric value...")
  }
  
  if(verbose){
    message("Annotating low expressed genes...")
  }
  
  # calculate mean expression across all samples for network genes
  genes_in_data <- intersect(node_metrics$name, rownames(expr_data))
  
  if(length(genes_in_data) == 0){
    stop("No network gene names matched rownames of expr_data — check that rownames are gene identifiers")
  }
  
  gene_mean_expr <- rowMeans(expr_data[genes_in_data, , drop = FALSE], na.rm = TRUE)
  
  # determine threshold based on percentile or absolute value
  if(low_expr_threshold > 0 && low_expr_threshold < 1){
    # interpret as percentile (e.g., 0.25 = bottom 25%)
    expr_cutoff <- quantile(gene_mean_expr, probs = low_expr_threshold, na.rm = TRUE)
    threshold_type <- "percentile"
  } else if(low_expr_threshold >= 1){
    # interpret as absolute expression value
    expr_cutoff <- low_expr_threshold
    threshold_type <- "absolute"
  } else {
    # low_expr_threshold == 0: no filtering
    expr_cutoff <- -Inf
    threshold_type <- "none"
  }
  
  # annotate genes (vectorized; genes absent from expr_data get NA)
  mean_expr_matched <- gene_mean_expr[match(node_metrics$name, names(gene_mean_expr))]
  is_low_expr <- mean_expr_matched <= expr_cutoff
  
  # add to node_metrics
  node_metrics$mean_expr_all_samples <- mean_expr_matched
  node_metrics$is_low_expression <- is_low_expr
  
  if(verbose){
    if(threshold_type == "percentile"){
      message(sprintf("  -> Using %g%% percentile threshold", low_expr_threshold * 100))
      message(sprintf("  -> Expression cutoff: %.3f", expr_cutoff))
    } else if(threshold_type == "absolute"){
      message(sprintf("  -> Using absolute expression threshold: %.3f", expr_cutoff))
    } else {
      message("  -> No low expression filtering applied")
    }
    
    n_low <- sum(is_low_expr, na.rm = TRUE)
    message(sprintf("  -> Flagged %d genes (%.1f%%) as low expression",
                    n_low, 100 * n_low / nrow(node_metrics)))
  }
  
  return(node_metrics)
  
}