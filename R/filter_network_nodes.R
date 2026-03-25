#' Filter Network Nodes by Expression and Tissue Specificity Criteria
#'
#' Applies one or more filters to a node metrics data frame to retain only
#' genes meeting the specified criteria. Filters can be applied independently
#' and in any combination.
#'
#' @param node_metrics A data frame of network node metrics as produced by
#'   the annotation pipeline. Must contain a \code{name} column.
#' @param tissue_filter Character string or \code{NULL}. Controls which genes
#'   are retained based on HPA tissue expression evidence. One of:
#'   \describe{
#'     \item{\code{"normal"}}{Retain genes with \code{hpa_normal_expressed == TRUE}.}
#'     \item{\code{"cancer"}}{Retain genes with \code{hpa_cancer_expressed == TRUE}.}
#'     \item{\code{"either"}}{Retain genes expressed in normal tissue OR cancer.}
#'     \item{\code{"both"}}{Retain genes expressed in normal tissue AND cancer.}
#'     \item{\code{NULL}}{No tissue filter applied (default).}
#'   }
#'   Requires columns from \code{derive_hpa_metrics()}.
#' @param remove_low_expressed Logical. If \code{TRUE}, genes flagged as
#'   low-expressed (\code{is_low_expression == TRUE}) are removed. Requires
#'   the \code{is_low_expression} column from
#'   \code{annotate_low_expressed_genes()}. Defaults to \code{FALSE}.
#' @param filter_anova_sig Logical. If \code{TRUE}, only genes with a
#'   significant ANOVA result (\code{anova_sig == TRUE}) are retained.
#'   Requires the \code{anova_sig} column from
#'   \code{run_anova_subtype_expression()}. Defaults to \code{FALSE}.
#' @param verbose Logical. If \code{TRUE} (default), messages showing how
#'   many genes pass each filter step are printed.
#'
#' @return A filtered version of \code{node_metrics} containing only genes
#'   that pass all active filters.
#'
#' @details
#' Filters are applied sequentially: tissue specificity -> low expression ->
#' ANOVA significance. Genes with \code{NA} in a filter column are treated
#' as not meeting the criterion and are removed.
#'
#' @seealso \code{\link{derive_hpa_metrics}}, \code{\link{annotate_low_expressed_genes}},
#'   \code{\link{run_anova_subtype_expression}}
#'
#' @export
filter_network_nodes <- function(node_metrics = NULL,
                                 tissue_filter = NULL,
                                 remove_low_expressed = FALSE,
                                 filter_anova_sig = FALSE,
                                 verbose = TRUE){
  
  # checks
  if(is.null(node_metrics)){
    stop("User must provide an incoming object with network node metrics...")
  }
  
  if(!"name" %in% colnames(node_metrics)){
    stop("node_metrics must contain a 'name' column...")
  }
  
  valid_tissue_options <- c("normal", "cancer", "either", "both")
  if(!is.null(tissue_filter) && !tissue_filter %in% valid_tissue_options){
    stop(sprintf(
      "tissue_filter must be one of: %s, or NULL",
      paste(sprintf('"%s"', valid_tissue_options), collapse = ", ")
    ))
  }
  
  if(!is.logical(remove_low_expressed) || length(remove_low_expressed) != 1){
    stop("remove_low_expressed must be a single logical value...")
  }
  
  if(!is.logical(filter_anova_sig) || length(filter_anova_sig) != 1){
    stop("filter_anova_sig must be a single logical value...")
  }
  
  n_start <- nrow(node_metrics)
  
  if(verbose){
    message("Filtering network nodes...")
    message(sprintf("  -> Starting with %d genes", n_start))
  }
  
  # filter 1: tissue specificity
  if(!is.null(tissue_filter)){
    
    if(tissue_filter %in% c("normal", "either", "both") &&
       !"hpa_normal_expressed" %in% colnames(node_metrics)){
      stop("tissue_filter requires 'hpa_normal_expressed' column. Run derive_hpa_metrics() first...")
    }
    
    if(tissue_filter %in% c("cancer", "either", "both") &&
       !"hpa_cancer_expressed" %in% colnames(node_metrics)){
      stop("tissue_filter requires 'hpa_cancer_expressed' column. Run derive_hpa_metrics() first...")
    }
    
    n_before <- nrow(node_metrics)
    
    node_metrics <- switch(tissue_filter,
                           normal = node_metrics[node_metrics$hpa_normal_expressed %in% TRUE, ],
                           cancer = node_metrics[node_metrics$hpa_cancer_expressed %in% TRUE, ],
                           either = node_metrics[node_metrics$hpa_normal_expressed %in% TRUE |
                                                   node_metrics$hpa_cancer_expressed %in% TRUE, ],
                           both   = node_metrics[node_metrics$hpa_normal_expressed %in% TRUE &
                                                   node_metrics$hpa_cancer_expressed %in% TRUE, ]
    )
    
    if(verbose){
      message(sprintf("  -> Tissue filter ('%s'): removed %d genes, %d remaining",
                      tissue_filter, n_before - nrow(node_metrics), nrow(node_metrics)))
    }
  }
  
  # filter 2: remove low expressed genes
  if(remove_low_expressed){
    
    if(!"is_low_expression" %in% colnames(node_metrics)){
      stop("remove_low_expressed requires 'is_low_expression' column. Run annotate_low_expressed_genes() first...")
    }
    
    n_before <- nrow(node_metrics)
    node_metrics <- node_metrics[!node_metrics$is_low_expression %in% TRUE, ]
    
    if(verbose){
      message(sprintf("  -> Low expression filter: removed %d genes, %d remaining",
                      n_before - nrow(node_metrics), nrow(node_metrics)))
    }
  }
  
  # filter 3: ANOVA significance
  if(filter_anova_sig){
    
    if(!"anova_sig" %in% colnames(node_metrics)){
      stop("filter_anova_sig requires 'anova_sig' column. Run run_anova_subtype_expression() first...")
    }
    
    n_before <- nrow(node_metrics)
    node_metrics <- node_metrics[node_metrics$anova_sig %in% TRUE, ]
    
    if(verbose){
      message(sprintf("  -> ANOVA significance filter: removed %d genes, %d remaining",
                      n_before - nrow(node_metrics), nrow(node_metrics)))
    }
  }
  
  if(verbose){
    message(sprintf("  -> Final: %d of %d genes retained (%.1f%%)",
                    nrow(node_metrics), n_start,
                    100 * nrow(node_metrics) / n_start))
  }
  
  return(node_metrics)
  
}
