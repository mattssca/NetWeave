#' Calculate Subtype-Specific Hub Scores
#'
#' Computes a hub score for each gene in each subtype by combining weighted
#' network centrality metrics with a subtype specificity score derived from
#' the mean expression columns already present in \code{node_metrics}. Subtypes
#' are detected automatically from \code{mean_expr_*} columns.
#'
#' @param node_metrics A data frame of network node metrics. Must contain a
#'   \code{name} column, centrality columns (\code{degree}, \code{betweenness},
#'   \code{eigenvector}), and at least two \code{mean_expr_<subtype>} columns
#'   as produced by \code{add_subtype_expression()}.
#' @param degree_weight Numeric. Weight applied to normalised degree centrality.
#'   Defaults to \code{0.4}.
#' @param betweenness_weight Numeric. Weight applied to normalised betweenness
#'   centrality. Defaults to \code{0.3}.
#' @param eigenvector_weight Numeric. Weight applied to eigenvector centrality.
#'   Defaults to \code{0.3}.
#'   The three weights are normalised to sum to 1 if they do not already.
#' @param specificity_method Character string controlling how subtype specificity
#'   is scored. One of:
#'   \describe{
#'     \item{\code{"expression_diff"}}{(Default) Difference between a gene's
#'       expression in the target subtype and the mean across all other subtypes,
#'       scaled by the gene's cross-subtype standard deviation.}
#'     \item{\code{"fold_change"}}{Log2 fold-change of target subtype expression
#'       over the mean of all other subtypes (pseudocount of 1 added).}
#'     \item{\code{"zscore"}}{Z-score of the target subtype expression relative
#'       to the gene's mean and standard deviation across all subtypes.}
#'   }
#'   In all methods the raw specificity score is linearly scaled to \code{[-1, 1]}
#'   across all genes before being multiplied by the centrality score.
#' @param verbose Logical. If \code{TRUE} (default), progress messages including
#'   the detected subtypes and per-subtype hub score ranges are printed.
#'
#' @return The \code{node_metrics} data frame with two additional columns per
#'   detected subtype:
#'   \describe{
#'     \item{\code{<subtype>_spec_score}}{Specificity score scaled to
#'       \code{[-1, 1]}. Positive values indicate higher expression in this
#'       subtype relative to others.}
#'     \item{\code{<subtype>_hub_score}}{Final hub score: weighted centrality
#'       multiplied by specificity score. Positive values indicate high-centrality
#'       genes upregulated in this subtype; negative values indicate
#'       high-centrality genes downregulated in this subtype.}
#'   }
#'
#' @details
#' \code{mean_expr_all_samples} is explicitly excluded from subtype detection.
#' Centrality metrics are min-max normalised to \code{[0, 1]} internally and
#' not written back to \code{node_metrics}. Genes with zero cross-subtype
#' standard deviation receive a floor of \code{0.01} to avoid division by zero.
#'
#' @seealso \code{\link{add_subtype_expression}}, \code{\link{recalculate_network_metrics}}
#'
#' @export
calculate_subtype_hub_scores <- function(node_metrics = NULL,
                                         degree_weight = 0.4,
                                         betweenness_weight = 0.3,
                                         eigenvector_weight = 0.3,
                                         specificity_method = c("expression_diff",
                                                                "fold_change",
                                                                "zscore"),
                                         verbose = TRUE){
  
  specificity_method <- match.arg(specificity_method)
  
  # checks
  if(is.null(node_metrics)){
    stop("User must provide an incoming object with network node metrics...")
  }
  
  if(!"name" %in% colnames(node_metrics)){
    stop("node_metrics must contain a 'name' column...")
  }
  
  required_centrality <- c("degree", "betweenness", "eigenvector")
  missing_centrality  <- required_centrality[!required_centrality %in% colnames(node_metrics)]
  if(length(missing_centrality) > 0){
    stop(sprintf(
      "node_metrics is missing centrality column(s): %s",
      paste(missing_centrality, collapse = ", ")
    ))
  }
  
  # auto-detect subtypes from mean_expr_* columns
  expr_cols <- grep("^mean_expr_", colnames(node_metrics), value = TRUE)
  expr_cols <- expr_cols[expr_cols != "mean_expr_all_samples"]
  subtypes  <- sub("^mean_expr_", "", expr_cols)
  
  if(length(expr_cols) < 2){
    stop("At least two mean_expr_<subtype> columns are required. Run add_subtype_expression() first...")
  }
  
  if(verbose){
    message("Calculating subtype-specific hub scores...")
    message(sprintf("  -> Detected subtypes: %s", paste(subtypes, collapse = ", ")))
    message(sprintf("  -> Specificity method: %s", specificity_method))
  }
  
  # normalise weights to sum to 1
  weight_sum <- degree_weight + betweenness_weight + eigenvector_weight
  if(abs(weight_sum - 1) > 0.01){
    warning(sprintf(
      "Weights sum to %.3f, not 1. Normalising automatically.", weight_sum
    ))
    degree_weight      <- degree_weight      / weight_sum
    betweenness_weight <- betweenness_weight / weight_sum
    eigenvector_weight <- eigenvector_weight / weight_sum
  }
  
  if(verbose){
    message(sprintf("  -> Weights — degree: %.2f, betweenness: %.2f, eigenvector: %.2f",
                    degree_weight, betweenness_weight, eigenvector_weight))
  }
  
  # helper: min-max normalise to [0, 1]
  norm_01 <- function(x){
    rng <- range(x, na.rm = TRUE)
    if(rng[2] == rng[1]) return(rep(0, length(x)))
    (x - rng[1]) / (rng[2] - rng[1])
  }
  
  # helper: scale vector to [-1, 1]
  norm_bipolar <- function(x){
    rng <- range(x, na.rm = TRUE)
    (x - rng[1]) / (rng[2] - rng[1] + 0.001) * 2 - 1
  }
  
  # compute weighted centrality score (shared across all subtypes)
  centrality_score <-
    degree_weight      * norm_01(node_metrics$degree) +
    betweenness_weight * norm_01(node_metrics$betweenness) +
    eigenvector_weight * norm_01(node_metrics$eigenvector)
  
  node_metrics[["centrality_score"]] <- centrality_score
  
  # compute global cross-subtype mean and sd per gene
  expr_matrix      <- as.matrix(node_metrics[, expr_cols, drop = FALSE])
  global_mean_expr <- rowMeans(expr_matrix, na.rm = TRUE)
  expr_sd          <- apply(expr_matrix, 1, sd, na.rm = TRUE)
  expr_sd[is.na(expr_sd) | expr_sd == 0] <- 0.01
  
  # compute per-subtype specificity and hub scores
  for(i in seq_along(subtypes)){
    
    subtype    <- subtypes[i]
    target_col <- expr_cols[i]
    other_cols <- expr_cols[-i]
    
    if(specificity_method == "expression_diff"){
      other_mean <- rowMeans(node_metrics[, other_cols, drop = FALSE], na.rm = TRUE)
      spec_raw   <- (node_metrics[[target_col]] - other_mean) / expr_sd
      
    } else if(specificity_method == "fold_change"){
      other_mean <- rowMeans(node_metrics[, other_cols, drop = FALSE], na.rm = TRUE)
      spec_raw   <- log2((node_metrics[[target_col]] + 1) / (other_mean + 1))
      
    } else if(specificity_method == "zscore"){
      spec_raw   <- (node_metrics[[target_col]] - global_mean_expr) / expr_sd
    }
    
    spec_score <- norm_bipolar(spec_raw)
    
    node_metrics[[paste0(subtype, "_spec_score")]] <- spec_score
    node_metrics[[paste0(subtype, "_hub_score")]]  <- centrality_score * spec_score
  }
  
  if(verbose){
    for(subtype in subtypes){
      col <- paste0(subtype, "_hub_score")
      message(sprintf("  -> %s hub score range: [%.3f, %.3f]",
                      subtype,
                      min(node_metrics[[col]], na.rm = TRUE),
                      max(node_metrics[[col]], na.rm = TRUE)))
    }
  }
  
  return(node_metrics)
  
}
