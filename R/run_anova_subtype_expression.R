#' Run ANOVA Analysis for Differential Expression Across Subtypes
#'
#' Performs a one-way ANOVA for each network gene to test whether mean expression
#' differs significantly across subtypes. Results including p-values, Bonferroni-
#' corrected p-values, eta-squared effect sizes, and a significance flag are
#' appended as new columns to \code{node_metrics}.
#'
#' @param node_metrics A data frame of network node metrics. Must contain a
#'   \code{name} column with gene identifiers matching the rownames of
#'   \code{expr_data}.
#' @param expr_data A numeric matrix or data frame of expression values with genes
#'   as rows (rownames = gene IDs) and samples as columns (colnames = sample IDs).
#' @param this_subtype_vector A named character vector mapping sample IDs (names)
#'   to their subtype labels (values). Must be named.
#' @param sig_threshold Numeric significance threshold applied to Bonferroni-
#'   corrected p-values for the \code{anova_sig} flag. Defaults to \code{0.05}.
#' @param verbose Logical. If \code{TRUE} (default), progress messages are printed.
#'
#' @return The \code{node_metrics} data frame with four additional columns:
#'   \describe{
#'     \item{anova_p_value}{Raw p-value from the one-way ANOVA.}
#'     \item{anova_bonf}{Bonferroni-corrected p-value (adjusted over all network genes).}
#'     \item{anova_effect_size}{Eta-squared effect size (SS_between / SS_total).}
#'     \item{anova_sig}{Logical flag; \code{TRUE} if \code{anova_bonf < sig_threshold}.}
#'   }
#'   Genes not present in \code{expr_data} or with insufficient subtype
#'   representation will have \code{NA} for all four columns.
#'
#' @details
#' Only samples present in both \code{colnames(expr_data)} and
#' \code{names(this_subtype_vector)} are used. For each gene, ANOVA is only run
#' when at least two subtypes have non-missing expression values. Bonferroni
#' correction is applied over all network genes (including untested ones), giving
#' a conservative correction.
#'
#' @importFrom stats aov p.adjust complete.cases
#'
#' @export
#' 
run_anova_subtype_expression = function(node_metrics = NULL,
                                        expr_data = NULL,
                                        this_subtype_vector = NULL,
                                        sig_threshold = 0.05,
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
  
  if(is.null(this_subtype_vector)){
    stop("User must provide a named vector with sample/subtype information...")
  }
  
  if(is.null(names(this_subtype_vector))){
    stop("this_subtype_vector must be a named vector where names are sample IDs...")
  }
  
  if(!is.numeric(sig_threshold) || length(sig_threshold) != 1 ||
     sig_threshold <= 0 || sig_threshold >= 1){
    stop("sig_threshold must be a single numeric value between 0 and 1...")
  }
  
  if(verbose){
    message("Running ANOVA analysis across subtypes...")
  }
  
  # identify testable genes and valid samples
  genes_in_data <- intersect(node_metrics$name, rownames(expr_data))
  valid_samples <- intersect(colnames(expr_data), names(this_subtype_vector))
  
  if(length(genes_in_data) == 0){
    stop("No network gene names matched rownames of expr_data — check that rownames are gene identifiers")
  }
  
  if(length(valid_samples) == 0){
    stop("No samples found in both expr_data and this_subtype_vector — check that sample IDs match")
  }
  
  # initialize result vectors
  anova_p_values  <- rep(NA_real_, nrow(node_metrics))
  anova_effect    <- rep(NA_real_, nrow(node_metrics))
  
  # run ANOVA for each gene
  for(i in seq_len(nrow(node_metrics))){
    gene_name <- node_metrics$name[i]
    
    if(!gene_name %in% genes_in_data) next
    
    gene_expr <- expr_data[gene_name, valid_samples]
    gene_subtypes <- this_subtype_vector[valid_samples]
    
    anova_df <- data.frame(
      expression = as.numeric(gene_expr),
      subtype    = factor(gene_subtypes)
    )
    anova_df <- anova_df[complete.cases(anova_df), ]
    
    if(length(unique(anova_df$subtype)) < 2 || nrow(anova_df) < 2) next
    
    anova_result <- tryCatch({
      aov(expression ~ subtype, data = anova_df)
    }, error = function(e) NULL)
    
    if(is.null(anova_result)) next
    
    anova_summary <- summary(anova_result)
    anova_p_values[i] <- anova_summary[[1]]["subtype", "Pr(>F)"]
    
    ss_between <- anova_summary[[1]]["subtype", "Sum Sq"]
    ss_total   <- sum(anova_summary[[1]][, "Sum Sq"])
    anova_effect[i] <- ss_between / ss_total
  }
  
  # apply Bonferroni correction over all network genes
  anova_bonf <- p.adjust(anova_p_values, method = "bonferroni")
  
  # significance flag
  anova_sig <- anova_bonf < sig_threshold
  anova_sig[is.na(anova_bonf)] <- NA
  
  # add columns to node_metrics
  node_metrics$anova_p_value    <- anova_p_values
  node_metrics$anova_bonf       <- anova_bonf
  node_metrics$anova_effect_size <- anova_effect
  node_metrics$anova_sig        <- anova_sig
  
  if(verbose){
    n_tested <- sum(!is.na(anova_p_values))
    n_sig    <- sum(anova_sig, na.rm = TRUE)
    
    message(sprintf("  -> Valid samples for analysis: %d", length(valid_samples)))
    message(sprintf("  -> Tested %d genes", n_tested))
    
    if(n_tested > 0){
      message(sprintf("  -> Significant genes (Bonferroni p < %g): %d (%.1f%%)",
                      sig_threshold, n_sig, 100 * n_sig / n_tested))
      
      if(n_sig > 0){
        effect_range <- range(anova_effect[anova_sig], na.rm = TRUE)
        message(sprintf("  -> Effect size range for significant genes: %.3f - %.3f",
                        effect_range[1], effect_range[2]))
      }
    } else {
      warning("No genes were tested — check sample/subtype matching")
    }
  }
  
  return(node_metrics)
  
}