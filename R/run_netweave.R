#' Run the Full NetWeave Analysis Pipeline
#'
#' A convenience wrapper that executes all eleven pipeline steps in sequence:
#' pathway gene retrieval, network construction and expansion, subtype expression
#' annotation, low-expression flagging, ANOVA, subtype ranking, HPA annotation,
#' HPA metric derivation, node filtering, network metric recalculation, and
#' subtype hub score calculation.
#'
#' @param seed_gene A single character string naming the seed gene of interest
#'   (e.g. \code{"ERBB2"}).
#' @param pathway_name A character string passed to \code{get_pathway_genes()}
#'   as \code{this_pathway}. Can be a Reactome pathway name or gene symbol.
#'   Ignored when \code{gene_list} is supplied.
#' @param gene_list An optional character vector of HGNC gene symbols. When
#'   provided, the Reactome lookup in step 1 is bypassed and this vector is
#'   used as the seed gene set directly. Takes precedence over
#'   \code{pathway_name}.
#' @param expr_data A numeric matrix of expression values (genes x samples) used
#'   for network construction, subtype expression, and ANOVA steps.
#' @param subtype_vector A named character vector mapping sample IDs to subtype
#'   labels. Used for subtype expression and ANOVA steps.
#' @param normal_hpa HPA normal tissue data frame. Passed to
#'   \code{add_hpa_annotations()}.
#' @param cancer_hpa HPA cancer data frame. Passed to
#'   \code{add_hpa_annotations()}.
#' @param this_tissue_normal Character string specifying the normal tissue to
#'   filter in HPA (e.g. \code{"Urinary bladder"}).
#' @param this_tissue_cancer Character string specifying the cancer type to
#'   filter in HPA (e.g. \code{"urothelial cancer"}).
#' @param low_expr_data A numeric matrix of expression values used specifically
#'   for the low-expression annotation step. Defaults to \code{expr_data} when
#'   \code{NULL}.
#' @param subtype_order An optional character vector specifying the desired
#'   subtype order. Passed to \code{add_subtype_expression()}. Defaults to
#'   \code{NULL}.
#' @param max_added_genes Integer. Maximum number of neighbour genes to add
#'   during network expansion. Defaults to \code{20}.
#' @param string_score_threshold Integer. Minimum STRING combined score (0-1000)
#'   to retain an edge. Defaults to \code{900}.
#' @param string_data_dir Character string path to a directory containing
#'   pre-downloaded STRING flat files. Defaults to \code{"data/string_db"}.
#' @param genes_blacklist Optional character vector of gene symbols to exclude
#'   from neighbour expansion. Defaults to \code{NULL}.
#' @param apply_blacklist Boolean flag, set to TRUE for removing ubiquitous cellular machinery 
#' (ribosomal proteins, histones, splicing factors, core translation factors). Default is FALSE.
#' @param low_expr_threshold Numeric threshold for low-expression annotation.
#'   Values in (0, 1) are treated as a percentile; values >= 1 as an absolute
#'   cutoff; \code{0} disables filtering. Defaults to \code{0.25}.
#' @param sig_threshold Numeric significance threshold for Bonferroni-corrected
#'   ANOVA p-values. Defaults to \code{0.05}.
#' @param tissue_filter Character string controlling HPA tissue filter in
#'   \code{filter_network_nodes()}. One of \code{"normal"}, \code{"cancer"},
#'   \code{"either"}, \code{"both"}, or \code{NULL}. Defaults to
#'   \code{"either"}.
#' @param remove_low_expressed Logical. If \code{TRUE}, low-expressed genes are
#'   removed during filtering. Defaults to \code{TRUE}.
#' @param filter_anova_sig Logical. If \code{TRUE}, only ANOVA-significant genes
#'   are retained during filtering. Defaults to \code{TRUE}.
#' @param degree_weight Numeric weight for degree centrality in hub score
#'   calculation. Defaults to \code{0.4}.
#' @param betweenness_weight Numeric weight for betweenness centrality. Defaults
#'   to \code{0.3}.
#' @param eigenvector_weight Numeric weight for eigenvector centrality. Defaults
#'   to \code{0.3}.
#' @param specificity_method Character string controlling subtype specificity
#'   scoring. One of \code{"expression_diff"}, \code{"fold_change"}, or
#'   \code{"zscore"}. Defaults to \code{"expression_diff"}.
#' @param verbose Logical. If \code{TRUE} (default), progress messages are
#'   printed at each step.
#'
#' @return A named list with three elements:
#'   \describe{
#'     \item{node_metrics}{Final annotated and filtered node metrics data frame
#'       with all pipeline columns appended.}
#'     \item{edge_data}{Edge data frame from \code{create_and_expand_network()},
#'       representing the full (pre-filter) network edges.}
#'     \item{pathway_genes}{Character vector of pathway genes retrieved in
#'       step 1.}
#'   }
#'
#' @seealso
#'   \code{\link{get_pathway_genes}},
#'   \code{\link{create_and_expand_network}},
#'   \code{\link{add_subtype_expression}},
#'   \code{\link{annotate_low_expressed_genes}},
#'   \code{\link{run_anova_subtype_expression}},
#'   \code{\link{rank_subtype_expression}},
#'   \code{\link{add_hpa_annotations}},
#'   \code{\link{derive_hpa_metrics}},
#'   \code{\link{filter_network_nodes}},
#'   \code{\link{recalculate_network_metrics}},
#'   \code{\link{calculate_subtype_hub_scores}}
#'
#' @export
run_netweave <- function(seed_gene,
                         pathway_name       = NULL,
                         gene_list          = NULL,
                         expr_data,
                         subtype_vector,
                         normal_hpa,
                         cancer_hpa,
                         this_tissue_normal,
                         this_tissue_cancer,
                         low_expr_data        = NULL,
                         subtype_order        = NULL,
                         max_added_genes      = 20,
                         string_score_threshold = 900,
                         string_data_dir      = "data/string_db",
                         genes_blacklist      = NULL,
                         apply_blacklist     = FALSE,
                         low_expr_threshold   = 0.25,
                         sig_threshold        = 0.05,
                         tissue_filter        = "either",
                         remove_low_expressed = TRUE,
                         filter_anova_sig     = TRUE,
                         degree_weight        = 0.4,
                         betweenness_weight   = 0.3,
                         eigenvector_weight   = 0.3,
                         specificity_method   = "expression_diff",
                         verbose              = TRUE) {
  
  if (is.null(low_expr_data)) low_expr_data <- expr_data
  
  # step 1 - get pathway genes
  if (verbose) message("=== Step 1: Retrieving pathway genes ===")
  pathway_genes <- get_pathway_genes(this_pathway = pathway_name,
                                     gene_list    = gene_list,
                                     verbose      = verbose)
  
  # step 2 - create and expand network
  if (verbose) message("=== Step 2: Building network ===")
  network <- create_and_expand_network(expr_data              = expr_data,
                                       seed_gene              = seed_gene,
                                       pathway_genes          = pathway_genes,
                                       max_added_genes        = max_added_genes,
                                       string_score_threshold = string_score_threshold,
                                       verbose                = verbose,
                                       genes_blacklist        = genes_blacklist,
                                       apply_blacklist        = apply_blacklist,
                                       string_data_dir        = string_data_dir)
  
  node_metrics <- network$node_metrics
  
  # step 3 - add subtype mean expressions
  if (verbose) message("=== Step 3: Adding subtype expression ===")
  node_metrics <- add_subtype_expression(this_subtype_vector = subtype_vector,
                                         expr_data           = expr_data,
                                         node_metrics        = node_metrics,
                                         subtype_order       = subtype_order,
                                         verbose             = verbose)
  
  # step 4 - annotate low expressed genes
  if (verbose) message("=== Step 4: Annotating low-expressed genes ===")
  node_metrics <- annotate_low_expressed_genes(node_metrics       = node_metrics,
                                               expr_data          = low_expr_data,
                                               low_expr_threshold = low_expr_threshold,
                                               verbose            = verbose)
  
  # step 5 - ANOVA for subtype-specific expression
  if (verbose) message("=== Step 5: Running ANOVA ===")
  node_metrics <- run_anova_subtype_expression(node_metrics        = node_metrics,
                                               expr_data           = expr_data,
                                               this_subtype_vector = subtype_vector,
                                               sig_threshold       = sig_threshold,
                                               verbose             = verbose)
  
  # step 6 - rank subtype expression
  if (verbose) message("=== Step 6: Ranking subtype expression ===")
  node_metrics <- rank_subtype_expression(node_metrics = node_metrics,
                                          verbose      = verbose)
  
  # step 7 - add HPA annotations
  if (verbose) message("=== Step 7: Adding HPA annotations ===")
  node_metrics <- add_hpa_annotations(node_metrics        = node_metrics,
                                      normal_hpa          = normal_hpa,
                                      cancer_hpa          = cancer_hpa,
                                      this_tissue_normal  = this_tissue_normal,
                                      this_tissue_cancer  = this_tissue_cancer,
                                      verbose             = verbose)
  
  # step 8 - derive HPA metrics
  if (verbose) message("=== Step 8: Deriving HPA metrics ===")
  node_metrics <- derive_hpa_metrics(node_metrics = node_metrics,
                                     verbose      = verbose)
  
  # step 9 - filter nodes
  if (verbose) message("=== Step 9: Filtering network nodes ===")
  node_metrics <- filter_network_nodes(node_metrics         = node_metrics,
                                       tissue_filter        = tissue_filter,
                                       remove_low_expressed = remove_low_expressed,
                                       filter_anova_sig     = filter_anova_sig,
                                       verbose              = verbose)
  
  # step 10 - recalculate network metrics on filtered set
  if (verbose) message("=== Step 10: Recalculating network metrics ===")
  node_metrics <- recalculate_network_metrics(node_metrics = node_metrics,
                                              edge_data    = network$edge_data,
                                              verbose      = verbose)
  
  # step 11 - calculate subtype hub scores
  if (verbose) message("=== Step 11: Calculating subtype hub scores ===")
  node_metrics <- calculate_subtype_hub_scores(node_metrics       = node_metrics,
                                               degree_weight      = degree_weight,
                                               betweenness_weight = betweenness_weight,
                                               eigenvector_weight = eigenvector_weight,
                                               specificity_method = specificity_method,
                                               verbose            = verbose)

  list(
    node_metrics  = node_metrics,
    edge_data     = network$edge_data,
    pathway_genes = pathway_genes
  )
}