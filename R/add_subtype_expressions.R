#' Add Subtype Mean Expression to Network Node Metrics
#'
#' Computes the mean expression for each subtype across network genes and appends
#' the results as new columns to a node metrics data frame. One column per subtype
#' is added, named \code{mean_expr_<subtype>}.
#'
#' @param this_subtype_vector A named character vector mapping sample IDs (names)
#'   to their subtype labels (values). Must be named.
#' @param subtype_order An optional character vector specifying the desired subtype
#'   order. Only subtypes present in \code{this_subtype_vector} are used. If
#'   \code{NULL}, subtypes are processed in the order returned by \code{unique()}.
#' @param expr_data A numeric matrix or data frame of expression values with genes
#'   as rows (rownames = gene IDs) and samples as columns (colnames = sample IDs).
#' @param node_metrics A data frame of network node metrics. Must contain a
#'   \code{name} column with gene identifiers matching the rownames of
#'   \code{expr_data}.
#' @param verbose Logical. If \code{TRUE} (default), progress messages are printed.
#'
#' @return The \code{node_metrics} data frame with additional columns
#'   \code{mean_expr_<subtype>} for each subtype, containing the mean expression
#'   value across samples of that subtype. Genes not found in \code{expr_data}
#'   will have \code{NA} for those columns.
#'
#' @details
#' Gene matching is performed by intersecting \code{node_metrics$name} with
#' \code{rownames(expr_data)}. Sample matching is performed by intersecting
#' the sample names for each subtype with \code{colnames(expr_data)}. If either
#' intersection is empty for a given subtype, a warning is issued and the column
#' is filled with \code{NA}.
#'
#' @export
#' 
add_subtype_expression = function(this_subtype_vector = NULL,
                                  subtype_order = NULL,
                                  expr_data = NULL,
                                  node_metrics = NULL,
                                  verbose = TRUE){
  
  #check data
  if(is.null(this_subtype_vector)){
    stop("User must provide a vector with sample/subtype information...")
  }
  
  if(is.null(names(this_subtype_vector))){
    stop("this_subtype_vector must be a named vector where names are sample IDs...")
  }
  
  if(is.null(expr_data)){
    stop("User must provide expression data...")
  }
  
  if(!is.matrix(expr_data) && !is.data.frame(expr_data)){
    stop("expr_data must be a matrix or data.frame...")
  }
  
  if(is.null(node_metrics)){
    stop("User must provide an incoming object with network node metrics...")
  }
  
  if(!"name" %in% colnames(node_metrics)){
    stop("node_metrics must contain a 'name' column...")
  }
  
  #message
  if(verbose){
    message("Adding subtype expression profile...")
  }
  
  #get unique subtypes present in data
  subtypes_present = unique(this_subtype_vector)
  
  #filter to only subtypes that exist in the data, maintaining desired order
  if(!is.null(subtype_order)){
    subtypes = subtype_order[subtype_order %in% subtypes_present]
  }else{
    subtypes = subtypes_present
  }
  
  if(verbose){
    message(sprintf("  -> Found %d subtypes: %s", 
                    length(subtypes), paste(subtypes, collapse = ", ")))
  }
  
  #loop through subtypes in specified order
  for(subtype in subtypes){
    
    col_name <- paste0("mean_expr_", subtype)
    
    if(col_name %in% colnames(node_metrics)){
      warning(sprintf("Column '%s' already exists in node_metrics and will be overwritten", col_name))
    }
    
    samples <- names(this_subtype_vector)[this_subtype_vector == subtype]
    valid_samples <- intersect(samples, colnames(expr_data))
    
    if(length(valid_samples) == 0){
      warning(sprintf("No valid samples found in expr_data for subtype '%s' — column will be NA", subtype))
      node_metrics[[col_name]] <- NA_real_
      next
    }
    
    genes_in_data <- intersect(node_metrics$name, rownames(expr_data))
    
    if(length(genes_in_data) == 0){
      warning("No network gene names matched rownames of expr_data — check that rownames are gene identifiers")
      node_metrics[[col_name]] <- NA_real_
      next
    }
    
    if(verbose){
      message(sprintf("  -> [%s] Using %d samples and %d genes",
                      subtype, length(valid_samples), length(genes_in_data)))
    }
    
    expr_sub <- expr_data[genes_in_data, valid_samples, drop = FALSE]
    mean_expr <- rowMeans(expr_sub, na.rm = TRUE)
    node_metrics[[col_name]] <- mean_expr[match(node_metrics$name, names(mean_expr))]
  }
  
  return(node_metrics)
  
}
