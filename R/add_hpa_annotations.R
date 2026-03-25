#' Add Human Protein Atlas (HPA) Annotations to Network Node Metrics
#'
#' Filters HPA normal tissue IHC and cancer expression data to the specified
#' tissues and left-joins the raw fields onto \code{node_metrics}.
#'
#' @param node_metrics A data frame of network node metrics. Must contain a
#'   \code{name} column with HGNC gene symbols.
#' @param this_tissue_normal A single character string specifying the normal
#'   tissue to filter on from the \code{Tissue} column of \code{normal_hpa}.
#' @param this_tissue_cancer A single character string specifying the cancer
#'   type to filter on from the \code{Cancer} column of \code{cancer_hpa}.
#' @param normal_hpa A data frame of HPA normal tissue IHC data. Expected
#'   columns include \code{Gene name}, \code{Gene}, \code{Tissue},
#'   \code{Cell type}, \code{Level}, \code{Reliability}, and
#'   \code{IHC tissue name}.
#' @param cancer_hpa A data frame of HPA cancer expression data. Expected
#'   columns include \code{Gene name}, \code{Gene}, \code{Cancer},
#'   \code{High}, \code{Medium}, \code{Low}, and \code{Not detected}.
#' @param verbose Logical. If \code{TRUE} (default), progress messages are
#'   printed.
#'
#' @return The \code{node_metrics} data frame with the following raw HPA
#'   columns appended:
#'   \describe{
#'     \item{hpa_normal_cell_type, hpa_normal_level,
#'       hpa_normal_reliability}{HPA normal tissue IHC fields.}
#'     \item{hpa_cancer_high, hpa_cancer_medium, hpa_cancer_low,
#'       hpa_cancer_not_detected}{HPA cancer patient count fields.}
#'   }
#'   Genes with multiple cell types in the normal data will produce multiple
#'   rows. Genes absent from HPA for the specified tissue will have \code{NA}
#'   for all added columns. Use \code{derive_hpa_metrics()} to compute derived
#'   annotations from these raw fields.
#'
#' @seealso \code{\link{derive_hpa_metrics}}
#'
#' @importFrom dplyr left_join
#'
#' @export
add_hpa_annotations <- function(node_metrics = NULL,
                                this_tissue_normal = NULL,
                                this_tissue_cancer = NULL,
                                normal_hpa = NULL,
                                cancer_hpa = NULL,
                                verbose = TRUE){
  
  # checks
  if(is.null(node_metrics)){
    stop("User must provide an incoming object with network node metrics...")
  }
  
  if(!"name" %in% colnames(node_metrics)){
    stop("node_metrics must contain a 'name' column...")
  }
  
  if(is.null(normal_hpa) || is.null(cancer_hpa)){
    stop("User must provide both normal_hpa and cancer_hpa data frames...")
  }
  
  if(!"Gene name" %in% colnames(normal_hpa)){
    stop(sprintf(
      "Expected column 'Gene name' not found in normal_hpa. Actual columns:\n  %s\n  Hint: re-read the TSV with check.names = FALSE",
      paste(colnames(normal_hpa), collapse = "\n  ")
    ))
  }
  
  if(!"Gene name" %in% colnames(cancer_hpa)){
    stop(sprintf(
      "Expected column 'Gene name' not found in cancer_hpa. Actual columns:\n  %s\n  Hint: re-read the TSV with check.names = FALSE",
      paste(colnames(cancer_hpa), collapse = "\n  ")
    ))
  }
  
  if(is.null(this_tissue_normal)){
    stop(sprintf(
      "User must provide this_tissue_normal. Available tissues:\n  %s",
      paste(get_hpa_tissues(normal_hpa), collapse = "\n  ")
    ))
  }
  
  if(is.null(this_tissue_cancer)){
    stop(sprintf(
      "User must provide this_tissue_cancer. Available cancer types:\n  %s",
      paste(get_hpa_tissues(cancer_hpa), collapse = "\n  ")
    ))
  }
  
  if(!this_tissue_normal %in% normal_hpa$Tissue){
    stop(sprintf(
      "Tissue '%s' not found in normal_hpa. Available tissues:\n  %s",
      this_tissue_normal, paste(get_hpa_tissues(normal_hpa), collapse = "\n  ")
    ))
  }
  
  if(!this_tissue_cancer %in% cancer_hpa$Cancer){
    stop(sprintf(
      "Cancer type '%s' not found in cancer_hpa. Available cancer types:\n  %s",
      this_tissue_cancer, paste(get_hpa_tissues(cancer_hpa), collapse = "\n  ")
    ))
  }
  
  if(verbose){
    message("Adding HPA annotations...")
    message(sprintf("  -> Normal tissue: %s", this_tissue_normal))
    message(sprintf("  -> Cancer type:   %s", this_tissue_cancer))
  }
  
  # filter and rename normal tissue data
  normal_hpa_tissue <- normal_hpa[normal_hpa$Tissue == this_tissue_normal, ]
  normal_hpa_tissue <- normal_hpa_tissue[, !colnames(normal_hpa_tissue) %in% c("Gene", "IHC tissue name", "Tissue")]
  
  names(normal_hpa_tissue)[names(normal_hpa_tissue) == "Gene name"]   <- "name"
  names(normal_hpa_tissue)[names(normal_hpa_tissue) == "Cell type"]   <- "hpa_normal_cell_type"
  names(normal_hpa_tissue)[names(normal_hpa_tissue) == "Level"]       <- "hpa_normal_level"
  names(normal_hpa_tissue)[names(normal_hpa_tissue) == "Reliability"] <- "hpa_normal_reliability"
  
  # filter and rename cancer data
  cancer_hpa_tissue <- cancer_hpa[cancer_hpa$Cancer == this_tissue_cancer, ]
  cancer_hpa_tissue <- cancer_hpa_tissue[, !colnames(cancer_hpa_tissue) %in% c("Gene", "Cancer")]
  
  names(cancer_hpa_tissue)[names(cancer_hpa_tissue) == "Gene name"]    <- "name"
  names(cancer_hpa_tissue)[names(cancer_hpa_tissue) == "High"]         <- "hpa_cancer_high"
  names(cancer_hpa_tissue)[names(cancer_hpa_tissue) == "Medium"]       <- "hpa_cancer_medium"
  names(cancer_hpa_tissue)[names(cancer_hpa_tissue) == "Low"]          <- "hpa_cancer_low"
  names(cancer_hpa_tissue)[names(cancer_hpa_tissue) == "Not detected"] <- "hpa_cancer_not_detected"
  
  if(verbose){
    message(sprintf("  -> Normal tissue rows: %d", nrow(normal_hpa_tissue)))
    message(sprintf("  -> Cancer rows:        %d", nrow(cancer_hpa_tissue)))
  }
  
  node_metrics %>%
    dplyr::left_join(normal_hpa_tissue, by = "name") %>%
    dplyr::left_join(cancer_hpa_tissue, by = "name")
  
}
