#' Retrieve Genes Associated with a Reactome Pathway
#'
#' Queries the \code{Ensembl2Reactome} object (must exist in the global
#' environment) to retrieve all HGNC gene symbols associated with a given
#' Reactome pathway name or gene symbol. Matching is attempted in order:
#' exact pathway name, partial pathway name, then gene symbol lookup via
#' \code{biomaRt}.
#'
#' Alternatively, supply a pre-defined gene list via \code{gene_list} to
#' bypass Reactome entirely.
#'
#' @param this_pathway A single character string specifying either a Reactome
#'   pathway name (exact or partial match, case-insensitive) or an HGNC gene
#'   symbol. If a gene symbol is provided, all Reactome pathways containing
#'   that gene are returned. Ignored when \code{gene_list} is supplied.
#' @param gene_list An optional character vector of HGNC gene symbols. When
#'   provided, the Reactome lookup is skipped entirely and this vector is
#'   returned (sorted, unique, empty strings removed) as the pathway gene set.
#'   Takes precedence over \code{this_pathway}.
#' @param verbose Logical. If \code{TRUE}, progress and match messages are
#'   printed. Defaults to \code{FALSE}.
#' @param export_data Logical. If \code{TRUE}, the resulting gene vector is
#'   saved to an \code{.Rdata} file at \code{out_path}. Defaults to
#'   \code{FALSE}.
#' @param out_path Character string specifying the file path (without
#'   extension) for the exported \code{.Rdata} file. Required when
#'   \code{export_data = TRUE}; ignored otherwise.
#'
#' @return A sorted character vector of unique HGNC gene symbols associated
#'   with the matched pathway(s). Empty HGNC entries are excluded.
#'
#' @details
#' When \code{gene_list} is supplied it takes precedence and the function
#' returns immediately without touching Reactome or \code{biomaRt}.
#'
#' Otherwise, the function requires the \code{Ensembl2Reactome} data frame to
#' exist in the global environment with at least the columns
#' \code{pathway_name} and \code{gene_id} (Ensembl IDs). This object is
#' typically loaded via \code{data(Ensembl2Reactome)}.
#'
#' Matching proceeds as follows:
#' \enumerate{
#'   \item Exact case-insensitive match on \code{pathway_name}.
#'   \item Partial case-insensitive match on \code{pathway_name} via
#'     \code{grepl}.
#'   \item If no pathway match is found, \code{this_pathway} is treated as an
#'     HGNC gene symbol. A \code{biomaRt} query is used to convert it to an
#'     Ensembl ID, and all pathways containing that gene are returned.
#' }
#'
#' Ensembl IDs from the matched pathways are converted to HGNC symbols via a
#' \code{biomaRt} query against the \code{hsapiens_gene_ensembl} dataset.
#' An internet connection is required whenever \code{biomaRt} is invoked
#' (partial/gene-symbol search paths).
#'
#' @importFrom dplyr filter
#' @importFrom biomaRt useMart getBM
#'
#' @export
get_pathway_genes = function(this_pathway = NULL,
                             gene_list = NULL,
                             verbose = FALSE,
                             export_data = FALSE,
                             out_path = NULL){
  
  # check data
  if(is.null(this_pathway) && is.null(gene_list)){
    stop("User must provide a query pathway/gene name (this_pathway) or a pre-defined gene list (gene_list)...")
  }
  
  if(export_data && is.null(out_path)){
    stop("out_path must be provided when export_data = TRUE...")
  }
  
  # bypass Reactome if user supplies a gene list directly
  if(!is.null(gene_list)){
    pathway_hgnc <- sort(unique(gene_list[gene_list != ""]))
    if(verbose){
      message(sprintf("  -> Using user-supplied gene list (%d genes)", length(pathway_hgnc)))
    }
    if(export_data){
      save(pathway_hgnc, file = paste0(out_path, ".Rdata"))
    }
    return(pathway_hgnc)
  }
  
  if(!exists("Ensembl2Reactome")){
    stop("Ensembl to Reactome data is missing...")
  }
  
  # message
  if(verbose){
    message("Retrieving pathway genes...")
  }
  
  # determine if input is a pathway name or gene
  # first, try exact pathway name match (case-insensitive)
  pathway_all <- Ensembl2Reactome %>%
    filter(tolower(pathway_name) == tolower(this_pathway))
  
  # if no exact match, try partial match
  if(nrow(pathway_all) == 0){
    pathway_all <- Ensembl2Reactome %>%
      filter(grepl(this_pathway, pathway_name, ignore.case = TRUE))
    
    if(nrow(pathway_all) > 0 && verbose){
      message(sprintf("  -> No exact pathway match for '%s'", this_pathway))
      message(sprintf("  -> Found %d pathway(s) with partial match",
                      length(unique(pathway_all$pathway_name))))
    }
  } else {
    if(verbose){
      message(sprintf("  -> Found pathway match: '%s'",
                      unique(pathway_all$pathway_name)[1]))
    }
  }
  
  # if still no results, try searching by gene symbol
  if(nrow(pathway_all) == 0){
    if(verbose){
      message(sprintf("  -> No pathway found matching '%s'", this_pathway))
      message("  -> Searching as gene name...")
    }
    
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    
    gene_ensembl <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          filters = "hgnc_symbol",
                          values = this_pathway,
                          mart = mart)
    
    if(nrow(gene_ensembl) == 0){
      stop(sprintf("No pathway or gene found matching '%s'", this_pathway))
    }
    
    pathway_all <- Ensembl2Reactome %>%
      filter(gene_id %in% gene_ensembl$ensembl_gene_id)
    
    if(nrow(pathway_all) == 0){
      stop(sprintf("Gene '%s' found but not present in any Reactome pathways", this_pathway))
    }
    
    if(verbose){
      message(sprintf("  -> Found gene '%s' in %d pathways",
                      this_pathway, length(unique(pathway_all$pathway_name))))
    }
  } else {
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  }
  
  if(verbose && length(unique(pathway_all$pathway_name)) > 1){
    message("Pathways identified:")
    for(i in seq_along(unique(pathway_all$pathway_name))){
      message(sprintf("  %d. %s", i, unique(pathway_all$pathway_name)[i]))
    }
  }
  
  # convert Ensembl IDs to HGNC symbols
  pathway_ensembl <- unique(pathway_all$gene_id)
  
  pathway_hgnc <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                        filters = "ensembl_gene_id",
                        values = pathway_ensembl,
                        mart = mart)
  
  pathway_hgnc <- sort(unique(pathway_hgnc$hgnc_symbol[pathway_hgnc$hgnc_symbol != ""]))
  
  if(verbose){
    message(sprintf("  -> Retrieved %d unique genes from pathway(s)", length(pathway_hgnc)))
  }
  
  # export if requested
  if(export_data){
    save(pathway_hgnc, file = paste0(out_path, ".Rdata"))
  }
  
  return(pathway_hgnc)
  
}