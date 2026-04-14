#' List All Reactome Pathways Containing a Gene
#'
#' A discovery helper that takes a single HGNC gene symbol and prints all
#' Reactome pathways in which that gene is present, along with the number of
#' genes in each pathway. Useful for selecting an appropriate pathway seed
#' before calling \code{get_pathway_genes()}.
#'
#' @param gene_symbol A single character string. HGNC gene symbol to query.
#' @param verbose Logical. If \code{TRUE} (default), pathway names and gene
#'   counts are printed to the console.
#'
#' @return A sorted character vector of matching Reactome pathway names,
#'   returned invisibly.
#'
#' @details
#' Requires the \code{Ensembl2Reactome} data frame to exist in the global
#' environment (load with \code{data(Ensembl2Reactome)}). Gene symbol to
#' Ensembl ID mapping is performed via \code{org.Hs.eg.db} across all three
#' Ensembl ID types (\code{ENSG}, \code{ENSP}, \code{ENST}) to match the
#' mixed ID types present in the Reactome file. No internet connection required.
#'
#' @importFrom AnnotationDbi select
#' @export
list_gene_pathways <- function(gene_symbol, verbose = TRUE){
  
  if(!is.character(gene_symbol) || length(gene_symbol) != 1){
    stop("gene_symbol must be a single character string...")
  }
  
  if(!exists("Ensembl2Reactome")){
    stop("Ensembl2Reactome data is missing. Load it with data(Ensembl2Reactome)...")
  }
  
  # look up all three Ensembl ID types for the gene symbol
  # Reactome gene_id column contains mixed ENSG, ENSP, ENST entries
  all_ids <- character(0)
  
  ensg_map <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                          keys = gene_symbol, columns = "ENSEMBL", keytype = "SYMBOL"),
    error = function(e) NULL
  )
  if(!is.null(ensg_map)) all_ids <- c(all_ids, ensg_map$ENSEMBL)
  
  ensp_map <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                          keys = gene_symbol, columns = "ENSEMBLPROT", keytype = "SYMBOL"),
    error = function(e) NULL
  )
  if(!is.null(ensp_map)) all_ids <- c(all_ids, ensp_map$ENSEMBLPROT)
  
  enst_map <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db::org.Hs.eg.db,
                          keys = gene_symbol, columns = "ENSEMBLTRANS", keytype = "SYMBOL"),
    error = function(e) NULL
  )
  if(!is.null(enst_map)) all_ids <- c(all_ids, enst_map$ENSEMBLTRANS)
  
  all_ids <- unique(all_ids[!is.na(all_ids) & all_ids != ""])
  
  if(length(all_ids) == 0){
    stop(sprintf("Gene symbol '%s' not found in org.Hs.eg.db...", gene_symbol))
  }
  
  matched <- Ensembl2Reactome[Ensembl2Reactome$gene_id %in% all_ids, ]
  
  if(nrow(matched) == 0){
    message(sprintf("Gene '%s' was found but is not present in any Reactome pathways.", gene_symbol))
    return(invisible(character(0)))
  }
  
  pathway_names <- sort(unique(matched$pathway_name))
  
  if(verbose){
    message(sprintf("Gene '%s' is present in %d Reactome pathway(s):",
                    gene_symbol, length(pathway_names)))
    for(i in seq_along(pathway_names)){
      n_genes <- length(unique(grep("^ENSG",
                                    Ensembl2Reactome$gene_id[Ensembl2Reactome$pathway_name == pathway_names[i]],
                                    value = TRUE
      )))
      message(sprintf("  %d. %s (%d genes)", i, pathway_names[i], n_genes))
    }
  }
  
  invisible(pathway_names)
}