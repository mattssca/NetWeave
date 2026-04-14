#' Build a dynamic gene blacklist from the network gene universe
#'
#' Identifies genes whose high STRING degree centrality reflects membership
#' in ubiquitous cellular machinery --- ribosomal complexes, histones,
#' splicing factors, and core translation factors --- rather than
#' subtype-specific signalling activity. Called internally by
#' \code{create_and_expand_network()} when \code{apply_blacklist = TRUE}.
#'
#' @param gene_universe Character vector of gene names present in the
#'   expanded network.
#' @param verbose Logical. If \code{TRUE}, prints the number and a preview
#'   of blacklisted genes identified.
#'
#' @return A character vector of gene symbols to exclude.
#' @keywords internal

get_blacklist <- function(gene_universe, verbose = TRUE) {
  
  # pattern-based families: genes whose high centrality reflects
  # translational machinery, chromatin structure, or RNA processing
  patterns <- c(
    "^RPL",   # Large ribosomal subunit
    "^RPS",   # Small ribosomal subunit
    "^RPLP",  # Ribosomal lateral stalk
    "^MRPL",  # Mitochondrial large subunit
    "^MRPS",  # Mitochondrial small subunit
    "^HIST",  # Canonical histones
    "^H1-",   # Histone H1 variants (new nomenclature)
    "^H2A",   # Histone H2A variants
    "^H2B",   # Histone H2B variants
    "^H3-",   # Histone H3 variants
    "^H4-",   # Histone H4 variants
    "^SNRP",  # Small nuclear ribonucleoproteins
    "^HNRNP", # Heterogeneous nuclear ribonucleoproteins
    "^PSM"   # Proteasome subunits
  )
  
  # explicit entries: genes not caught by patterns but known to cause
  # hub inflation in genome-wide PPI analyses
  explicit <- c(
    "UBC",    "UBB",    "UBA52",  "UBA80",   # Ubiquitin precursors
    "EEF1A1", "EEF1A2", "EEF2",              # Elongation factors
    "EIF4A1", "EIF4A2",                       # Initiation factors
    "GAPDH",  "PKM",    "ENO1",   "LDHA",    # Metabolic housekeeping
    "ACTB",   "ACTG1"                         # Ubiquitous cytoskeletal
  )
  
  # match patterns against the actual gene universe
  pattern_hits <- unlist(lapply(patterns, function(p) {
    gene_universe[grepl(p, gene_universe)]
  }))
  
  explicit_hits <- explicit[explicit %in% gene_universe]
  
  blacklist <- sort(unique(c(pattern_hits, explicit_hits)))
  
  if(verbose && length(blacklist) > 0){
    message(sprintf(
      "  -> Blacklist: %d gene(s) identified (%s%s)",
      length(blacklist),
      paste(head(blacklist, 5), collapse = ", "),
      ifelse(length(blacklist) > 5,
             sprintf(" ... and %d more", length(blacklist) - 5),
             "")
    ))
  }
  
  return(blacklist)
}