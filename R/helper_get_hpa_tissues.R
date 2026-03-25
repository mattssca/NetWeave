#' List Available Tissues in HPA Data
#'
#' Returns the unique tissue or cancer type labels present in an HPA data
#' frame. Useful for discovering valid values for the \code{this_tissue_normal}
#' and \code{this_tissue_cancer} arguments of \code{add_hpa_annotations()}.
#'
#' @param hpa_data A data frame of HPA data. Either a normal tissue IHC data
#'   frame (containing a \code{Tissue} column) or a cancer data frame
#'   (containing a \code{Cancer} column).
#'
#' @return A sorted character vector of unique tissue or cancer type labels.
#'
#' @seealso \code{\link{add_hpa_annotations}}
#'
#' @examples
#' \dontrun{
#' get_hpa_tissues(hpa_normal_data)
#' get_hpa_tissues(hpa_cancer_data)
#' }
#'
#' @export
get_hpa_tissues <- function(hpa_data = NULL){
  
  if(is.null(hpa_data)){
    stop("User must provide an HPA data frame...")
  }
  
  if("Tissue" %in% colnames(hpa_data)){
    return(sort(unique(hpa_data$Tissue)))
  }
  
  if("Cancer" %in% colnames(hpa_data)){
    return(sort(unique(hpa_data$Cancer)))
  }
  
  stop("hpa_data must contain either a 'Tissue' or 'Cancer' column...")
  
}