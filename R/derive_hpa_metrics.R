#' Derive HPA Expression Metrics from Raw Annotation Columns
#'
#' Computes derived expression classifications, confidence scores, filter flags,
#' and prioritisation flags from the raw HPA columns added by
#' \code{add_hpa_annotations()}.
#'
#' @param node_metrics A data frame containing the raw HPA columns produced by
#'   \code{add_hpa_annotations()}: \code{hpa_normal_level},
#'   \code{hpa_normal_reliability}, \code{hpa_cancer_high},
#'   \code{hpa_cancer_medium}, \code{hpa_cancer_low}, and
#'   \code{hpa_cancer_not_detected}.
#' @param verbose Logical. If \code{TRUE} (default), a summary of annotated
#'   genes and filter pass rate is printed.
#'
#' @return The input data frame with the following columns appended:
#'   \describe{
#'     \item{hpa_normal_expressed}{Logical; \code{TRUE} if normal level is
#'       High or Medium.}
#'     \item{hpa_normal_reliable}{Reliability tier: \code{"High"},
#'       \code{"Medium"}, \code{"Low"}, or \code{NA}.}
#'     \item{hpa_cancer_expressed}{Logical; \code{TRUE} if any High or Medium
#'       cancer staining is recorded.}
#'     \item{hpa_cancer_confidence}{\code{"High"}, \code{"Medium"},
#'       \code{"Low"}, or \code{"None"} based on combined High+Medium count.}
#'     \item{hpa_expression_pattern}{One of \code{"Constitutive"},
#'       \code{"Lost_in_cancer"}, \code{"Cancer_upregulated"},
#'       \code{"Not_expressed"}, \code{"Insufficient_data"}, or
#'       \code{"Unknown"}.}
#'     \item{hpa_data_quality}{\code{"High"}, \code{"Medium"}, \code{"Low"},
#'       or \code{"None"}.}
#'     \item{hpa_filter_has_data, hpa_passes_filter}{Logical filter flags.}
#'     \item{hpa_high_confidence, hpa_therapeutic_target,
#'       hpa_potential_suppressor}{Logical prioritisation flags.}
#'   }
#'
#' @seealso \code{\link{add_hpa_annotations}}
#'
#' @importFrom dplyr mutate case_when
#'
#' @export
derive_hpa_metrics <- function(node_metrics = NULL,
                               verbose = TRUE){
  
  # checks
  if(is.null(node_metrics)){
    stop("User must provide an incoming object with network node metrics...")
  }
  
  required_cols <- c("hpa_normal_level", "hpa_normal_reliability",
                     "hpa_cancer_high", "hpa_cancer_medium",
                     "hpa_cancer_low", "hpa_cancer_not_detected")
  missing_cols <- required_cols[!required_cols %in% colnames(node_metrics)]
  
  if(length(missing_cols) > 0){
    stop(sprintf(
      "Missing required HPA column(s): %s\n  Run add_hpa_annotations() first.",
      paste(missing_cols, collapse = ", ")
    ))
  }
  
  if(verbose) message("Deriving HPA expression metrics...")
  
  result <- node_metrics %>%
    mutate(
      
      # classify normal expression
      hpa_normal_expressed = case_when(
        is.na(hpa_normal_level)                        ~ NA,
        hpa_normal_level %in% c("High", "Medium")      ~ TRUE,
        hpa_normal_level %in% c("Low", "Not detected") ~ FALSE,
        TRUE                                           ~ NA
      ),
      
      # classify normal reliability
      hpa_normal_reliable = case_when(
        is.na(hpa_normal_reliability)                           ~ NA_character_,
        hpa_normal_reliability %in% c("Enhanced", "Supported") ~ "High",
        hpa_normal_reliability == "Approved"                   ~ "Medium",
        hpa_normal_reliability == "Uncertain"                  ~ "Low",
        TRUE                                                   ~ NA_character_
      ),
      
      # classify cancer expression
      hpa_cancer_expressed = case_when(
        is.na(hpa_cancer_high)                              ~ NA,
        (hpa_cancer_high + hpa_cancer_medium) > 0          ~ TRUE,
        hpa_cancer_low > 0 | hpa_cancer_not_detected >= 10 ~ FALSE,
        TRUE                                               ~ NA
      ),
      
      # cancer expression confidence
      hpa_cancer_confidence = case_when(
        is.na(hpa_cancer_high)                     ~ "None",
        (hpa_cancer_high + hpa_cancer_medium) >= 6 ~ "High",
        (hpa_cancer_high + hpa_cancer_medium) >= 3 ~ "Medium",
        (hpa_cancer_high + hpa_cancer_medium) >= 1 ~ "Low",
        TRUE                                       ~ "None"
      ),
      
      # overall expression pattern
      hpa_expression_pattern = case_when(
        is.na(hpa_normal_expressed) | is.na(hpa_cancer_expressed) ~ "Insufficient_data",
        hpa_normal_expressed  & hpa_cancer_expressed              ~ "Constitutive",
        hpa_normal_expressed  & !hpa_cancer_expressed             ~ "Lost_in_cancer",
        !hpa_normal_expressed & hpa_cancer_expressed              ~ "Cancer_upregulated",
        !hpa_normal_expressed & !hpa_cancer_expressed             ~ "Not_expressed",
        TRUE                                                      ~ "Unknown"
      ),
      
      # overall data quality
      hpa_data_quality = case_when(
        !is.na(hpa_normal_level) & !is.na(hpa_cancer_high) &
          hpa_normal_reliable %in% c("High", "Medium") ~ "High",
        !is.na(hpa_normal_level) & !is.na(hpa_cancer_high) ~ "Medium",
        is.na(hpa_normal_level)  | is.na(hpa_cancer_high)  ~ "Low",
        TRUE                                               ~ "None"
      ),
      
      # filter flags
      hpa_filter_has_data  = hpa_data_quality %in% c("High", "Medium"),
      hpa_passes_filter    = hpa_filter_has_data & (
        (hpa_normal_expressed %in% TRUE & hpa_normal_reliable %in% c("High", "Medium")) |
          hpa_cancer_confidence %in% c("High", "Medium")
      ),
      
      # prioritisation flags
      hpa_high_confidence = (
        hpa_data_quality == "High" &
          hpa_cancer_confidence == "High" &
          hpa_normal_reliable == "High"
      ),
      hpa_therapeutic_target = (
        hpa_expression_pattern == "Cancer_upregulated" &
          hpa_cancer_confidence %in% c("High", "Medium")
      ),
      hpa_potential_suppressor = (
        hpa_expression_pattern == "Lost_in_cancer" &
          hpa_normal_expressed %in% TRUE
      )
    )
  
  if(verbose){
    n_normal_expressed  <- sum(result$hpa_normal_expressed %in% TRUE, na.rm = TRUE)
    n_cancer_expressed  <- sum(result$hpa_cancer_expressed %in% TRUE, na.rm = TRUE)
    n_both_expressed    <- sum(result$hpa_normal_expressed %in% TRUE & 
                                 result$hpa_cancer_expressed %in% TRUE, na.rm = TRUE)
    n_pass              <- sum(result$hpa_passes_filter, na.rm = TRUE)
    
    message(sprintf("  -> Genes expressed in normal tissue:  %d", n_normal_expressed))
    message(sprintf("  -> Genes expressed in tumour tissue:  %d", n_cancer_expressed))
    message(sprintf("  -> Genes expressed in both:           %d", n_both_expressed))
  }
  
  return(result)
  
}