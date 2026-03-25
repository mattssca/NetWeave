#' Sjodahl 2017 Bladder Cancer Expression Data
#'
#' Microarray gene expression data from the Sjodahl 2017 bladder cancer cohort.
#' Rows are HGNC gene symbols, columns are sample identifiers.
#'
#' @format A numeric matrix with genes as rows (rownames = HGNC gene symbols)
#'   and samples as columns (colnames = sample IDs ending in \code{.CEL}).
#'
#' @source Sjodahl G et al. (2017). A molecular pathologic framework for
#'   risk stratification of stage T1 urothelial carcinoma.
#'   \emph{European Urology}, 71(6), 901–910.
#'   \doi{10.1016/j.eururo.2016.11.012}
#'
#' @seealso \code{\link{sjodahl_2017_subtypes}}
"sjodahl_2017"


#' Sjodahl 2017 Bladder Cancer Subtype Vector
#'
#' A named character vector mapping sample IDs to their molecular subtype
#' classification from the Sjodahl 2017 bladder cancer cohort.
#'
#' @format A named character vector where:
#'   \describe{
#'     \item{names}{Sample IDs matching \code{colnames(sjodahl_2017)}.}
#'     \item{values}{Molecular subtype labels (e.g. \code{"Uro"},
#'       \code{"GU"}, \code{"Mes"}, etc.).}
#'   }
#'
#' @seealso \code{\link{sjodahl_2017}}
"sjodahl_2017_subtypes"
