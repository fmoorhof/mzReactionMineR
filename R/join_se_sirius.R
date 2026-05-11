#' join_se_sirius
#'
#' A convenience function that joins a SIRIUS output file to a
#' SummarizedExperiment object.
#'
#' @importFrom SummarizedExperiment rowData SummarizedExperiment
#' @importFrom dplyr left_join
#' @importFrom utils read.delim
#'
#' @param object a SummarizedExperiment object
#' @param path_to_sirius path to the SIRIUS output file that should be joined
#' @param id_col character. The name of the column in the rowData that contains the feature ids. Default is "id".
#'
#' @returns a SummarizedExperiment object with the SIRIUS data joined to the
#' rowData
#' @export
#'
join_se_sirius <- function(
    object = NULL,
    path_to_sirius = NULL,
    id_col = "id"
) {
  canopus_structure_summary <- utils::read.delim(
    path_to_sirius
  ) %>%
    dplyr::rename(
      id_col = .data$mappingFeatureId
    ) %>%
    dplyr::mutate(
      id_col = as.character(.data[[id_col]])
    )
  new_rowData <- dplyr::left_join(
    base::as.data.frame(
      SummarizedExperiment::rowData(object)
    ),
    canopus_structure_summary
  )
  if(dim(SummarizedExperiment::rowData(object))[1] < dim(new_rowData)[1]){
    print(
      "Duplicate rows in the SIRIUS output file were found. Removing duplicates...",
    )
    SummarizedExperiment::rowData(object) <- new_rowData[!duplicated(new_rowData$id),]
  } else {
    SummarizedExperiment::rowData(object) <- new_rowData
  }
  return(object)
}
utils::globalVariables(c(".",".data"))
