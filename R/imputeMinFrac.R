#' imputeMinFrac
#'
#' Function to impute missing values in a SummarizedExperiment assay with a
#' fraction of the minimum value for each row.
#'
#' @param object A SummarizedExperiment object
#' @param assay The assay to impute
#' @param fraction The fraction of the minimum value to impute (default is 5)
#' @param new_assay_name The name of the new assay to store the imputed values
#'        (default is "imputed")
#'
#' @returns a SummarizedExperiment object with a new assay containing the imputed values
#' @export
#'
imputeMinFrac <- function(
    object,
    assay,
    fraction = 5,
    new_assay_name = "imputed"
) {

  new_object <- object

  assays(new_object)[[new_assay_name]] <- t(apply(
    assays(new_object)[[assay]], 1, function(row) {
      replace(row, is.na(row), min(row, na.rm = TRUE)/5)
    }))

  return(new_object)

}
