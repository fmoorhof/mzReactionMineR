#' removeFeatureMz
#'
#' A function that removes features from a SummarizedExperiment object based on
#' m/z values. Useful to remove background signals (i.e. columns bleed).
#'
#' @importFrom SummarizedExperiment rowData SummarizedExperiment
#'
#' @param object a SummarizedExperiment object
#' @param mz a vector of m/z values to remove
#' @param mz_tolerance Numeric vector of length 2. Absolute and relative (ppm)
#'        m/z tolerance for finding the internal standard.
#' @param mz_col Character. The name of the column in rowData that contains the m/z values.
#' @param id_col Character. The name of the column in rowData that contains the ids.
#'
#' @returns a SummarizedExperiment object with the features removed
#' @export
#'
removeFeatureMz <- function(
    object,
    mz,
    mz_tolerance = c(0.005, 5),
    mz_col = "mz",
    id_col = "id"
) {

  # returns a vector of ids to remove
  id_to_remove <- sapply(seq_along(mz), function(i) {
    # function from utils.R
    get_id(
      object = object,
      rt = NULL,
      mz = mz[i],
      mz_range = calc_mz_range(mz[i], mz_tolerance),
      id_col = id_col,
      mz_col = mz_col
    )

  })

  id_to_remove <- unlist(id_to_remove)

  new_object <- object[!rowData(object)[[id_col]] %in% id_to_remove, ]

  return(new_object)

}
