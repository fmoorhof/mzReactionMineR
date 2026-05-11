#' normalizeIS
#'
#' A function to normalize the intensities of an assay in a SummarizedExperiment
#'  object based on the intensities of an internal standard. The internal
#'  standard can be specified by its id, or by its rt and mz values. The
#'  function will find the closest feature to the provided rt and mz values, and
#'  use its intensities for normalization.
#'
#'  Normalization is defined as:
#'
#'    \deqn{\hat{Y_ij} = Y_ij\frac{\overline{S}{S_j}}}
#'
#'  Where \eqn{\hat{Y_ij}} is the normalized value of peak i in sample j,
#'  \eqn{Y_ij} is the original value  of peak i in sample j, \eqn{S_j} is the
#'  internal standard intensity in sample j and \eqn{\overline{S}} is the mean
#'  or median of the internal standard across samples.
#'
#' @importFrom dplyr %>% select mutate filter slice_min pull between across
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData assays<-
#' @importFrom stats median
#'
#' @param object A SummarizedExperiment object
#' @param assay Character. The name of the assay to be normalized.
#' @param rt Numeric. The theoretical retention time of the internal standard.
#' @param mz Numeric. The theoretical m/z of the internal standard.
#' @param id Character. The id of the internal standard. If "none", the function
#'        will find the closest feature to the provided rt and mz values.
#'        Default is "none".
#' @param type Either "mean" oder "median". How to calculate the reference IS
#'        value for normalization. Default is "median".
#' @param mz_tolerance Numeric vector of length 2. Absolute and relative (ppm)
#'        m/z tolerance for finding the internal standard.
#' @param rt_tolerance Numeric. Absolute retention time tolerance for finding
#'        the internal standard.
#' @param new_assay_name Character. The name of the new assay that will contain
#'        the normalized intensities. Default is "is_normalized".
#' @param remove Logical. Whether to remove the internal standard from the
#'        object after normalization. Default is TRUE.
#' @param id_col Character. The name of the column in rowData that contains the feature ids.
#' @param rt_col Character. The name of the column in rowData that contains the retention times.
#' @param mz_col Character. The name of the column in rowData that contains the m/z values.
#'
#' @returns A SummarizedExperiment object with a new assay containing the
#'        normalized intensities.
#' @export
normalizeIS <- function(
    object,
    assay = NULL,
    rt = NULL,
    mz = NULL,
    id = "none",
    type = "median",
    mz_tolerance = c(0.005, 10),
    rt_tolerance = 0.1,
    new_assay_name = "is_normalized",
    remove = TRUE,
    id_col = "id",
    rt_col = "rt",
    mz_col = "mz"
) {

  if(id != "none") {
    is_id <- id
  } else {

    is_id <- get_id(
      object = object,
      rt = rt,
      mz = mz,
      mz_range = calc_mz_range(mz, mz_tolerance),
      rt_range = calc_rt_range(rt, rt_tolerance),
      id_col = id_col,
      rt_col = rt_col,
      mz_col = mz_col
    )

  }

  is_intensities <- get_intensities_id(
    object = object,
    assay = assay,
    id = is_id,
    id_col = id_col
  )

  if(any(is.na(is_intensities))) {

    stop("Internal standard absent in some samples.")

  }

  correction_factor <- switch(
    type,
    mean = is_intensities/mean(is_intensities),
    median = is_intensities/median(is_intensities)
  )

  new_object <- divide_by_feature(
    object = object,
    assay = assay,
    vector = correction_factor,
    new_assay_name = new_assay_name
    )

  if(remove) {
    new_object <- new_object[-which(rowData(new_object)[[id_col]] == is_id), ]
  }

  return(new_object)

}

utils::globalVariables(c(".data", ".", "mz_diff", "rt_diff"))
