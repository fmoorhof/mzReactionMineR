#' modified cosine.
#'
#' A wrapper function to calculate all pairwise modified cosine scores based on
#'     the Spectra package.
#'
#'@importFrom Spectra precursorMz peaksData acquisitionNum
#'@importFrom MsCoreUtils gnps join_gnps
#'
#' @param sps Spectra object. Must contain multiple MS/MS spectra
#' @param tolerance numeric. Absolute mz tolerance
#' @param ppm numeric. PPM tolerance for matching precursor m/z values
#'
#' @returns An upper triangular matrix of modified cosine scores
#' @export
#'
modified_cosine <- function(sps, tolerance = 0, ppm = 10) {
  # generate empty matrix
  modified_cosine <- matrix(
    NA_real_,
    nrow = length(sps),
    ncol = length(sps)
  )
  colnames(modified_cosine) <- acquisitionNum(sps)
  rownames(modified_cosine) <- acquisitionNum(sps)
  # loop over all spectra
  for(i in 1:length(sps)) {
    # define precursor mz and peak data for i
    spectrum_i <- peaksData(sps[i])[[1]]
    precursor_mz_i <- precursorMz(sps[i])
    # compare spectrum i against all other spectra
    for(j in 1:length(sps)) {
      # only generate scores for upper triangular matrix
      if(i < j) {
        # define precursor mz and peak data for j
        spectrum_j <- peaksData(sps[j])[[1]]
        precursor_mz_j <- precursorMz(sps[j])
        # use Spectra function to join MS/MS data
        map <- join_gnps(spectrum_i[, 1], spectrum_j[, 1], precursor_mz_i, precursor_mz_j, ppm = ppm)
        # calculate modified cosine
        modified_cosine[i,j] <- gnps(spectrum_i[map[[1]], ], spectrum_j[map[[2]], ])
      }
    }
  }
  return(modified_cosine)
}
