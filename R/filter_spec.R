#' filter_spec
#' A function that filters a spectrum object based on the number of peaks and
#'     intensity threshold.
#'
#' @importFrom Spectra filterIntensity applyProcessing mz intensity filterMzValues
#'     filterPrecursorPeaks
#' @param sps spectrum object. Input data
#' @param intensity_threshold numeric. minimum intensity threshold. Defaults to 0
#'     only meaningful for non-normalized data.
#' @param intensity_threshold_pct numeric. minimum intensity in % of max peak
#'     Peaks below will be removed.
#' @param min_peaks Integer. Minimum number of peaks after thresholding.
#' @param max_peaks Integer. Maximum number of peaks after thresholding.
#'     Choses the most intense peaks if more than max_peaks are present.
#' @param max_pct numeric. Maximum percentage of peaks to keep, if n_peaks > max_peaks.
#' @param remove_precursor logical. If TRUE, precursor peaks (+/- 10 mz) will be removed.
#' @param tolerance numeric. Tolerance for precursor peak removal in Da.
#' @param ppm numeric. ppm for precursor peak removal.
#' @returns A spectrum object with filtered peaks.
#' @export
filter_spec <- function(sps,
                        intensity_threshold = 0,
                        intensity_threshold_pct = 1,
                        min_peaks = 4,
                        max_peaks = 50,
                        max_pct = 0.95,
                        remove_precursor = TRUE,
                        tolerance = 0.005,
                        ppm = 20) {
  # copy input object
  tmp <- sps
  # generte list of filtered spectra
  output_list <- lapply(1:length(tmp), function(i) {
    # copy spectrum_i
    spectrum <- tmp[i]
    #
    if(remove_precursor) {
      # remove precursor peak
      spectrum <- filterPrecursorPeaks(spectrum,
                                       tolerance = tolerance,
                                       ppm = ppm,
                                       mz = ">=")
      #spectrum <- applyProcessing(spectrum)
    }
    # get intensities vector
    intensities <- unlist(intensity(spectrum))
    # filter for intensity
    min_intensity <- max(
      intensity_threshold_pct/100 * max(intensities),
      intensity_threshold
    )
    spectrum <- filterIntensity(
      spectrum,
      intensity = c(min_intensity,Inf)
    )
    # get number of peak in spectrum
    n_peaks <- length(unlist(mz(spectrum)))
    # remove if n_peaks < min_peaks
    if(n_peaks < min_peaks){
      return(NULL)
    }
    #
    if(n_peaks > max_peaks) {
      # get intensities vector
      intensities <- unlist(intensity(spectrum))
      # Get the indices of the most intense peaks
      idx <- order(intensities, decreasing = TRUE)[1:max_peaks]
      # get peaks to keep based on max_pct
      idx_new <- idx[
        1:which.max(
          cumsum(intensities[idx]/sum(intensities[idx])) > max_pct
        )
      ]
      # Subset the peaksData to keep only those rows
      spectrum <- filterMzValues(
        spectrum,
        mz = unlist(mz(spectrum))[idx_new],
        ppm = 0
      )
    }
    # return results
    return(applyProcessing(spectrum))
  })
  # remove empty entries
  output_list <- Filter(Negate(is.null), output_list)
  # combine list into single spectrum object
  output <- do.call(c, output_list)
  # return
  return(output)
}
