#' blankSubtractionSE
#'
#' Use blank samples to filter features in a SummarizedExperiment object based on a specified ratio of sample signal to blank signal.
#'
#' @importFrom dplyr %>% mutate group_by
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
#'
#' @param object A SummarizedExperiment object containing the data to be filtered.
#' @param assay  Character. what assay to use. (i.e. "height")
#' @param sample_col Character. Column name in colData that contains sample names (i.e. "filename")
#' @param blanks Character vector. Names of blank samples in the sample_col column of colData
#' @param ratio  Numeric. Minimum ratio of sample signal to blank signal for a feature to be retained. Default is 3.
#' @param ratio_type Character. Method to calculate the blank signal for each feature. Options are "maximum", "mean", or "median". Default is "maximum".
#' @param min_detection_blank Integer. Minimum number of blank samples in which a feature must be detected (i.e. non-NA value) to be considered for filtering. Default is 1.
#' @param remove Logical. Whether to remove the blank samples from the resulting SummarizedExperiment object. Default is TRUE.
#' @param id_col Character. Column name in rowData that contains unique feature identifiers (i.e. "id")
#' @param sample_col Character. Column name in rowData that contains unique feature identifiers (i.e. "sample")
#'
#' @returns A SummarizedExperiment object with features filtered based on the specified ratio.
#' @export
#'

blankSubtractionSE <- function(
    object = NULL,
    assay = NULL,
    blanks = NULL,
    ratio_type = "maximum",
    ratio = 3,
    min_detection_blank = 1,
    remove = TRUE,
    id_col = "id",
    sample_col = "filename"
) {


  # get column names of rowData

  id_cols <- names(as.data.frame(rowData(object)))

  # make long data frame of blank samples and calculate the blank value for each feature

  data_blank <- cbind(
    as.data.frame(rowData(object)),
    as.data.frame(assays(object)[[assay]])
  ) %>%
    pivot_longer(
      cols = -all_of(id_cols),
      names_to = sample_col,
      values_to = "Value"
    ) %>%
    inner_join(
      as.data.frame(colData(object))
    ) %>%
    filter(
      .data[[sample_col]] %in% blanks,
      !is.na(.data[["Value"]])
    ) %>%
    group_by(.data[[id_col]]) %>%
    filter(n() >= min_detection_blank)


  # calculate ratio type

  if(ratio_type == "maximum"){
    data_blank <- data_blank %>%
      summarize(Value = max(.data$Value, na.rm = TRUE))
  } else if(ratio_type == "mean"){
    data_blank <- data_blank %>%
      summarize(Value = mean(.data$Value, na.rm = TRUE))
  } else if(ratio_type == "median"){
    data_blank <- data_blank %>%
      summarize(Value = median(.data$Value, na.rm = TRUE))
  } else {
    stop("ratio_type must be one of 'maximum', 'mean', or 'median'")
  }
    # make long data frame of all samples and filter based on blank values


    data_samples <- cbind(
      as.data.frame(rowData(object)),
      as.data.frame(assays(object)[[assay]])
    ) %>%
      pivot_longer(
        cols = -all_of(id_cols),
        names_to = sample_col,
        values_to = "Value"
      ) %>%
      inner_join(
        as.data.frame(colData(object))
      ) %>%
      filter(
        !(.data[[sample_col]] %in% blanks),
        !is.na(.data[["Value"]]),
        .data[[id_col]] %in% data_blank[[id_col]]
      ) %>%
      group_by(.data[[id_col]])


    # calculate ratio type

    if(ratio_type == "maximum"){
      data_samples <- data_samples %>%
        summarize(Value = max(.data$Value, na.rm = TRUE))
    } else if(ratio_type == "mean"){
      data_samples <- data_samples %>%
        summarize(Value = mean(.data$Value, na.rm = TRUE))
    } else if(ratio_type == "median"){
      data_samples <- data_samples %>%
        summarize(Value = median(.data$Value, na.rm = TRUE))
    } else {
      stop("ratio_type must be one of 'maximum', 'mean', or 'median'")
    }

    # filter based on ratio

    ids_to_remove <- data_blank %>%
      left_join(
        data_samples,
        by = id_col,
        suffix = c("_blank", "_sample")
      ) %>%
      filter(.data$Value_sample < (ratio * .data$Value_blank) | is.na(.data$Value_sample)) %>%
      pull(id_col)

    if(length(ids_to_remove) > 0){
      result <- object[!rowData(object)[[id_col]] %in% ids_to_remove,]
    } else {
      warning("No features were filtered based on the provided criteria. Returning original object.")
      result <- object
    }

    if(remove) {
      result <- result[,!colData(result)[[sample_col]] %in% blanks]
    }

    return(result)


  } # end

 utils::globalVariables(c(".",".data","Value"))
