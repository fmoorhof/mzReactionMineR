#' metaboscape_to_se
#' A function that converts the output of a metaboscape analysis into a
#'        SummarizedExperiment Object.
#'
#' @importFrom utils read.csv
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @param path_to_file chr. Path to the csv file containing the metaboscape
#'        output.
#' @param sample_meta_data data.frame. Sample meta data used as colData.
#' @param assay_name chr. Name of the assay to be created from the data. Default is "intensity".
#' @param filename_column chr. The name of the column in the sample_meta_data
#'        that contains the column names of the metaboscape output. Default is "filename".
#'
#' @returns a SummarizedExperiment Object
#' @export
#'
metaboscape_to_se <- function(
    path_to_file,
    sample_meta_data,
    assay_name = "intensity",
    filename_column = "filename"
) {

  # parse sample_meta_data file from environment

  if(!filename_column %in% colnames(sample_meta_data)) {
    stop("filenames column not found in sample_meta_data.
         Please ensure that the column name specified in the 'filename_column'
         argument matches a column in your sample_meta_data.")
  }

  filenames <- sample_meta_data[[filename_column]]

  # ensure filenames are basenames without extension
  # bruker uses a *.d extension

  filenames <- gsub(pattern = "\\.d$", replacement = "", basename(filenames))

  # make names, because R might change filenames in the column of the
  # metaboscape output

  filenames <- make.names(filenames)

  sample_meta_data[[filename_column]] <- filenames

  # read in metaboscape output
  # data is stored in a simple csv file wit rowdata + filenames as column

  metaboscape_output <- read.csv(path_to_file)

  # row data is everything except the filename column

  row_data <- metaboscape_output[
    , !colnames(metaboscape_output) %in% filenames
  ]

  # assay data is contained in columns that match the filenames
  # calling it like that should also keep order

  if(!all(filenames %in% colnames(metaboscape_output))) {
    warning(
      paste0(
      "Not all filenames from sample_meta_data are present as columns in the
      metaboscape output.
      Please check for any missing entries in your sample_meta_data."
      ) # end paste
    ) # end warning

    filenames <- filenames[filenames %in% colnames(metaboscape_output)]
    sample_meta_data <- sample_meta_data[
      sample_meta_data[[filename_column]] %in% filenames,
    ]
    #

  } # end if

  assay_data <- metaboscape_output[,filenames]

  if(!all(colnames(assay_data) %in% filenames)) {
    warning(
      paste0(
        "Not all columns in the assay data match the filenames from  the
        sample_meta_data. Missing entries will be removed from the data. Please
        check for any unwanted discrepancies between the column names in the
        metaboscape output and the filenames in your sample_meta_data.\n
        Removing files:", print(!colnames(assay_data) %in% filenames)
      ) # end paste
    ) # end warning
  } # end if

  # generate a SummarizedExperiment object

  assay_list <- list(
    as.matrix(assay_data)
  )
  names(assay_list) <- assay_name

  output_se <- SummarizedExperiment(
    assays = assay_list,
    rowData = row_data,
    colData = sample_meta_data
  )

  return(output_se)

}
