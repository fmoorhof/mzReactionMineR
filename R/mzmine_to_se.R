#' mzmine_to_se
#'
#' Convert mzMine feature table to a SummarizedExperiment object. Assumes the
#'     "Export to CSV (modular)" was used.
#'
#' @importFrom utils read.csv
#' @importFrom dplyr %>% mutate group_by
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
#' @param path_to_file path to the mzmine feature table
#' @param sep character. separator used in the mzmine .csv feature table
#' @param sample_meta_data data.frame. sample meta data to become colData
#'     the first columnn has to be the sample names
#' @param filenames character. name of the column in the meta data that
#'     corresponds to the sample names in the feature table
#' @param assays character vector. what result type to use. (i.e. "height")
#'
#' @returns a summarizedExperiment Object
#' @export
#'

mzmine_to_se <- function(
    path_to_file,
    sep = ",",
    sample_meta_data = NULL,
    assays = NULL,
    filenames = "filename"
) {

  # load feature data

  print("Loading mzMine feature table")

  features <- read.csv(
    path_to_file,
    sep = sep
  ) %>%
    mutate(
      id = as.character(id) # convert id from int to character 
    )

  # extract columns with specified assay name

  print("Processing feature table")


  if (is.null(assays)) {

    print("Assays not specified, using default: all")

    assays <- names(features[, grepl("datafile[.]", names(features))]) %>%
      sub(".*\\.", "", .) %>%
      unique()
    assays <- assays[!assays %in% c("min","max")]

  } else {

    tmp <- names(features[, grepl("datafile[.]", names(features))]) %>%
      sub(".*\\.", "", .) %>%
      unique()

    # if assays are manually specified, check if they are in the feature table
    if(!all(assays %in% tmp)) {
      stop("Assay names not found in feature table")
    } else {

      assays <- assays

    }
  }

  # generate list of data matrices

  assays_list <- lapply(assays, function(x) {
    tmp <- features[, grepl(paste0("[.]", x,"$"), names(features))]
    names(tmp) <- gsub(paste0("[.]",x), "", names(tmp))
    names(tmp) <- gsub(paste0("datafile[.]"), "", names(tmp))
    names(tmp) <- make.names(names(tmp))
    return(as.matrix(tmp))
  })
  names(assays_list) <- assays

  # test for correct sample names

  samples_equal <- all(sapply(assays_list[-1], function(df) {
    identical(colnames(df), colnames(assays_list[[1]]))
  }))

  if(!samples_equal) {

    stop("Sample names in the feature table are not equal across assays")

  } else {

    samples <- colnames(assays_list[[1]])

  }

  # rename to fit the file names in meta_data

  sample_meta_data[[filenames]] <- make.names(sample_meta_data[[filenames]])

  if(any(!samples %in% sample_meta_data[[filenames]])) {
    warning("Sample names in the feature table are not equal to sample names in the meta data")

    print(
      paste0(
        "Removing sample not present in meta data: ",
        samples[!samples %in% sample_meta_data[[filenames]]]
      )
    )

    sample_meta_data <- sample_meta_data[
      sample_meta_data[[filenames]] %in% samples,
    ]

    assays_list <- lapply(assays_list, function(x) {
      x[,sample_meta_data[[filenames]]]
    })

    samples <- colnames(assays_list[[1]])

  }

  sample_meta_data <- sample_meta_data[
    match(samples, sample_meta_data[[filenames]]),
  ]

  rownames(sample_meta_data) <- sample_meta_data[[filenames]]

  # generate summarized experiment

  print("Creating SummarizedExperiment object")

  rowData_cols <- names(features)[!grepl("datafile[.]",names(features))]

  se <- SummarizedExperiment::SummarizedExperiment(
    rowData = features[,rowData_cols],
    assays = assays_list,
    colData = sample_meta_data
  )

  return(se)

}


