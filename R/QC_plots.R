#' QC_plots
#'
#' This function generates quality control plots from mzMine feature tables.
#' Inspired by xcms::plotQC.R
#'
#' @importFrom dplyr %>% mutate select filter group_by summarize_at
#' @importFrom utils read.csv
#' @importFrom SummarizedExperiment SummarizedExperiment assays rowData colData
#' @importFrom mzReactionMineR se_to_long
#' @importFrom methods setClass new
#' @importFrom grDevices nclass.Sturges
#' @import ggplot2
#' @param path_to_file character. Path to the mzMine feature table.
#' @param sample_meta_data data.frame. Sample metadata.
#'    The first column should be the sample names.
#' @param what Any of "rt_dev", "mz_dev", "rt_dev_sample", "mz_dev_sample".
#'     What plots to produce. Default is all.
#' @param return Logical. Whether to return the plotting data
#'
#' @returns a list object of class QCPlots that contains the summarized
#'     experiment and a list of ggplots
#' @export
#'
QC_plots <- function(
    path_to_file,
    sample_meta_data,
    sep = ",",
    what = c(
      "rt_dev",
      "mz_dev",
      "rt_dev_sample",
      "mz_dev_sample"
    ),
    return = FALSE
){

  # load feature data

  print("Loading mzMine feature table")

  features <- utils::read.csv(
    path_to_file,
    sep = sep
  ) %>%
    mutate(
      id = as.character(id) # convert id from int to character
    )

  # extract datafile names

  print("Processing feature table")

  sample_names <- names(features)[grepl("^datafile.*\\.mz$",names(features))]

  # rename to fit the file names in meta_data

  sample_names <- gsub(paste0("[.]mz"), "", sample_names)
  sample_names <- gsub(paste0("datafile[.]"), "", sample_names)

  sample_meta_data[,1] <- make.names(sample_meta_data[,1])

  samples <- intersect(
    sample_names,
    sample_meta_data[,1]
  )

  # create a list of necessary assays.

  assay_list <- list()

  # rt deviation

  if(any(c("rt_dev", "rt_dev_sample") %in% what)){
    assay_list[["rt_dev"]] <- features$rt - features[
      ,grepl("^datafile.*\\.rt$", colnames(features))
    ]

    names(assay_list[["rt_dev"]]) <-
      gsub(paste0("[.]rt"), "", names(assay_list[["rt_dev"]]))
    names(assay_list[["rt_dev"]]) <-
      gsub(paste0("datafile[.]"), "", names(assay_list[["rt_dev"]]))

    assay_list[["rt_dev"]] <- assay_list[["rt_dev"]][samples]

  }

  # mz deviation

  if(any(c("mz_dev", "mz_dev_sample") %in% what)){
    assay_list[["mz_dev"]] <- features$mz - features[
      ,grepl("^datafile.*\\.mz$", colnames(features))
    ]

    names(assay_list[["mz_dev"]]) <-
      gsub(paste0("[.]mz"), "", names(assay_list[["mz_dev"]]))
    names(assay_list[["mz_dev"]]) <-
      gsub(paste0("datafile[.]"), "", names(assay_list[["mz_dev"]]))

    assay_list[["mz_dev"]] <- assay_list[["mz_dev"]][samples]

  }


  # generate summarized experiment

  print("Creating SummarizedExperiment object")

  rowData_cols <- names(features)[!grepl("datafile[.]",names(features))]

  se <- SummarizedExperiment::SummarizedExperiment(
    rowData = features[,rowData_cols],
    assays = assay_list,
    colData = sample_meta_data[
      match(samples, sample_meta_data[,1]),
    ]
  )

  # generate list of specified ggplots

  plots <- list()

  if("mz_dev" %in% what){

    tmp <- mzReactionMineR::se_to_long(se, "mz_dev") %>%
      dplyr::filter(
        !is.na(.data$value)
      )

    breaks <- pretty(range(tmp$value),
                     n = nclass.Sturges(tmp$value),
                     min.n = 1)

    plots[["mz_dev"]] <- tmp %>%
      ggplot(
        aes(
          x = value
        )
      ) +
      geom_histogram(
        aes(y = after_stat(density)),
        alpha = 0.25,
        fill = "white",
        color = "black",
        breaks = breaks
      ) +
      labs(
        x = "m/z deviation",
        y = "Density",
        title = "m/z deviation"
      ) +
      theme_classic()

    print(plots[["mz_dev"]])

  }

  if("rt_dev" %in% what){
    tmp <- mzReactionMineR::se_to_long(se, "rt_dev") %>%
      dplyr::filter(
        !is.na(.data$value)
      )

    breaks <- pretty(range(tmp$value),
                     n = nclass.Sturges(tmp$value),
                     min.n = 1)

    plots[["rt_dev"]] <- tmp %>%
      ggplot(
        aes(
          x = value
        )
      ) +
      geom_histogram(
        aes(y = after_stat(density)),
        alpha = 0.25,
        fill = "white",
        color = "black",
        breaks = breaks
      ) +
      labs(
        x = "rt deviation",
        y = "Density",
        title = "rt deviation"
      ) +
      theme_classic()

    print(plots[["rt_dev"]])

  }


  if("rt_dev_sample" %in% what){

    plots[["rt_dev_sample"]] <- se_to_long(se, "rt_dev") %>%
      dplyr::group_by(
        filename
      ) %>%
      dplyr::summarize_at(
        "value",
        .funs = list(median),
        na.rm = TRUE
      ) %>%
      ggplot(
        aes(
          x = filename,
          y = value,
          fill = value
        )
      ) +
      geom_col() +
      labs(
        x = "Sample",
        y = "Median rt deviation",
        title = "rt deviation by sample"
      ) +
      scale_fill_viridis_c(option = "cividis") +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      )

    print(plots[["rt_dev_sample"]])

  }

  if("mz_dev_sample" %in% what){

    plots[["mz_dev_sample"]] <- se_to_long(se, "mz_dev") %>%
      group_by(
        filename
      ) %>%
      summarize_at(
        "value",
        .funs = list(median),
        na.rm = TRUE
      ) %>%
      ggplot(
        aes(
          x = filename,
          y = value,
          fill = value
        )
      ) +
      geom_col() +
      labs(
        x = "Sample",
        y = "Median m/z deviation",
        title = "m/z deviation by sample"
      ) +
      scale_fill_viridis_c(option = "cividis") +
      theme_classic() +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.x = element_blank()
      )

    print(plots[["mz_dev_sample"]])

  }

  output <- list(
    experiment = se,
    plots = plots
  )

  if(return == TRUE){
    return(output)
  }

}

utils::globalVariables(c(".data", ".", "..density..", "density", "filename",
                         "value"))
