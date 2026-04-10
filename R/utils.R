# A set of helper functions used in the main exported functions.


# function to find the closest id of a given rt/mz pair ------------------------

#' get_id
#'
#' Finds the id of the closest feature to the provided rt and mz values, within
#' the specified ranges. If rt is not provided, it will find the closest feature
#'  based on mz values only.
#'
#' @param object a SummarizedExperiment object
#' @param rt Numeric. Retention time value
#' @param mz Numeric. mz value
#' @param mz_range Numeric vector of length 2. The range of mz values to consider when finding the closest feature.
#' @param rt_range Numeric vector of length 2. The range of rt values to consider when finding the closest feature.umer
#' @param id_col Character. The name of the column in rowData that contains the feature ids.
#' @param rt_col Character. The name of the column in rowData that contains the retention times.
#' @param mz_col Character. The name of the column in rowData that contains the m/z values.
#'
#' @returns the id of the closest feature to the provided rt and mz values
#'
get_id <- function(
    object,
    rt = NULL,
    mz = NULL,
    mz_range = c(0, Inf),
    rt_range = c(0, Inf),
    id_col = "id",
    rt_col = "rt",
    mz_col = "mz"
) {

  if(is.null(rt)) {

    # find id
    result_id <- get_rowData(object) %>%
      filter(
        between(.data[[mz_col]], mz_range[1], mz_range[2])
      ) %>%
      mutate(
        mz_diff = abs(.data[[mz_col]]-mz)
      ) %>%
      slice_min(
        across(c(mz_diff)), n=1
      ) %>%
      pull(.data[[id_col]])

  } else {

    # find id
    result_id <- get_rowData(object) %>%
      filter(
        between(.data[[mz_col]], mz_range[1], mz_range[2]),
        between(.data[[rt_col]], rt_range[1], rt_range[2])
      ) %>%
      mutate(
        mz_diff = abs(.data[[mz_col]]-mz),
        rt_diff = abs(.data[[rt_col]]-rt)
      ) %>%
      slice_min(
        across(c(rt_diff,mz_diff)), n=1
      ) %>%
      pull(.data[[id_col]])

  }

  return(result_id)

}

# function to get an intensity vector from a specified assay and id ------------

get_intensities_id <- function(
    object,
    assay = NULL,
    id,
    id_col = "id"
) {

  values <- assays(object[which(rowData(object)[[id_col]] == id), ])[[assay]]
  return(as.numeric(values))

}

# functions to divide an assay by a vector of intensities ----------------------

divide_by_feature <- function(
    object,
    assay = NULL,
    vector = NULL,
    new_assay_name = NULL
) {

  new_object <- object
  assays(new_object, withDimnames = FALSE)[[new_assay_name]] <- t(apply(
    t(assays(new_object)[[assay]]),2,function(col) {
      col/vector
    }))

  return(new_object)

}

# function to calculate mz range based on absolute tolerance and ppm -----------

calc_mz_range <- function(
    mz,
    tolerance = c(0.005, 5)
) {

  mz_tolerance <- max(tolerance[1], (tolerance[2] * (mz/1e6)))
  mz_range <- c(mz - mz_tolerance, mz + mz_tolerance)
  return(mz_range)

}

# function to calculate rt range based on absolute tolerance and ppm -----------

calc_rt_range <- function(
    rt,
    tolerance = 0.1
) {

  rt_range <- c(rt - tolerance, rt + tolerance)
  return(rt_range)

}

# function to calculate cosine distance between rows of a matrix ---------------

calc_cosine_distance <- function(mat) {

  similarity <- (mat %*% t(mat)) / # calculate pairwise dotproduct
    # divide by the product of the length of the vectors
    outer(sqrt(rowSums(mat**2)),sqrt(rowSums(mat**2)))

  distance <- 1 - similarity

  return(distance)

}

# function to make a knn graph from a distance matrix --------------------------

#' make_knn_graph
#'
#' @importFrom igraph graph_from_adjacency_matrix
#' @param distance_mat distance matrix
#' @param k number of neighbors
#'
#' @returns an igraph object representing the knn graph
#'
make_knn_graph <- function(distance_mat, k) {

  # make 0 matrix to store knn graph
  knn_graph <- matrix(0, nrow = nrow(distance_mat), ncol = ncol(distance_mat))

  # loop over each row and find the k nearest neighbors of each point, excluding itself
  for(i in 1:nrow(knn_graph)) {

    knn_graph[i, order(distance_mat[i, ], decreasing = FALSE)[2:(k+1)]] <- 1

  }

  # make symmetric
  knn_graph <- pmax(knn_graph, t(knn_graph))

  # make igraph object from adjacency matrix
  knn_graph <- graph_from_adjacency_matrix(knn_graph, mode = "undirected" )

  return(knn_graph)

}
