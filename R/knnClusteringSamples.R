#' knnClusteringSamples
#'
#' A function that performs clustering (community detection) on samples in a
#' SummarizedExperiment object based on a k-nearest neighbor graph.
#'
#' @importFrom SummarizedExperiment rowData assays
#' @importFrom MatrixGenerics rowVars rowMaxs
#' @importFrom igraph cluster_louvain membership
#' @importFrom stats prcomp dist
#'
#' @param object a SummarizedExperiment object
#' @param assay Character. The name of the assay to be used for clustering.
#' @param filter_type Character. The method to select features for clustering.
#'        Must be either "intsnity" or "variance". Default is "area".
#' @param n_top Integer. The number of features to select based on
#'        filter_type.
#' @param calc_PCA Logical. Wether to perform PCA before clustering.
#'        Default is TRUE.
#' @param min_PC Integer. The minimum number of PCs to retain. Default is 5.
#' @param PC_var Numeric. The minimum cumulative variance that the retained PCs should explain.
#' @param k Integer. The number of neighbors to use for the k-nearest neighbor graph. Default is 5.
#' @param resolution Numeric. Resolution parameter for the Louvain algorithm. Default is 1.
#' @param metric Either "euclidean", "manhattan" or "cosine". The distance
#'        metric to use for calculating the distance matrix. Default is "cosine".
#' @param scale Logical. Whether to scale the data before calculating the
#'        distance matrix. Default is TRUE.
#' @param id_col Character. The name of the column in rowData(object) that
#'        contains the feature ids. Default is "id".
#'
#' @returns a data.frame with two columns: "sample" and "cluster".
#' @export
#'
knnClusteringSamples <- function(
  object,
  assay,
  filter_type = "intensity",
  n_top = nrow(object),
  calc_PCA = TRUE,
  min_PC = 5,
  PC_var = 0.8,
  k = 5,
  resolution = 1,
  metric = "cosine",
  scale = TRUE,
  id_col = "id"
) {

  # select ids based on rowwise variance

  print(paste0("Select top, ", n_top, " features based on", filter_type))

  if(filter_type == "intensity") {
    top_ids <- get_rowData(object)[
      order(rowMaxs(assays(object)[[assay]]), decreasing = TRUE)[1:n_top],
    ][[id_col]]

  } else if(filter_type == "variance") {

    top_ids <- get_rowData(object)[
      order(rowVars(assays(object)[[assay]]), decreasing = TRUE)[1:n_top],
    ][[id_col]]

  } else {

    stop("Invalid value for 'filter_type'. Must be either 'intensity' or 'variance'.")

  }


  input_data <- assays(object[rowData(object)[[id_col]] %in% top_ids, ])[[assay]]

  if(scale) {

    print("Scaling data.")

    input_data <- scale(t(input_data))
  } else {

    input_data <- t(input_data)

  }

  if(calc_PCA) {

    print("Calculating PCA.")

    pca_res <- prcomp(input_data)

    var_explained <- cumsum((pca_res$sdev^2) / sum(pca_res$sdev^2))

    num_PCs <- max(
      which.max(var_explained > PC_var),
      min_PC
    )

    print(paste0(
      "Retaining ", num_PCs, " PCs, explaining ",
      round(var_explained[num_PCs] * 100, 2), "% of the total variance."
    ))

    input_data <- pca_res$x[, 1:num_PCs]

  }

  # calculate distance matrices

  print("Calculating distance matrix.")

  if(metric %in% c("euclidean", "manhattan")) {

    distance_mat <- as.matrix(dist(input_data, method = metric))

  } else if(metric == "cosine") {

    # cosine similarity
    distance_mat <- calc_cosine_distance(input_data)

  } else {

    stop("Invalid value for 'metric'. Must be either 'euclidean', 'manhattan' or 'cosine'.")
}
  # generate knn graph

  print(paste0("Generating knn (k=", k, ") graph."))

  knn_graph <- make_knn_graph(distance_mat, k)

  # clustering

  print(paste0("Execute clustering using Louvain algorithm with resolution ", resolution, "."))

  knn_communities <- cluster_louvain(knn_graph, resolution = resolution)
  #
  knn_clusters <- data.frame(
    sample = row.names(distance_mat),
    cluster = as.character(membership(knn_communities))
  )

  return(knn_clusters)

}
