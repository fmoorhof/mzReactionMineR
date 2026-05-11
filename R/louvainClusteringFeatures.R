#' knnClusteringSamples
#'
#' A function that performs clustering (community detection) on features in a
#' SummarizedExperiment object based on correlation.
#'
#' @importFrom SummarizedExperiment rowData assays
#' @importFrom MatrixGenerics rowVars rowMaxs
#' @importFrom igraph cluster_louvain membership vertex_attr vertex_attr<- V
#'        degree delete_vertices subgraph
#' @importFrom stats prcomp cor
#' @importFrom dplyr right_join
#'
#' @param object a SummarizedExperiment object
#' @param assay Character. The name of the assay to be used for clustering.
#' @param R_trsh Numeric. Correlation threshold for an edge to be drawn.
#' @param resolution Numeric. Resolution parameter for the Louvain algorithm.
#'        Default is 1.
#' @param min_degree Integer. minimum number of edges for a feature to be kept.
#' @param min_cluster_size Integer. Minimum cluster size to be kept after clustering.
#' @param metric Either "pearson", "kendall", or "spearman". The correlation
#'        metric to use for calculating the correlation matrix. Default is "pearson".
#' @param type Either "pos", "neg", or "both". What kinds of correlation will be
#'        used for the graph: only positive, only negative or both.
#' @param filter_type Character. The method to select features for clustering.
#'        Must be either "intsnity" or "variance". Default is "area".
#' @param n_top Integer. The number of features to select based on
#'        filter_type.
#' @param calc_PCA Logical. Wether to perform PCA before clustering.
#'        Default is TRUE.
#' @param min_PC Integer. The minimum number of PCs to retain. Default is 5.
#' @param PC_var Numeric. The minimum cumulative variance that the retained PCs should explain.
#' @param scale Logical. Whether to scale the data before calculating the
#'        distance matrix. Default is TRUE.
#' @param id_col Character. The name of the column in rowData(object) that
#'        contains the feature ids. Default is "id".
#'
#' @returns a data.frame with two columns: "sample" and "cluster".
#' @export
#'
louvainClusteringFeatures <- function(
    object,
    assay,
    R_trsh = 0.8,
    resolution = 1,
    min_degree = 20,
    min_cluster_size = 50,
    metric = "pearson",
    type = "pos",
    filter_type = "intensity",
    n_top = nrow(object),
    return_removed = TRUE,
    calc_PCA = FALSE,
    min_PC = 5,
    PC_var = 0.8,
    scale = TRUE,
    id_col = "id"
) {

  # select ids based on rowwise variance

  print(paste0("Select top, ", n_top, " features based on ", filter_type))

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

    input_data <- (scale(t(input_data)))

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

  # generate graph -------------------------------------------------------------

  # calculate correlation matrix

  print("Calculating correlation matrix.")

  if(!metric %in% c("pearson", "kendall", "spearman")) {

    stop("Invalid value for 'metric'.
         Must be either 'pearson', 'kendall' or 'spearman'.")

  } else {

    cor_mat <- cor(input_data, method = metric)

  }

  # generate adjacency matrix

  adj_mat <- switch(
    type,
    pos = (cor_mat > R_trsh) * 1,
    neg = (cor_mat < R_trsh) * 1,
    both = (abs(cor_mat) > R_trsh) * 1
  )

  rm(cor_mat)

  print(paste0("Generating graph."))

  G <- graph_from_adjacency_matrix(adj_mat, mode = "undirected")
  vertex_attr(G)$name <- rowData(object)[[id_col]]

  rm(adj_mat)

  # graph clean-up -------------------------------------------------------------

  # low degree

  features_to_remove <- V(G)[degree(G) < min_degree]$name
  G_clean <- delete_vertices(G, features_to_remove)

  print(paste0(
    "Removed ", length(features_to_remove), " features with degree < ", min_degree)
  )

  # clustering -----------------------------------------------------------------

  print(paste0("Execute clustering using Louvain algorithm with resolution ", resolution, "."))

  communities <- cluster_louvain(G_clean, resolution = resolution)
  #
  clusters <- data.frame(
    V(G_clean)$name,
    as.character(membership(communities))
  )
  colnames(clusters) <- c(id_col, "cluster")

  # remove small clusters ------------------------------------------------------

  clusters_to_remove <- clusters %>%
    group_by(cluster) %>%
    summarize(n = n()) %>%
    filter(
      n < min_cluster_size
    )

  print(paste0(
    "Removed ", length(clusters_to_remove$cluster), " cluster with < ",
    min_cluster_size, " features, totalling ", sum(clusters_to_remove$n),
    " features."
  ))

  clusters_clean <- clusters %>%
    filter(!cluster %in% clusters_to_remove$cluster)

  # rename clusters

  clusters_clean <- clusters_clean %>%
    group_by(cluster) %>%
    summarize(n = n()) %>%
    ungroup() %>%
    arrange(desc(n)) %>%
    mutate(
      cluster_new = 1:nrow(.)
    ) %>%
    right_join(clusters_clean) %>%
    select(all_of(id_col), cluster_new) %>%
    dplyr::rename("cluster" = cluster_new)

  # calculate degree within community ------------------------------------------

  cluster_connectivity <- list()

  for(i in unique(clusters_clean$cluster)) {

    vertices_to_keep <- clusters_clean[clusters_clean$cluster == i,][[id_col]]

    connectivity <- degree(subgraph(G_clean, vertices_to_keep))

    cluster_connectivity[[i]] <- data.frame(
      names(connectivity),
      i,
      connectivity
    )
    names(cluster_connectivity[[i]]) <- c(id_col, "cluster", "cluster_connectivity")
  }

  clusters_clean <- bind_rows(cluster_connectivity)
  row.names(clusters_clean) <- NULL

  if(return_removed) {

    print("Adding removed features to output as cluster '-1'")

    removed_features <- data.frame(
      c(
        features_to_remove,
        clusters[clusters$cluster %in% clusters_to_remove$cluster,id_col]
      ),
      -1,
      NA
    )
    names(removed_features) <- c(id_col, "cluster", "cluster_connectivity")

    clusters_clean <- rbind(
      clusters_clean,
      removed_features
    )

  }

  return(clusters_clean)

}


utils::globalVariables(c("cluster", "desc", "cluster_new"))
