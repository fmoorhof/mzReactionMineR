#' anovaLimma
#'
#' A wrapper function that performs an ANOVA test on a summarized experiment using limma.
#'
#' @importFrom dplyr %>% select
#' @importFrom limma lmFit eBayes topTable makeContrasts
#' @importFrom stats contr.sum as.formula model.matrix contrasts<-
#' @param object A summarized experiment object.
#' @param assay Character. What is the name of the assay the ANOVA should be
#'     executed on.
#' @param blocking_variables Character vector. The names of the columns in
#'     colData(object) that should be used as blocking factors
#' @param test_variables Character vector. The names of the columns in
#'     colData(object) that should be used as test variables.
#' @param padj_method Character. The method used to adjust the p-values.
#'     Default is "fdr".
#' @param return_colums Character vector. The columns from rowData(object) that are reported in results. Default is c("id", "rt", "mz").
#' @return Returns a data.frame that contains the ANOVA table.
#' @export

anovaLimma <- function(object = NULL,
                       assay = NULL,
                       blocking_variables = NULL,
                       test_variables = NULL,
                       padj_method = "fdr",
                       return_colums = c("id", "rt", "mz")) {

  sample_data <- as.data.frame(SummarizedExperiment::colData(object))

  formula_string <- paste(
    " ~ 0 +",
    paste(test_variables, collapse = " + ")
  )

  X <- model.matrix(
    data = sample_data,
    as.formula(formula_string)
  )

  for(name in test_variables){
    colnames(X) <- gsub(name, "", colnames(X))
  }

  if(!is.null(blocking_variables)) {

    Z <- list()

    for(block_name in blocking_variables) {

      block <- as.factor(sample_data[,block_name])

      contrasts(block) <- contr.sum(levels(block))

      Z[[block_name]] <- model.matrix(~block)[,-1,drop=FALSE]

      colnames(Z[[block_name]]) <- levels(block)[-1]

    }

    do.call(cbind, Z) -> Z

      X <- cbind(Z, X)

  }


  fit <- lmFit(
    assays(object)[[assay]] %>% replace(is.na(.),0),
    design = X
  )

  fit2 <- eBayes(fit)

  anova_res <- cbind(
    as.data.frame(SummarizedExperiment::rowData(object)),
    topTable(
      fit2,
      number=Inf,
      sort.by="none",
      adjust.method = "fdr",
      coef = colnames(X)[length(unique(sample_data[,unlist(blocking_variables)])):(length(colnames(X)))]
    )
  ) %>%
    select(return_colums, 'F', P.Value, adj.P.Val)

  return(anova_res)

}

utils::globalVariables(c(".data", ".", "id", "rt", "mz", "P.Value", "adj.P.Val"))
