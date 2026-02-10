# main.R - Example workflow for mzReactionMineR
# This script demonstrates the usage of the mzReactionMineR package

# Import of libraries required here in main file
library(mzReactionMineR)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)

# ---- 1. Data Import ----
# Replace with the path to your mock or real data files
features_file <- "mock_data/peaklist2_mzmine.csv"
sample_meta_data_file <- "mock_data/plate5.csv"  # first column needs to be called filename fo QC plotting (hardcoded)

# Import data (mzmine_to_se), make sure to specify the correct field separator of .csv on error
sample_meta_data <- read.csv(sample_meta_data_file, sep = ";")
se <- mzmine_to_se(features_file,
                  sep = ";",
                  sample_meta_data = sample_meta_data,
)  # assays = "area")
cat("Imported SummarizedExperiment:\n")
print(se)

# ---- 2. Data Processing ----
# Filter SummarizedExperiment (filterSe)
se_filtered <- filterSe(object = se,
                        assay = "area",
                        sample_col = colnames(colData(se))[1], # use first colData column as sample_col
                        group_col = "none",
                        not_in = "none",
                        min_abundance = 10,
                        min_pct = 0.8,
                        min_n = 1L,
                        mz_range = c(100, 250),
                        rt_range = c(1.0, 2.5)
                        )
cat("\nFiltered SummarizedExperiment:\n")
print(se_filtered)

# ---- 3. Normalization ----
# If normalizePQN returns a SummarizedExperiment, use it directly
norm_matrix <- normalizePQN(assay(se_filtered, "area"), measure = "median")
assay(se_filtered, "area") <- norm_matrix
se_norm <- se_filtered
cat("\nNormalized SummarizedExperiment:\n")
print(se_norm)

# ---- 4. Statistical Analysis ----
# Example: ANOVA using Limma (anovaLimma)
if ("Group" %in% colnames(colData(se_norm))) {
  anova_res <- anovaLimma(
    object = se_norm,
    assay = "area",
    blocking_variables = NULL,
    test_variables = c("Group")
  )
  cat("\nANOVA Results:\n")
  print(head(anova_res))
} else {
  cat("\nNo 'Group' column found in colData. Skipping ANOVA.\n")
}

# ---- 5. Visualization ----
# QC Plots (QC_plots)
qc_result <- QC_plots(
  path_to_file = features_file, # path to the mzMine feature table
  sample_meta_data = sample_meta_data,
  sep = ";",
  what = c("rt_dev", "rt_dev_sample"),  # mz fails on: Error in FUN(left, right) : non-numeric argument to binary operator
  return = TRUE
)
# Save each ggplot in the returned list to a file
if (!dir.exists("./qc_plots_R/")) dir.create("./qc_plots_R/")
for (plot_name in names(qc_result$plots)) {
  ggplot2::ggsave(
    filename = paste0("./qc_plots_R/", plot_name, ".pdf"),
    plot = qc_result$plots[[plot_name]],
    width = 7, height = 5
  )
}
cat("\nQC plots saved to ./qc_plots_R/\n")
