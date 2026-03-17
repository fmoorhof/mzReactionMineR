# main.R - Example workflow for mzReactionMineR
# This script demonstrates the usage of the mzReactionMineR package

# Import of libraries required here in main file
library(mzReactionMineR)
library(SummarizedExperiment)
library(dplyr)
library(tidyr)

# ---- 1. Data Import ----
# Replace with the path to your mock or real data files
features_file <- "mock_data/peaklist52_mzmine.csv"  # peaklist5_mzmine.csv"
sample_meta_data_file <- "mock_data/plate5.csv"  # first column needs to be called filename fo QC plotting (hardcoded)

# Import data (mzmine_to_se), make sure to specify the correct field separator of .csv on error
sample_meta_data <- read.csv(sample_meta_data_file, sep = ",")
se <- mzmine_to_se(features_file,
                  sep = ",",
                  sample_meta_data = sample_meta_data,
)  # assays = "area")
cat("Imported SummarizedExperiment:\n")
print(se)


# ---- Blank filter ----
# Use blankSubtractionSE to filter features based on blank samples (last column in colData)
blank_sample <- tail(rownames(colData(se)), 1)
se_blanked <- blankSubtractionSE(
  object = se,
  assay = "area",
  sample_col = colnames(colData(se))[1],
  blanks = blank_sample,
  ratio = 3,
  ratio_type = "maximum",
  min_detection_blank = 1,
  id_col = "id"
)
cat("\nBlank filtering applied.\n")

# ---- 2. Data Processing ----
# Filter SummarizedExperiment (filterSe)
se_filtered <- filterSe(object = se_blanked,
                        assay = "area",
                        sample_col = colnames(colData(se_blanked))[1], # use first colData column as sample_col
                        # group_col = "none",
                        # not_in = tail(rownames(colData(se_blanked)), 1),
                        min_abundance = 10000,
                        min_pct = 0.0001,
                        min_n = 1L,
                        mz_range = c(200, 600),
                        rt_range = c(2.0, 6.0)
                        )
cat("\nFiltered SummarizedExperiment:\n")
print(se_filtered)

# ---- Export Filtered Data ----
# Export the filtered assay matrix (area) with rowData for mzMine re-import
filtered_export <- cbind(as.data.frame(rowData(se_filtered)), as.data.frame(assay(se_filtered, "area")))
write.csv(filtered_export, file = "filtered_for_mzmine.csv", row.names = FALSE)
cat("\nFiltered data exported to filtered_for_mzmine.csv\n")

# ---- 3. Normalization ----
# If normalizePQN returns a SummarizedExperiment, use it directly
norm_matrix <- normalizePQN(assay(se_filtered, "area"), measure = "median")
assay(se_filtered, "area") <- norm_matrix
se_norm <- se_filtered
cat("\nNormalized SummarizedExperiment:\n")
print(se_norm)

# ---- Export Normalized Data ----
# Export the normalized assay matrix (area) with rowData for mzMine re-import
norm_export <- cbind(as.data.frame(rowData(se_norm)), as.data.frame(assay(se_norm, "area")))
write.csv(norm_export, file = "normalized_for_mzmine.csv", row.names = FALSE)
cat("\nNormalized data exported to normalized_for_mzmine.csv\n")

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
  sep = ",",
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

# !!!!!maybe skip all the last part and invest in peak annotation approaches that i can integrate later
