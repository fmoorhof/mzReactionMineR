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


# ---- Blank Subtraction ----
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
cat("\nBlank subtraction applied.\n")

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

# ---- 3. Execute Python Script ----
# Run match_pools.py script which is processing filtered_for_mzmine.csv and creates matched_pool_areas.xlsx
system("python Python/match_pools.py")