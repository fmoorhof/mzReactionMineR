# main.R - Example workflow for mzReactionMineR
# This script demonstrates the usage of the mzReactionMineR package

# Load required packages
# Source all R scripts from the local R/ directory
r_scripts <- list.files("R", pattern = "\\.R$", full.names = TRUE)
sapply(r_scripts, source)
# when existing import paths should remain they need to be imported here in main such as library(dplyr), else dplyr:: import changes required?

# ---- 1. Data Import ----
# Replace with the path to your mock or real data files
features_file <- "mock_data/peaklist2_mzmine.csv"
sample_meta_data_file <- "mock_data/plate5.csv"

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
                        sample_col = "sampleNames",
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
