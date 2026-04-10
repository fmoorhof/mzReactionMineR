R package for downstream analysis of mass spectrometry data after peak picking with external software.

Funded by the German Federal Ministry of Education and Research and the state of Saxony-Anhalt as part of the project DiP-NA-WIR.

This package contains various functions used during analysis of 
    different metabolomics experiments. The main use is to simplify statistical
    workflows of output from other software like MZmine.

## Supported functionalities
```mermaid
graph TD
    A[Input] -->|MZmine data| B[Data Import]
    B --> C[Data Processing]
    C --> C1[Filtering]
    C --> C2[Normalization]
    C --> C3[Alignment]
    C --> C4[Statistical Analysis]
    C --> C5[Spectra Processing]
    C1 --> D1[filterSe]
    C1 --> D2[filterSe_ims]
    C1 --> D3[filter_spec]
    C2 --> D4[normalizePQN]
    C3 --> D5[join_aligner]
    C3 --> D6[join_aligner_ims]
    C3 --> D7[join_se_sirius]
    C4 --> D8[anovaLimma]
    C4 --> D9[contrastLimma]
    C5 --> D10[filterIntensity]
    C5 --> D11[filterMzValues]
    C5 --> D12[applyProcessing]
    C --> E[Output]
    E --> E1[QC_plots]
    E --> E2[se_to_long]
    E --> E3[mass_difference_network]
```

## Install
Requires R language
```
git clone https://github.com/DavidRuescher95/mzReactionMineR.git
cd mzReactionMineR

R
> install.packages("BiocManager")
> devtools::install(".")
```

## Usage
See `main.R` for an example minimal workflow of how to use and combine functionalities of this project.

## Troubleshooting install:
No R on system and no sudo access -> use R within conda:
```
conda create -n r -c conda-forge r-base r-languageserver r-devtools
conda activate r

R
> install.packages("BiocManager")
> devtools::install(".")
```
no dev tool package installed or location cwd (".") not found
```
install.packages("devtools")

devtools::install(
  "path/to/your/mzReactionMineR",
  dependencies = TRUE,
)
```

