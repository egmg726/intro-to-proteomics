---
title: Setup
---


## Data Sets

<!--
FIXME: place any data you want learners to use in `episodes/data` and then use
       a relative link ( [data zip file](data/lesson-data.zip) ) to provide a
       link to it, replacing the example.com link.
-->

We will be using this dataset.

https://www.ebi.ac.uk/pride/archive/projects/PXD047585

Read the paper here: https://dx.doi.org/10.3390/BIOMEDICINES12020333





:::: prereq

Some knowledge of R is assumed, no proteomics knowledge is assumed.

You should review introductory materials here.

This lesson assumes you have R and RStudio installed on your computer.


::::


## Software Setup

::::::::::::::::::::::::::::::::::::::: discussion

If you not have R and RStudio already installed, please download them here:

[Download and install the latest version of R using the UniMelb mirror](https://cran.ms.unimelb.edu.au/).
[Download and install RStudio](https://posit.co/download/rstudio-desktop/#download).


:::::::::::::::::::::::::::::::::::::::::::::::::::


## Install Libraries

```r


# Packages from CRAN
cran_packages <- c(
  "limpa",           # Proteomics data processing and DE analysis
  "dplyr",           # Data manipulation
  "readxl",          # Read Excel files
  "curl",            # Download files from URLs
  "pheatmap",        # Heatmap visualization
  "EnhancedVolcano", # Volcano plots
  "STRINGdb",        # Protein-protein interaction network visualization
  "arrow"            # Dependency for .parquet reading in limpa
)

# Packages from Bioconductor
bioc_packages <- c(
  "clusterProfiler", # Functional enrichment analysis
  "org.Hs.eg.db",    # Human gene annotation (for GO/KEGG)
  "rpx",             # Interface to the ProteomeXchange Repository
  "STRINGdb"         #  Interface to the STRING protein-protein interactions database
)

# Function to install missing CRAN packages
install_if_missing <- function(pkgs, repo = "https://cloud.r-project.org") {
  to_install <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(to_install) > 0) {
    install.packages(to_install, repos = repo)
  }
}

# 1. Install CRAN packages
install_if_missing(cran_packages)

# 2. Install Bioconductor manager if needed
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# 3. Install Bioconductor packages
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = TRUE)
  }
}



```



