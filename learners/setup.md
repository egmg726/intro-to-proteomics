---
title: Setup
---

## Workshop overview

This workshop will introduce participants to proteomic data processing, cleaning, and analysis, 
with a focus on **label-free, bottom-up, DIA, mass-spectrometry based quantitative proteomics**. 

We will begin with a brief overview of the proteomics workflow and introduction to the dataset used in this tutorial. 
This will be followed by an introduction to DIA-NN as an open-source DIA-focused proteomic data processing software. 
Participants will then be introduced to the *limpa* package for protein quantification and guided through key steps for 
cleaning and converting proteomics data to an analysis-ready format, including removing contaminants and quality filtering. 
Finally, we will introduce some basic statistical analyses for interpreting proteomics data, 
including differential expression analysis, enrichment analysis, and protein-protein interaction network analysis. 

<br>

:::: prereq
# LEARNING OBJECTIVES

By the end of this workshop, you will be able to: 

- Describe the proteomics workflow and understand the types of information that can be acquired from proteomic testing. 
- List the steps required to process and clean proteomics data, and justify which methods are most appropriate for your data. 
- Apply a variety of statistical analyses to a proteomics dataset and interpret their meaning. 

::::

<br>

### Timing

The anticipated workshop duration when delivered to a group of participants is **4 hours**. 

<br>

### Skill level

This introductory workshop assumes some knowledge of R and basic biology. No proteomics knowledge is assumed.

## Prior to the workshop

:::: checklist

## REQUIRED KNOWLEDGE

- This workshop assumes participants have a basic understanding of R or have previously attended 
an Intro to R workshop. Please review introductory materials [here](https://melbournebioinformatics.github.io/intro-to-r/).

::::

:::: checklist

## REQUIRED SOFTWARE

Attendees are required to bring their own laptop computers. Please ensure you have installed:

- [Chrome](https://www.google.com/chrome/) or [FireFox](https://www.mozilla.org/en-US/)
- [R](https://cran.ms.unimelb.edu.au/) (Download and install the latest version of R using the UniMelb mirror)
- [RStudio](https://posit.co/download/rstudio-desktop/#download)
- R packages required for this workshop (see below)

::::

### Installing required R packages

Please copy and run the below code to install the required R packages prior to the workshop.

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

<br>

:::: checklist
## REQUIRED DATA

In this workshop, we will be using [this dataset](https://www.ebi.ac.uk/pride/archive/projects/PXD047585) downloaded from PRIDE, 
a public repository for mass spectrometry-based proteomics data.

You can read the associated paper [here](https://dx.doi.org/10.3390/BIOMEDICINES12020333).

You do not need to download any data prior to the workshop.

::::


<!--
FIXME: place any data you want learners to use in `episodes/data` and then use
       a relative link ( [data zip file](data/lesson-data.zip) ) to provide a
       link to it, replacing the example.com link.
-->
