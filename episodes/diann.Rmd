---
title: 'Working with DIA-NN'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- What is DIANN?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Learn what DIANN is

::::::::::::::::::::::::::::::::::::::::::::::::

## What is DIA-NN?

DIA-NN is designed to analyze **DIA / SWATH-MS** data — a type of mass spectrometry that captures fragment ions from all peptides within a sample, providing comprehensive proteome coverage.

Key features include:
- **Neural network–based scoring:** Improves identification confidence.
- **Interference correction:** Enhances quantification accuracy.
- **Library-free analysis:** Can work without a spectral library (generating one directly from DIA data).
- **High speed and efficiency:** Suitable for large-scale datasets.
- **Cross-run normalization:** Ensures consistent quantification across experiments.

For background information, see the original publication:  
*Demichev, V., et al. (2020). DIA-NN: neural networks and interference correction enable deep proteome coverage in high throughput. Nature Methods.*  
[https://doi.org/10.1038/s41592-019-0638-x](https://doi.org/10.1038/s41592-019-0638-x)

---

## General Workflow

The DIA-NN workflow typically consists of the following steps:

1. **Prepare input data**
   - Raw DIA data files (`.raw`, `.wiff`, or converted `.mzML` files).
   - (Optional) A spectral library (`.tsv` or `.speclib`).

2. **Run DIA-NN**
   - Can be run via the **graphical user interface (GUI)** or **command line**.
   - Specify inputs, output directory, and parameters (e.g., FDR cutoff, protein grouping, library mode).

3. **Output files**
   - Quantification tables for proteins and peptides (`report.tsv` or `.parquet`).
   - Run summaries and performance metrics.
   - Optional normalized intensity data for downstream statistical analysis.

4. **Post-processing**
   - Further analysis can be done in R (e.g., `tidyverse`, `MSstats`, `limma`), Python, or visualization tools.
   - Data can also be imported into platforms like **Perseus**, **Spectronaut**, or **Skyline** for downstream exploration.



## Example Command-Line Usage

Below is an example of running DIA-NN in command-line mode (adjust paths for your setup):

```bash
DIA-NN.exe --f "path/to/mzML_files/" \
           --out "path/to/output_folder/" \
           --lib "path/to/library.tsv" \
           --threads 8 \
           --fasta "path/to/fasta_file.fasta" \
           --qvalue 0.01 \
           --verbose 1
```

This command tells DIA-NN to:

Analyze all .mzML files in a folder

Use a specified spectral library and FASTA file

Run on 8 threads

Output results with a 1% false discovery rate cutoff

## Important Note

Due to the amount of time it will take to process the data, we will be using processed data in this workshop.

::::::::::::::::::::::::::::::::::::: keypoints 

- DIANN is fun and is theoretically easy to use!

::::::::::::::::::::::::::::::::::::::::::::::::

