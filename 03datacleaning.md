---
title: 'Data cleaning with R'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- How do we convert our processed proteomics data (DIA-NN output) into an analysis-ready format?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Understand the steps required to clean proteomics data and generate analysis-ready protein values using `limpa`.

::::::::::::::::::::::::::::::::::::::::::::::::

## Load processed data and sample annotation 

After processing raw proteomics data, further cleaning and analysis can be conducted using **R**, **Python**, or by importing data to other platforms like **Perseus** or **Skyline**.

In this tutorial, we will be teaching you one method of cleaning and analysing proteomics data using R. **To follow along the remainder of this tutorial, you should open a new R script in RStudio and run the example code provided.**

:::::callout
# The limpa package

In this tutorial, we are using the `limpa` package to load and clean our data. This package was published in 2025 by the [Smyth Lab](https://github.com/SmythLab/limpa), based at the Walter and Eliza Hall Institute of Medical Research (WEHI). 

For background information, see the original publications:

    - Li, M., Cobbold, S. A., & Smyth, G. K. (2025). Quantification and differential analysis of mass spectrometry proteomics data with probabilistic recovery of information from missing values. bioRxiv, 2025.2004.2028.651125. https://doi.org/10.1101/2025.04.28.651125 

    - Li, M., & Smyth, G. K. (2023). Neither random nor censored: estimating intensity-dependent probabilities for missing values in label-free proteomics. Bioinformatics, 39(5). https://doi.org/10.1093/bioinformatics/btad200 

You can also find further information in the [vignette](https://bioconductor.org/packages/release/bioc/vignettes/limpa/inst/doc/limpa.html) or [reference manual](https://bioconductor.org/packages/release/bioc/manuals/limpa/man/limpa.pdf) accessed via the [Bioconductor page](https://bioconductor.org/packages/release/bioc/html/limpa.html).

:::::

### Setup

First, load the required R packages.


``` r
library(limpa)
library(readxl)
library(pheatmap)
library(EnhancedVolcano)
library(STRINGdb)
library(clusterProfiler)
library(org.Hs.eg.db)
library(arrow)
library(rpx)
library(dplyr)
library(stringr)
```

### Load DIA-NN output

**Load the peptide-level data** processed by DIA-NN (`.parquet` file), **filtering observations** by quality metrics as recommended in the [limpa vignette](https://bioconductor.org/packages/release/bioc/vignettes/limpa/inst/doc/limpa.html). The peptide values will automatically be **converted to log2** upon import.


``` r
y.peptide <- readDIANN(file='data/MBIntroToProteomics.parquet',format="parquet",
                       q.columns = c("Q.Value","Lib.Q.Value","Lib.PG.Q.Value"),
                       q.cutoffs = 0.01)
```

``` output
Length of q-value columns does not match with length of q-value cutoffs. Use q.cutoffs[1] for all columns.
```

:::::challenge

These quality filtering metrics are recommended in the [limpa vignette](https://bioconductor.org/packages/release/bioc/vignettes/limpa/inst/doc/limpa.html) for data searched in DIA-NN **with match-between-runs**.

What `q.columns` are recommended for data searched **without MBR**?

:::solution


``` r
y.peptide <- readDIANN(file='data/MBIntroToProteomics.parquet',format="parquet",
                       q.columns = c("Q.Value","Global.Q.Value","Global.PG.Q.Value"),
                       q.cutoffs = 0.01)
```
:::
:::::

Let's look at the structure of the data.


``` r
names(y.peptide)
```

``` output
[1] "E"     "genes"
```


``` r
head(y.peptide$E)
head(y.peptide$genes)
```

As you can see, our data has been imported as a limma EList object, with components **'E' (log2 peptide matrix)** and **'genes' (feature annotation)** containing information about each peptide.

<br>

### Load sample annotation

We also need to **import information about our samples** for later analyses.


``` r
samples_stool <- as.data.frame(readxl::read_excel('data/SampleAnnotation.xlsx'))
head(samples_stool)
```

``` output
  file_order Sample_id-new Batch Class
1         33       B1_29_c    B1   CDr
2         93       B1_87_e    B1   UCr
3         71       B1_66_e    B1   UCr
4          9       B1_08_b    B1   aCD
5         26       B1_22_c    B1   CDr
6         31       B1_27_c    B1   CDr
                                                                                                                                         file_name
1 D:\\Elmira\\IBD\\01-MSconverter\\00-DIA_NN input-3 batches\\mix.mzMLpeakpicking(batch 1,2)\\Mix\\20201016_UdSjfb_20190115_freshstool_067-29.mzML
2    D:\\Elmira\\IBD\\01-MSconverter\\00-DIA_NN input-3 batches\\mix.mzMLpeakpicking(batch 1,2)\\Mix\\20201016_UdSjfb_20190115_freshstool_187.mzML
3 D:\\Elmira\\IBD\\01-MSconverter\\00-DIA_NN input-3 batches\\mix.mzMLpeakpicking(batch 1,2)\\Mix\\20201016_UdSjfb_20190115_freshstool_143-66.mzML
4  D:\\Elmira\\IBD\\01-MSconverter\\00-DIA_NN input-3 batches\\mix.mzMLpeakpicking(batch 1,2)\\Mix\\20201016_UdSjfb_20190115_freshstool_019-8.mzML
5 D:\\Elmira\\IBD\\01-MSconverter\\00-DIA_NN input-3 batches\\mix.mzMLpeakpicking(batch 1,2)\\Mix\\20201016_UdSjfb_20190115_freshstool_053-22.mzML
6 D:\\Elmira\\IBD\\01-MSconverter\\00-DIA_NN input-3 batches\\mix.mzMLpeakpicking(batch 1,2)\\Mix\\20201016_UdSjfb_20190115_freshstool_063-27.mzML
```

We need to clean our sample annotation data to add an identifying column which **matches the sample column names in our peptide matrix.**

Given we are only using a subset of the main dataset in this tutorial, we can also filter the sample annotation dataframe to **only include rows pertaining to the samples in our subset.**


``` r
# Create sample_name column that matches peptide matrix column names
samples_stool$fixed_file_name <- gsub("\\\\", "/", samples_stool$file_name)
samples_stool$fixed_file_name <- basename(samples_stool$fixed_file_name)
samples_stool$sample_name <- gsub(".mzML", "",samples_stool$fixed_file_name)

# Filter to only samples present in our dataset
samples_stool <- samples_stool %>% filter(samples_stool$sample_name %in% colnames(y.peptide$E))
```

We can attach our sample annotation data to the EList object for later analyses.


``` r
# Ensure sample row order matches column order of the peptide matrix
sample_info <- samples_stool[samples_stool$sample_name %in% colnames(y.peptide$E),]
rownames(sample_info) <- sample_info$sample_name
y.peptide$E <- y.peptide$E[,rownames(sample_info)]

# Add experimental metadata to EList
y.peptide$targets <- sample_info[, c('Class', 'Batch')]
```

## Quality filtering

After reading in the data, it is common to conduct additional filtering to improve the quality of your data and assist with downstream interpretation of results.

- **Filtering out non-proteotypic peptides:** Removes peptides that belong to more than one protein.
- **Filtering out compound proteins:** Removes peptides belonging to compound protein groups consisting of multiple proteins separated by ";" delimiters.
- **Filtering out singleton peptides:** Removes peptides that are the only peptide belonging to that protein (i.e. keeps only proteins with multiple peptides).


``` r
y.peptide <- filterNonProteotypicPeptides(y.peptide)
y.peptide <- filterCompoundProteins(y.peptide)
y.peptide <- filterSingletonPeptides(y.peptide, min.n.peptides = 2)
```

These filtering steps are optional; if they are not run, `limpa` will still be able to conduct protein quantification. For more information about these peptide-level filters and whether they are appropriate for your experiment, please consult the [limpa documentation](https://bioconductor.org/packages/release/bioc/html/limpa.html).

## Removing contaminants

Contaminant background signals can influence mass spectrometry-based proteomics data. It is essential to **remove contaminant proteins** from the sample preparation process to avoid confounding results.

[Frankenfield et al. (2022)](https://doi.org/10.1021/acs.jproteome.2c00145) created a new set of sample-type specific and universal contaminant FASTA files which they showed reduced false identifications, increased protein identifications and did not influence protein quantification for DIA workflows. These FASTA files updated upon widely used contaminant protein lists from MaxQuant and the common Repository of Adventitious Proteins (cRAP).

In the previous lesson, we incorporated the universal contaminants FASTA into our spectral library. This FASTA is structured such that all contaminant proteins are named beginning with **'Cont_'**, which means we can simply remove any proteins from our analysis containing this string text in their name.


``` r
# Remove contaminant proteins from the genes dataframe
y.peptide$genes <- y.peptide$genes %>% filter(!str_detect(Protein.Group, "Cont_"))

# Match peptides in the matrix to only those remaining in the dataframe
y.peptide$E <- y.peptide$E[rownames(y.peptide$genes), ]
```

:::::callout
# More about contaminants

You can **download the contaminant FASTA files from [GitHub](https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics).**

If you are analysing data that has been processed without a contaminant FASTA, you can also download an **Excel spreadsheet contaminant list** from the same link and use this to identify and remove contaminant proteins from your dataset.

For more information, see the original publication:

Frankenfield, A. M., Ni, J., Ahmed, M., & Hao, L. (2022). Protein Contaminants Matter: Building Universal Protein Contaminant Libraries for DDA and DIA Proteomics. Journal of proteome research, 21(9), 2104–2113. https://doi.org/10.1021/acs.jproteome.2c00145

:::::

## Protein quantification

Missing data has been a longstanding issue impacting the accuracy and reproducibility of mass spectrometry-based proteomic analyses, complicating downstream analyses. Most established imputation methods are not designed to account for the specific attributes of proteomics data, and there is no consensus in the literature for how missing values should be managed in proteomics experiments. This presents an issue with using the protein values directly generated by platforms like DIA-NN, which contain many missing values.

In this tutorial, we are using the `dpcQuant()` method of protein quantification from the `limpa` package, or *Linear Models for Proteomics Data* [(Li, Cobbold & Smyth, 2025)](https://doi.org/10.1101/2025.04.28.651125). The **detection probability curve (DPC)** models the relationship between peptide intensity and missingness, showing that almost all missing values in label-free proteomics data are not missing at random and should be accounted for [(Li & Smyth, 2023)](https://doi.org/10.1093/bioinformatics/btad200). The `dpcQuant()` uses Bayesian statistical methods based on peptide intensities and the DPC coefficients for a particular dataset to generate a complete matrix of protein values without the need for imputation. This method outputs a linear model object suitable for downstream analysis with the `limma` package.

To learn more about the detection probability curve and the `limpa` package, please read the original publications.

### Using limpa and the DPC for protein quantification

First, we will generate and visualise our **detection probability curve**.

The DPC coefficients are related to the dataset and software used for processing. A steeper slope means there is more statistical information in the pattern of missing data. A slope closer to 0 means there is more data missing at random.



``` r
# Estimate the DPC for our data
dpcfit <- dpcON(y.peptide, robust=TRUE)

# Print the intercept and slope of our DPC
dpcfit$dpc
```

``` output
     beta0      beta1 
-7.2678650  0.6033267 
```

``` r
# Plot the DPC
plotDPC(dpcfit)
```

<img src="fig/03datacleaning-rendered-unnamed-chunk-11-1.png" style="display: block; margin: auto;" />

Next, we can use the DPC to generate our log2 protein values. This function requires the following input:

- a limma EList object, or a numeric matrix, of peptide-level log2 expression values, with samples as columns and peptides as rows (`y.peptide`)
- vector of protein IDs (`"Protein.Names"`)
- numeric vector providing the slope and intercept of the DPC (`dpcfit`)


``` r
y.protein <- dpcQuant(y.peptide, "Protein.Names", dpc=dpcfit)
```

``` output
Estimating hyperparameters ...
```

``` output
Quantifying proteins ...
```

``` output
Proteins: 280 Peptides: 2672
```

*Note: our dataset is quite small, so this process should only take a few seconds. For larger datasets, it may take minutes or hours for this command to run.*

:::callout
# Alternate versions of the DPC

The base `dpc()` assumes observed log-intensities are normally distributed.

A modified version, used in this tutorial, `dpcON(robust=TRUE)` downweighs outliers.

Alternatively, `dpcCN()` assumes the complete set of log-intensities are normally distributed. Using this method may underestimate the DPC slope if there is a non-small proportion of peptides missing by chance.

You can also choose your own slope by passing the argument `dpc.slope=` to `dpcQuant()` instead of `dpc=dpcfit`. Using this method, an intercept value will automatically be generated. You can check this value by using `estimateDPCIntercept()`. The authors recommend a slope of 0.6 - 0.8 is usually appropriate for data generated using DIA-NN.

You may wish to conduct further statistical analyses to determine which method is most appropriate for your experiment; however, the authors suggest that downstream analysis should not be significantly affected by differing DPC coefficients.

You can read more about the differences between each version of the DPC in the [paper](https://www.biorxiv.org/content/10.1101/2025.04.28.651125v1) or the [limpa manual](https://bioconductor.org/packages/release/bioc/manuals/limpa/man/limpa.pdf).
:::

:::::challenge

The limpa package also provides a function to impute missing peptide values and generate a complete peptide matrix using the DPC. Can you search the [limpa documentation](https://bioconductor.org/packages/release/bioc/html/limpa.html) to identify the command used to impute and quantify missing values in a rowwise manner?

:::solution

`dpcImpute()`

:::
::::::

### dpcQuant() output

We have now generated an EList object for our protein values. Let's investigate the other components of this object:


``` r
names(y.protein)
```

``` output
[1] "E"       "genes"   "other"   "targets" "dpc"    
```

``` r
names(y.protein$other)
```

``` output
[1] "n.observations" "standard.error"
```

``` r
names(y.protein$genes)
```

``` output
[1] "Protein.Group" "Protein.Names" "Genes"         "NPeptides"    
[5] "PropObs"      
```

| EList component   | Description                                                                         |
| :-----:           | ------------------                                                                        |
| **E**             | A complete matrix of log2 protein values with no missing data, with samples as columns and proteins as rows |
| **genes**         | Feature annotation dataframe, with columns for: UniProt IDs (Protein.Group), proteins (Protein.Names), gene names (Genes), number of peptides used to quantify each protein (NPeptides), and the proportion of possible observations used to generate each protein value (PropObs) |
| **other**         | List of two matrices with columns and rows corresponding to protein matrix E: number of peptides observed for each protein in each sample (n.observations), and the standard error associated with each protein value for each sample (standard.error) |
| **targets**       | Sample annotation dataframe |
| **dpc**           | DPC intercept and slope coefficient used for protein quantification |

:::::challenge

Compare the first few rows of the protein (`E`), `n.observations` and `standard.error` matrices. Take a look at the `NPeptides` and `PropObs` values for those proteins also. What do you notice about the relationship between these quality metrics?

:::solution

Proteins with smaller `NPeptides` and `PropObs` values, and smaller `n.observations` for individual samples, have higher standard error. 

The more peptide data available for a particular protein in a particular sample, the greater the confidence in the protein value, and therefore the lower the standard error.

If you are interested in a particular protein, it is important to look at these quality metrics for each of your samples.

:::
:::::

## Dimensionality Reduction and QC

Let's take a look at our dataset. The `plotMDSUsingSEs()` function, or *Multidimensional Scaling Plot of Gene Expression Profiles, Using Standard Errors*, is similar to `plotMDS()` in the `limma` package, but takes account of the standard errors generated by `dpcQuant()`. This is essentially a PCA plot for the newly generated protein values that takes into consideration their standard error.


``` r
plotMDSUsingSEs(y.protein)
```

<img src="fig/03datacleaning-rendered-unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

That plot is extremely messy and doesn't tell us much! Let's clean it up.


``` r
# Factorise and set colors for class and batch
Class <- factor(y.peptide$targets$Class)
Class.color <- Class
levels(Class.color) <- hcl.colors(nlevels(Class), palette = "cividis")

Batch <- factor(y.peptide$targets$Batch)
Batch.color <- Batch
levels(Batch.color) <- hcl.colors(nlevels(Batch), palette = "Set 2")

# Class visualisation
plotMDSUsingSEs(y.protein, pch=16, col=as.character(Class.color))
```

<img src="fig/03datacleaning-rendered-unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

``` r
# Batch visualisation
plotMDSUsingSEs(y.protein, pch=16, col=as.character(Batch.color))
```

<img src="fig/03datacleaning-rendered-unnamed-chunk-15-2.png" style="display: block; margin: auto;" />

``` r
# Class (labels) and batch (colors) visualisation
plotMDSUsingSEs(y.protein, labels=Class, col=as.character(Batch.color))
```

<img src="fig/03datacleaning-rendered-unnamed-chunk-15-3.png" style="display: block; margin: auto;" />

We can already see a distinction between the aCD and Ctrl groups, which will be explored further in the next lesson.

It also looks like there are some differences between the batches. **We should add batch as a covariate to our analyses in the next lesson.**

:::callout
# Managing batch effects

Explaining the different methods of identifying, accounting and correcting for batch effects in proteomics data is an entire tutorial itself! For further information, we recommend the [Melbourne Integrative Genomics Workshop: Managing batch effects in biological studies](https://github.com/melbintgen/Batch-effect-management?tab=readme-ov-file).

If you are concerned about batch effects in your data, University of Melbourne staff and students or affiliates can also [contact us](https://mdhs.unimelb.edu.au/melbournebioinformatics/training-and-support/expert-advice) for a free consultation.

:::

### Normalisation

Normalisation is another key element of the proteomics data cleaning workflow. Normalisation techniques aim to remove non-biological variation in the data, for example, variation due to sample loading or batch effects.

The default RT-dependent normalisation applied during DIA-NN processing is quite effective and often sufficient for DIA data; however, each experimental environment and techniques will have a different optimal approach. You may wish to apply an additional method of normalisation to your data.

To learn more about normalisation, keep an eye on the [Melbourne Bioinformatics Eventbrite page](https://www.eventbrite.com.au/o/melbourne-bioinformatics-13058846490) for a normalisation workshop coming in 2026.

:::spoiler
# More about normalisation

To learn more about different normalisation techniques for proteomics data, we recommend the following papers:

Peng, H., Wang, H., Kong, W., Li, J., & Goh, W. W. B. (2024). Optimizing differential expression analysis for proteomics data via high-performing rules and ensemble inference. Nature communications, 15(1), 3922. https://doi.org/10.1038/s41467-024-47899-w

Tseng, C. Y., Salguero, J. A., Breidenbach, J. D., Solomon, E., Sanders, C. K., Harvey, T., Thornhill, M. G., Palmisano, S. J., Sasiene, Z. J., Blackwell, B. R., McBride, E. M., Luchini, K. A., LeBrun, E. S., Alvarez, M., Mach, P. M., Rivera, E. S., & Glaros, T. G. (2025). Evaluation of normalization strategies for mass spectrometry-based multi-omics datasets. Metabolomics : Official journal of the Metabolomic Society, 21(4), 98. https://doi.org/10.1007/s11306-025-02297-1

Välikangas, T., Suomi, T., & Elo, L. L. (2018). A systematic evaluation of normalization methods in quantitative label-free proteomics. Briefings in bioinformatics, 19(1), 1–11. https://doi.org/10.1093/bib/bbw095

:::


Good work! We now have a clean dataset ready for analysis.


::::::::::::::::::::::::::::::::::::: keypoints 

- Key elements of the proteomics data cleaning workflow include quality filtering, removing contaminants, protein quantification, managing batch effects, and normalisation.
- We can use the *limpa* package in R to clean our peptide data and generate a complete, analysis-ready protein matrix with corresponding standard error values.

::::::::::::::::::::::::::::::::::::::::::::::::

