---
title: 'Data cleaning using R'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- What is DE analysis

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- To do DE analysis

::::::::::::::::::::::::::::::::::::::::::::::::

## Load processed data and sample annotation 

After processing raw proteomics data, further cleaning and analysis can be conducted using R, Python, or by importing data to other platforms like **Perseus** or **Skyline**.

In this workshop, we will be teaching you one method of cleaning and analysing proteomics data using R. To follow along the remainder of this workshop, you should **open a new R script in RStudio** and **run the example code provided.**

:::::callout

# The limpa package

In this workshop, we are using the *limpa* package to load and clean our data. This package was published in 2025 by the [Smyth Lab](https://github.com/SmythLab/limpa), based at the Walter and Eliza Hall Institute of Medical Research (WEHI). 

For background information, see the original publications:

Li, M., Cobbold, S. A., & Smyth, G. K. (2025). Quantification and differential analysis of mass spectrometry proteomics data with probabilistic recovery of information from missing values. bioRxiv, 2025.2004.2028.651125. https://doi.org/10.1101/2025.04.28.651125 

Li, M., & Smyth, G. K. (2023). Neither random nor censored: estimating intensity-dependent probabilities for missing values in label-free proteomics. Bioinformatics, 39(5). https://doi.org/10.1093/bioinformatics/btad200 

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

Load the peptide-level data processed by DIA-NN (`.parquet` file), filtering observations by quality metrics as recommended in the [limpa vignette](https://bioconductor.org/packages/release/bioc/vignettes/limpa/inst/doc/limpa.html).


``` r
y.peptide <- readDIANN(file='data/MBIntroToProteomics.parquet',format="parquet",
                       q.columns = c("Q.Value","Lib.Q.Value","Lib.PG.Q.Value"),
                       q.cutoffs = 0.01)
```

``` error
Error in readDIANN(file = "data/MBIntroToProteomics.parquet", format = "parquet", : could not find function "readDIANN"
```

:::::challenge

These quality filtering metrics are recommended in the [limpa vignette](https://bioconductor.org/packages/release/bioc/vignettes/limpa/inst/doc/limpa.html) for data searched in DIA-NN **with match-between-runs**.

What q.columns are recommended for data searched **without MBR**?

:::solution


``` r
y.peptide <- readDIANN(file='data/MBIntroToProteomics.parquet',format="parquet",
                       q.columns = c("Q.Value","Global.Q.Value","Global.PG.Q.Value"),
                       q.cutoffs = 0.01)
```

``` error
Error in readDIANN(file = "data/MBIntroToProteomics.parquet", format = "parquet", : could not find function "readDIANN"
```
:::
:::::

Let's look at the structure of the data.


``` r
names(y.peptide)
```

``` error
Error: object 'y.peptide' not found
```

``` r
head(y.peptide$E)
```

``` error
Error: object 'y.peptide' not found
```

``` r
head(y.peptide$genes)
```

``` error
Error: object 'y.peptide' not found
```

As you can see, our data has been imported as a limma EList object, with components **'E' (log2 peptide matrix)** and **'genes' (feature annotation)** containing information about each peptide.

### Load sample data

We also need to import information about our samples for later analyses.


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

We need to clean our sample annotation data to add an identifying column which matches the sample column names in our peptide matrix.

Given we are only using a subset of the main dataset in this workshop, we can also filter the sample annotation dataframe to only include rows pertaining to the samples in our subset.


``` r
# Create sample_name column that matches peptide matrix column names
samples_stool$fixed_file_name <- gsub("\\\\", "/", samples_stool$file_name)
samples_stool$fixed_file_name <- basename(samples_stool$fixed_file_name)
samples_stool$sample_name <- gsub(".mzML", "",samples_stool$fixed_file_name)

# Filter to only samples present in our dataset
samples_stool <- samples_stool %>% filter(samples_stool$sample_name %in% colnames(y.peptide$E))
```

``` error
Error in samples_stool %>% filter(samples_stool$sample_name %in% colnames(y.peptide$E)): could not find function "%>%"
```

We can attach our sample annotation data to the EList object for later analyses.


``` r
# Ensure sample row order matches column order of the peptide matrix
sample_info <- samples_stool[samples_stool$sample_name %in% colnames(y.peptide$E),]
```

``` error
Error: object 'y.peptide' not found
```

``` r
rownames(sample_info) <- sample_info$sample_name
```

``` error
Error: object 'sample_info' not found
```

``` r
y.peptide$E <- y.peptide$E[,rownames(sample_info)]
```

``` error
Error: object 'y.peptide' not found
```

``` r
# Add experimental metadata to EList
y.peptide$targets <- sample_info[, c('Class', 'Batch')]
```

``` error
Error: object 'sample_info' not found
```

``` r
# Set colors for later plotting
Class <- factor(y.peptide$targets$Class)
```

``` error
Error: object 'y.peptide' not found
```

``` r
Class.color <- Class
```

``` error
Error: object 'Class' not found
```

``` r
levels(Class.color) <- hcl.colors(nlevels(Class), palette = "Set 2")
```

``` error
Error: object 'Class' not found
```

``` r
Batch <- factor(y.peptide$targets$Batch)
```

``` error
Error: object 'y.peptide' not found
```

``` r
Batch.color <- Batch
```

``` error
Error: object 'Batch' not found
```

``` r
levels(Batch.color) <- hcl.colors(nlevels(Batch), palette = "cividis")
```

``` error
Error: object 'Batch' not found
```


## Quality filtering

After reading in the data, it is common to conduct additional filtering to improve the quality of your data and assist with downstream interpretation of results.

- **Filtering out non-proteotypic peptides:** Removes peptides that belong to more than one protein.
- **Filtering out compound proteins:** Removes peptides belonging to compound protein groups consisting of multiple proteins separated by ";" delimiters.
- **Filtering out singleton peptides:** Removes peptides that are the only peptide belonging to that protein (i.e. keeps only proteins with multiple peptides).


``` r
y.peptide <- filterNonProteotypicPeptides(y.peptide)
```

``` error
Error in filterNonProteotypicPeptides(y.peptide): could not find function "filterNonProteotypicPeptides"
```

``` r
y.peptide <- filterCompoundProteins(y.peptide)
```

``` error
Error in filterCompoundProteins(y.peptide): could not find function "filterCompoundProteins"
```

``` r
y.peptide <- filterSingletonPeptides(y.peptide, min.n.peptides = 2)
```

``` error
Error in filterSingletonPeptides(y.peptide, min.n.peptides = 2): could not find function "filterSingletonPeptides"
```

These filtering steps are optional; limpa will still work if they are not run. For more information about these peptide-level filters and whether they are appropriate for your experiment, please consult the limpa documentation.

## Removing contaminants

Contaminant background signals can influence mass spectrometry-based proteomics data. It is essential to remove contaminant proteins from the sample preparation process to avoid confounding results.

Frankenfield et al. (2022) updated widely used contaminant protein lists from MaxQuant and the common Repository of Adventitious Proteins (cRAP) to create a new set of sample-type specific and universal contaminant FASTA files which reduced false identifications, increased protein identifications and did not influence protein quantification for DIA workflows. 

In the previous lesson, we incorporated the universal contaminants fasta into our spectral library. This fasta is structured such that all contaminant proteins are named beginning with **'Cont_'**, which means we can simply remove any proteins from our analysis containing this string text.


``` r
# Remove contaminant proteins from the genes dataframe
y.peptide$genes <- y.peptide$genes %>% filter(!str_detect(Protein.Group, "Cont_"))
```

``` error
Error in y.peptide$genes %>% filter(!str_detect(Protein.Group, "Cont_")): could not find function "%>%"
```

``` r
# Match peptides in the matrix to only those remaining in the dataframe
y.peptide$E <- y.peptide$E[rownames(y.peptide$genes), ]
```

``` error
Error: object 'y.peptide' not found
```

:::::callout
# More about contaminants

You can **download the contaminant fasta files from [GitHub](https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics).**

If you are analysing data that has been processed without the contaminant fasta, you can also download an **Excel spreadsheet contaminant list** from the same link and use this to identify and remove contaminant proteins from your dataset.

For more information, see the original publication:

Frankenfield AM, Ni J, Ahmed M, Hao L. Protein contaminants matter: building universal protein contaminant libraries for DDA and DIA proteomics. Journal of proteome research. 2022 Jul 6;21(9):2104-13. https://pubs.acs.org/doi/full/10.1021/acs.jproteome.2c00145

:::::

## Protein quantification

Missing data has been a longstanding issue impacting the accuracy and reproducibility of mass spectrometry-based proteomic analyses, complicating downstream analyses. Most imputation methods are not designed to account for the specific attributes of proteomics data, and there is no consensus in the literature for how missing values should be managed in proteomics experiments.

Li & Smyth (2023) published the detection probability curve (DPC), which models the relationship between peptide intensity and missingness. They show that almost all missing values in label-free proteomics data are not missing at random and should be accounted for.

Li, Cobbold & Smyth (2025) present a solution to this problem through *limpa* package, or *Linear Models for Proteomics Data*. Their method uses Bayesian statistical methods to generate a complete matrix of protein values without the need for imputation, based on peptide intensities and the DPC coefficients. This method outputs a linear model object suitable for downstream analysis with the *limma* package.

To learn more about the detection probability curve and the limpa package, see the original publications referenced at the beginning of this lesson.

### Using limpa and the DPC for protein quantification

First, we will generate and visualise our detection probability curve.

The DPC coefficients are related to the dataset and software used for processing. A steeper slope means there is more statistical information in the pattern of missing data. A slope closer to 0 means there is more data missing at random.



``` r
dpcfit <- dpcON(y.peptide, robust=TRUE)
```

``` error
Error in dpcON(y.peptide, robust = TRUE): could not find function "dpcON"
```

``` r
# Print the intercept and slope of our DPC
dpcfit$dpc
```

``` error
Error: object 'dpcfit' not found
```

``` r
# Plot the DPC
plotDPC(dpcfit)
```

``` error
Error in plotDPC(dpcfit): could not find function "plotDPC"
```

Now we can use the DPC to generate our log2 protein values. This function requires the following input:

- limma EList object, or a numeric matrix, of peptide-level log2 expression values, with samples as columns and peptides as rows
- vector of protein IDs
- numeric vector providing the slope and intercept of the DPC


``` r
y.protein <- dpcQuant(y.peptide, "Protein.Names", dpc=dpcfit)
```

``` error
Error in dpcQuant(y.peptide, "Protein.Names", dpc = dpcfit): could not find function "dpcQuant"
```

*Note: our dataset is quite small, so this process should only take a few seconds. For larger datasets, it may take minutes or hours for this command to run.*

:::callout
# Alternate versions of the DPC

The base `dpc()` assumes observed log-intensities are normally distributed.

A modified version `dpcON(robust=TRUE)` downweighs outliers.

Alternatively, `dpcCN()` assumes the complete set of log-intensities are normally distributed. Using this method may underestimate the DPC slope if there is a non-small proportion of peptides missing by chance.

You can also choose your own slope by passing the argument `dpc.slope=` to `dpcQuant()`. Using this method, an intercept value will automatically be generated. You can check this value by using `estimateDPCIntercept()`. The authors recommend a slope of 0.6 - 0.8 is usually appropriate for data generated using DIA-NN.

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

``` error
Error: object 'y.protein' not found
```

``` r
names(y.protein$other)
```

``` error
Error: object 'y.protein' not found
```

``` r
names(y.protein$genes)
```

``` error
Error: object 'y.protein' not found
```

| EList Component     | Description                                                                         |
| ------              | ------------------                                                                        |
| **E**               | A complete matrix of log2 protein values with no missing data, with samples as columns and proteins as rows |
| **genes**           | Feature annotation dataframe, including UniProt IDs (Protein.Group), proteins (Protein.Names), gene names (Genes), number of peptides used to quantify each protein (NPeptides), and the proportion of possible observations used to generate each protein value (PropObs) |
| **other**           | Two matrices with columns and rows corresponding to E: number of peptides observed in each sample for each protein (n.observations) and the standard error associated with each protein value for each sample (standard.error) |
| **targets**         | Sample annotation dataframe |
| **dpc**             | DPC intercept and slope coefficient used for protein quantification |

:::::challenge

Compare the first few rows of the protein (E), n.observations and standard.error matrices. Take a look at the NPeptides and PropObs values for those proteins also. What do you notice about the relationship between these values?

:::solution

Proteins with smaller NPeptides and PropObs values, and smaller n.observations for individual samples, have higher standard error. 

The more peptide data available for a particular protein in a particular sample, the greater the confidence in the protein value, and therefore the lower the standard error.

If you are interested in a particular protein, it is important to look at these quality metrics for each of your samples.

:::
:::::



### Dimensionality Reduction and QC






That plot is extremely messy and doesn't tell us much! Let's clean it up.





``` r
plotMDSUsingSEs(y.protein)
```

``` error
Error in plotMDSUsingSEs(y.protein): could not find function "plotMDSUsingSEs"
```

``` r
# Class visualisation

plotMDSUsingSEs(y.protein, pch=16, col=as.character(Class.color))
```

``` error
Error in plotMDSUsingSEs(y.protein, pch = 16, col = as.character(Class.color)): could not find function "plotMDSUsingSEs"
```

Hmmm, not a lot of differences. Let's see if there are any batch effects.


``` r
# Batch visualisation

plotMDSUsingSEs(y.protein, pch=16, col=as.character(Batch.color))
```

``` error
Error in plotMDSUsingSEs(y.protein, pch = 16, col = as.character(Batch.color)): could not find function "plotMDSUsingSEs"
```

It looks like there are some differences between the batches. Let's remove the batch effects and see how it looks.



``` r
# Correct for batch effects
y.protein.rbe <- removeBatchEffect(y.protein,batch = y.protein$targets$Batch)
```

``` error
Error in removeBatchEffect(y.protein, batch = y.protein$targets$Batch): could not find function "removeBatchEffect"
```

``` r
plotMDS(y.protein.rbe, pch=16, col=as.character(Batch.color))
```

``` error
Error in plotMDS(y.protein.rbe, pch = 16, col = as.character(Batch.color)): could not find function "plotMDS"
```

## Subset Comparison: Ctrl vs aUC


``` r
cdr_bool <- y.peptide$targets$Class %in% c('Ctrl','aUC')
```

``` error
Error: object 'y.peptide' not found
```

``` r
plotMDSUsingSEs(y.protein[, cdr_bool], 
                pch=16, col=as.character(Class.color[cdr_bool]))
```

``` error
Error in plotMDSUsingSEs(y.protein[, cdr_bool], pch = 16, col = as.character(Class.color[cdr_bool])): could not find function "plotMDSUsingSEs"
```

``` r
plotMDS(y.protein[, cdr_bool], 
         pch=16, col=as.character(Class.color[cdr_bool]))
```

``` error
Error in plotMDS(y.protein[, cdr_bool], pch = 16, col = as.character(Class.color[cdr_bool])): could not find function "plotMDS"
```

``` r
plotMDS(y.protein.rbe[, cdr_bool], 
         pch=16, col=as.character(Class.color[cdr_bool]))
```

``` error
Error in plotMDS(y.protein.rbe[, cdr_bool], pch = 16, col = as.character(Class.color[cdr_bool])): could not find function "plotMDS"
```

We should add batch as a covariate to our analysis.

## Differential Expression Analysis

We now fit a linear model that includes both Class and Batch effects.


``` r
design <- model.matrix(~0 + Class + Batch)
```

``` error
Error in eval(predvars, data, env): object 'Class' not found
```

``` r
fit <- dpcDE(y.protein, design, plot = TRUE)
```

``` error
Error in dpcDE(y.protein, design, plot = TRUE): could not find function "dpcDE"
```

``` r
fit <- eBayes(fit)
```

``` error
Error in eBayes(fit): could not find function "eBayes"
```


``` r
topTable(fit, coef = "ClassCtrl")
```

``` error
Error in topTable(fit, coef = "ClassCtrl"): could not find function "topTable"
```

## Create and Explore Contrasts

We can define custom contrasts (e.g., comparing aCD vs Ctrl) and explore differential proteins.


``` r
# Define contrasts
contrasts_fit <- makeContrasts(Ctrl_aCD = ClassaCD - ClassCtrl, 
                               levels = colnames(design))
```

``` error
Error in makeContrasts(Ctrl_aCD = ClassaCD - ClassCtrl, levels = colnames(design)): could not find function "makeContrasts"
```

``` r
# Apply contrasts
fit2 <- contrasts.fit(fit, contrasts = contrasts_fit)
```

``` error
Error in contrasts.fit(fit, contrasts = contrasts_fit): could not find function "contrasts.fit"
```

``` r
fit2 <- eBayes(fit2)
```

``` error
Error in eBayes(fit2): could not find function "eBayes"
```

``` r
# Summaries
topTable(fit2, coef = 1)
```

``` error
Error in topTable(fit2, coef = 1): could not find function "topTable"
```


``` r
# Compare multiple coefficients
topTable(fit, coef = 1)
```

``` error
Error in topTable(fit, coef = 1): could not find function "topTable"
```

``` r
topTable(fit, coef = 3)
```

``` error
Error in topTable(fit, coef = 3): could not find function "topTable"
```

``` r
topTable(fit, coef = 4)
```

``` error
Error in topTable(fit, coef = 4): could not find function "topTable"
```

``` r
topTable(fit, coef = 6)
```

``` error
Error in topTable(fit, coef = 6): could not find function "topTable"
```

``` r
# Example diagnostic plots
plotMD(fit, coef = 1)
```

``` error
Error in plotMD(fit, coef = 1): could not find function "plotMD"
```

``` r
plotMD(fit, coef = 2)
```

``` error
Error in plotMD(fit, coef = 2): could not find function "plotMD"
```

## Example: aUC vs Ctrl Comparison


``` r
results <- topTable(fit, coef = 3, number = Inf)
```

``` error
Error in topTable(fit, coef = 3, number = Inf): could not find function "topTable"
```

``` r
# Visualize a specific protein
plotProtein(y.protein, "S10A9_HUMAN", col = as.character(Class.color))
```

``` error
Error in plotProtein(y.protein, "S10A9_HUMAN", col = as.character(Class.color)): could not find function "plotProtein"
```

``` r
legend('topleft', legend = levels(Class), fill = levels(Class.color))
```

``` error
Error: object 'Class' not found
```

## Visualize Top Significant Proteins

We identify the top 50 most significant and variable proteins, then visualize via a clustered heatmap.


``` r
sig_proteins <- results %>%
  filter(adj.P.Val < 0.05) %>%
  top_n(50, wt = abs(logFC)) %>%
  pull(Protein.Names)
```

``` error
Error in results %>% filter(adj.P.Val < 0.05) %>% top_n(50, wt = abs(logFC)) %>% : could not find function "%>%"
```

``` r
expr_matrix <- y.protein$E[sig_proteins, ]
```

``` error
Error: object 'y.protein' not found
```

``` r
scaled_expr <- t(scale(t(expr_matrix)))
```

``` error
Error: object 'expr_matrix' not found
```

``` r
col_annotation <- data.frame(Class = Class,
                             row.names = colnames(y.protein$E))
```

``` error
Error: object 'y.protein' not found
```

``` r
pheatmap(scaled_expr,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = col_annotation)
```

``` error
Error in pheatmap(scaled_expr, cluster_rows = TRUE, cluster_cols = TRUE, : could not find function "pheatmap"
```


## Volcano plot


``` r
EnhancedVolcano(results,
                lab = results$Protein.Names,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 3.0,
                labSize = 3.0)
```

``` error
Error in EnhancedVolcano(results, lab = results$Protein.Names, x = "logFC", : could not find function "EnhancedVolcano"
```


## Network and Enrichment Analysis

### STRING Protein–Protein Interaction Network

Note: this part does not seem to be working


``` r
string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400)
mapped <- string_db$map(results, "Protein.Names", removeUnmappedRows = TRUE)
de_proteins <- mapped %>% filter(adj.P.Val < 0.05)

# Plot network of significant proteins
string_db$plot_network(de_proteins$STRING_id)
```

### Functional Enrichment (GO and KEGG)


``` r
de_proteins <- results %>% filter(adj.P.Val < 0.05)
```

``` error
Error in results %>% filter(adj.P.Val < 0.05): could not find function "%>%"
```

``` r
# Convert UniProt IDs to Entrez IDs for enrichment
converted <- bitr(de_proteins$Protein.Group,
                  fromType = "UNIPROT",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
```

``` error
Error in bitr(de_proteins$Protein.Group, fromType = "UNIPROT", toType = "ENTREZID", : could not find function "bitr"
```

``` r
# GO Biological Process enrichment
ego <- enrichGO(gene = converted$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)
```

``` error
Error in enrichGO(gene = converted$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", : could not find function "enrichGO"
```

``` r
dotplot(ego, showCategory = 10)
```

``` error
Error in dotplot(ego, showCategory = 10): could not find function "dotplot"
```

``` r
# KEGG pathway enrichment
ekegg <- enrichKEGG(gene = converted$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
```

``` error
Error in enrichKEGG(gene = converted$ENTREZID, organism = "hsa", pvalueCutoff = 0.05): could not find function "enrichKEGG"
```

``` r
dotplot(ekegg, showCategory = 10)
```

``` error
Error in dotplot(ekegg, showCategory = 10): could not find function "dotplot"
```


## Potential next steps

* Compare coefficients between classes as an exercise.

* Visualize selected proteins of interest across experimental groups.

* Extend enrichment analyses using Reactome or GSEA approaches.

* Integrate sample metadata for clinical covariates or longitudinal modeling.




::::::::::::::::::::::::::::::::::::: keypoints 

- Data!

::::::::::::::::::::::::::::::::::::::::::::::::

