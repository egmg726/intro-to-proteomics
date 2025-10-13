---
title: 'Exploring the dataset on PRIDE'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- What does it mean to clean data?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Load and clean the data!

::::::::::::::::::::::::::::::::::::::::::::::::




## Load Libraries



``` r
library(limpa)
```

``` output
Loading required package: limma
```

``` r
library(EnhancedVolcano)
```

``` output
Loading required package: ggplot2
```

``` output
Loading required package: ggrepel
```

``` r
library(pheatmap)
library(rpx)
library(dplyr)
```

``` output

Attaching package: 'dplyr'
```

``` output
The following objects are masked from 'package:stats':

    filter, lag
```

``` output
The following objects are masked from 'package:base':

    intersect, setdiff, setequal, union
```

## The dataset we're using

The data comes from a study done exploring a non-invasive diagnostic approach for Inflammatory Bowel Disease (IBD) using SWATH Mass Spectrometry and machine learning to analyse the stool proteome. Traditional methods like colonoscopy are invasive, prompting the need for alternatives. The researchers collected 123 stool samples and identified 48 differentially expressed proteins, narrowing them down to 7 key proteins using feature selection. A Support Vector Machine (SVM) model was developed, achieving 96% sensitivity and 76% specificity in distinguishing active IBD patients from symptomatic non-IBD patients. This approach demonstrates the potential for accurate, non-invasive IBD diagnosis, improving patient management and reducing the need for invasive procedures. 


Shajari, E.; Gagné, D.; Malick, M.; Roy, P.; Noël, J.-F.; Gagnon, H.; Brunet, M.A.; Delisle, M.; Boisvert, F.-M.; Beaulieu, J.-F. Application of SWATH Mass Spectrometry and Machine Learning in the Diagnosis of Inflammatory Bowel Disease Based on the Stool Proteome. Biomedicines 2024, 12, 333. https://doi.org/10.3390/biomedicines12020333.


:::::: challenge

Can you see where to find the publicly available data in the paper?


::::::


## Explore the dataset with `rpx`

PRIDE is part of the **ProteomeXchange Consortium**, which provides open access to proteomics data.

For more background on navigating PRIDE data, you can look through this course: [PRIDE Quick Tour (EMBL-EBI Training)](https://www.ebi.ac.uk/training/online/courses/pride-quick-tour/)

The `rpx` package gives us programmatic access to the ProteomeXchange without having to leave R.


``` r
px_id <- 'PXD047585'
px <- PXDataset(px_id)
```

``` output
Querying ProteomeXchange for PXD047585.
```

``` r
# Dataset title
pxtitle(px)
```

``` output
[1] "Application of SWATH Mass Spectrometry and Machine Learning in Diagnosis of Inflammatory Bowel Disease Based on Stool Proteome"
```

``` r
# Link to the dataset on PRIDE
pxurl(px)
```

``` output
[1] "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/02/PXD047585"
```

``` r
# Taxonomy (organism) information
pxtax(px)
```

``` output
[1] "Homo sapiens (human)"
```

## Exploring Dataset Files

Let’s see which files are included in this dataset.


``` r
# Retrieve list of dataset files
px_files <- pxfiles(px)
```

``` output
Project PXD047585 files (255):
 [remote] 20201016_UdSjfb_20190115_freshstool_003.mzML
 [remote] 20201016_UdSjfb_20190115_freshstool_003.wiff
 [remote] 20201016_UdSjfb_20190115_freshstool_003.wiff.scan
 [remote] 20201016_UdSjfb_20190115_freshstool_005.mzML
 [remote] 20201016_UdSjfb_20190115_freshstool_005.wiff
 [remote] 20201016_UdSjfb_20190115_freshstool_005.wiff.scan
 [remote] 20201016_UdSjfb_20190115_freshstool_007-3.mzML
 [remote] 20201016_UdSjfb_20190115_freshstool_007.wiff
 [remote] 20201016_UdSjfb_20190115_freshstool_007.wiff.scan
 [remote] 20201016_UdSjfb_20190115_freshstool_011.mzML
 ...
```

``` r
# Display the first few files
head(px_files)
```

``` output
[1] "20201016_UdSjfb_20190115_freshstool_003.mzML"     
[2] "20201016_UdSjfb_20190115_freshstool_003.wiff"     
[3] "20201016_UdSjfb_20190115_freshstool_003.wiff.scan"
[4] "20201016_UdSjfb_20190115_freshstool_005.mzML"     
[5] "20201016_UdSjfb_20190115_freshstool_005.wiff"     
[6] "20201016_UdSjfb_20190115_freshstool_005.wiff.scan"
```



PRIDE datasets often include a mix of raw data, search results, processed outputs, and metadata.

To understand the data composition, we’ll examine file extensions.



``` r
# Count the number of files by extension
table(sapply(strsplit(px_files, "\\."), tail, 1))
```

``` output

    fas    mzML     pdf    scan speclib     tsv     txt    wiff    xlsx 
      1      78       1      78       2      13       1      78       3 
```

Alternatively, we can get the same result with the tools package:


``` r
library(tools)
table(file_ext(px_files))
```

``` output

    fas    mzML     pdf    scan speclib     tsv     txt    wiff    xlsx 
      1      78       1      78       2      13       1      78       3 
```


Understanding File Formats

Each file extension represents a specific proteomics data type.
For example, .raw may indicate vendor instrument output and .mzML indicates an open format for MS data

For more information:

[PRIDE File Formats Documentation](https://www.ebi.ac.uk/pride/markdownpage/pridefileformats)

[Paper on maass spec file formats (PMC3518119)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3518119/)



:::::: challenge

Try visiting the dataset page on the PRIDE website and compare:

The number of files listed there versus those returned by rpx.

The types of files (e.g., .raw, .mzML, .txt, etc.)

Questions to consider:

Why might the PRIDE web interface and R output differ?

What are the advantages of using programmatic access via rpx?


::::::


Note: add solution



::::::::::::::::::::::::::::::::::::: keypoints 

- Garbage in, garbage out!

::::::::::::::::::::::::::::::::::::::::::::::::

