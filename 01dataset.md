---
title: 'Exploring the dataset'
teaching: 10
exercises: 5
---

:::::::::::::::::::::::::::::::::::::: questions 

- What is the data we are using in this tutorial?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Understand the dataset we are using in this tutorial.
- Learn how to download publicly available proteomics data from within R.

::::::::::::::::::::::::::::::::::::::::::::::::

## The data

The data we are using in this tutorial comes from a study done exploring a non-invasive diagnostic approach for Inflammatory Bowel Disease (IBD) using mass spectrometry and machine learning to analyse the stool proteome. Traditional diagnostic methods like colonoscopies are invasive, prompting the need for alternatives. The researchers collected 123 stool samples and identified 48 differentially expressed proteins, narrowing them down to 7 key proteins using feature selection. A Support Vector Machine (SVM) model was developed, achieving 96% sensitivity and 76% specificity in distinguishing active IBD patients from symptomatic non-IBD patients. This approach demonstrates the potential for accurate, non-invasive IBD diagnosis, improving patient management and reducing the need for invasive procedures.

**For the purpose of this tutorial**, we will be using a subset of the total data presented in the paper, including:

- 11 control (Ctrl) samples

    - 6 from batch one + 5 from batch three
    
- 9 active Crohn’s disease (aCD) samples

    - 5 from batch one + 4 from batch three

**Reference**
Shajari, E.; Gagné, D.; Malick, M.; Roy, P.; Noël, J.-F.; Gagnon, H.; Brunet, M.A.; Delisle, M.; Boisvert, F.-M.; Beaulieu, J.-F. Application of SWATH Mass Spectrometry and Machine Learning in the Diagnosis of Inflammatory Bowel Disease Based on the Stool Proteome. Biomedicines 2024, 12, 333. https://doi.org/10.3390/biomedicines12020333.


:::::: challenge

Can you see where to find the publicly available data in the paper?

::::::

## Loading the data

Let's start by loading the packages required to access data stored on PRIDE. 

::: spoiler

## Learn more about PRIDE

PRIDE is part of the **ProteomeXchange Consortium**, which provides open access to proteomics data.

For more background on navigating PRIDE data, you can look through this course: [PRIDE Quick Tour (EMBL-EBI Training)](https://www.ebi.ac.uk/training/online/courses/pride-quick-tour/)

:::


``` r
library(rpx)
library(dplyr)
```

The `rpx` package gives us programmatic access to the ProteomeXchange without having to leave R.


``` r
# Set the dataset identifier and load
px_id <- 'PXD047585' 
px <- PXDataset(px_id)
```



We can then perform several functions to extract the contents of the `px` object:


``` r
# Print the dataset title
pxtitle(px)
```

``` output
[1] "Application of SWATH Mass Spectrometry and Machine Learning in Diagnosis of Inflammatory Bowel Disease Based on Stool Proteome"
```

``` r
# Print a link to the dataset on PRIDE
pxurl(px)
```

``` output
[1] "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2024/02/PXD047585"
```

``` r
# Print the taxonomy (organism) information for the dataset
pxtax(px)
```

``` output
[1] "Homo sapiens (human)"
```


### Exploring dataset files

Let’s see which files are included in this dataset.



``` r
# Retrieve a list of dataset files
px_files <- pxfiles(px)
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

Each file extension represents a specific proteomics data type. 
For example, `.raw` may indicate vendor instrument output and `.mzML` indicates an open format for MS data.

::: spoiler

## Learn more about proteomics file formats

For more information about different proteomics file formats, check out:

- [PRIDE File Formats Documentation](https://www.ebi.ac.uk/pride/markdownpage/pridefileformats)

- [Paper on mass spec file formats (PMC3518119)](https://pmc.ncbi.nlm.nih.gov/articles/PMC3518119/)

:::

:::::: challenge

Try visiting the dataset page on the PRIDE website and compare:

- The number of files listed there versus those returned by rpx.

- The types of files (e.g., .raw, .mzML, .txt, etc.)

Questions to consider:

- Why might the PRIDE web interface and R output differ?

- What are the advantages of using programmatic access via rpx?

::: solution

Programmatic access allows you to automate each step of your analysis, making your workflow fully reproducible.

:::
::::::

### Downloading required files

Good news! We have already downloaded and processed the data for you, as described in the next lesson. 
For future reference, you can use the **Curl** package to download files as demonstrated below.


``` r
library(curl)

curl::curl_download(url='https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/02/PXD047585/SampleAnnotation.xlsx',destfile = 'data/SampleAnnotation.xlsx')
```

- We retrieved the `url` in an earlier step using `pxurl()`. Note:
    
    - We have replaced `ftp://` at the start with `https://`
    - We have specified the name of the file we wish to download (`SampleAnnotation.xlsx`), from the list we sourced earlier using `pxfiles()`

- We set `destfile` to our desired name and location for the downloaded file.


::::::::::::::::::::::::::::::::::::: keypoints 

- The dataset used in this tutorial includes stool samples used to investigate non-invasive diagnostic methods for Inflammatory Bowel Disease.
- We can use the `rpx` package to download and explore data from PRIDE or the ProteomeXchange Consortium without leaving R.

::::::::::::::::::::::::::::::::::::::::::::::::

