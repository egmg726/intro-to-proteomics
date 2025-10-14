---
title: 'Data cleaning and DE analysis'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- What is DE analysis

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- To do DE analysis

::::::::::::::::::::::::::::::::::::::::::::::::



## Load processed data and additional files








This data was processed using DIA-NN in the method described in the paper.


``` r
# Getting the subset of the data

# curl::curl_download('https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/02/PXD047585/report.tsv','data/report.tsv')
# report_tsv <- read.csv('data/report.tsv', sep='\t')
# write.csv(stool_tsv_test,'~data/report.csv')
# parquetize::csv_to_parquet('data/report.csv',path_to_parquet = 'data/report.parquet')
# 
# report_full <- arrow::read_parquet('data/report.parquet')
# 
# proteins_keep <- c("RED2_HUMAN', 'MINT_HUMAN', 'Cationic', 'OVCH2_HUMAN', 'ITA4_HUMAN', 'Trypsin', 'G3PT_HUMAN', 'ZN451_HUMAN', 'DCD_HUMAN', 'CDHR2_HUMAN', 'DSG1_HUMAN', 'BRCA2_HUMAN', 'SIAE_HUMAN', 'CTRC_HUMAN', 'CEL3A_HUMAN', 'A1AT_HUMAN', 'CYTA_HUMAN', 'CADH1_HUMAN', 'TRY2_HUMAN', 'CEAM7_HUMAN', 'ITLN1_HUMAN', 'CEL3B_HUMAN', 'DPEP1_HUMAN', 'CUZD1_HUMAN', 'ZG16_HUMAN', 'ENPP7_HUMAN', 'HV601_HUMAN', 'SMR3B_HUMAN', 'PA21B_HUMAN', 'FCGBP_HUMAN', 'LEG4_HUMAN', 'RN170_HUMAN', 'IGLL5_HUMAN', 'PPBI_HUMAN', 'CLD34_HUMAN', 'CEAM5_HUMAN', 'CEAM6_HUMAN', 'MGA_HUMAN', 'TFF3_HUMAN', 'ENTK_HUMAN', 'IGKC_HUMAN', 'HV307_HUMAN', 'IGJ_HUMAN', 'CBPB1_HUMAN', 'KV401_HUMAN', 'CTRB2_HUMAN', 'CBPA1_HUMAN', 'AMYP_HUMAN', 'SODC_HUMAN', 'CEL2A_HUMAN")
# report_small <- report_full %>%
#   filter(Protein.Group %in% proteins_keep | runif(n()) < 0.10)
# arrow::write_parquet(report_small, "data/report_subset.parquet")
```





``` r
# Load peptide-level data processed by DIA-NN
y.peptide <- readDIANN(file='data/report_subset.parquet',format="parquet")

names(y.peptide)
```

``` output
[1] "E"     "genes"
```

## Annotating and filtering our dataset

In order to do differential expression analysis, we need information about what samples were used.


``` r
# The SampleAnnotation file was downloaded directly from PRIDE using this link>
# curl::curl_download(url='https://ftp.pride.ebi.ac.uk/pride/data/archive/2024/02/PXD047585/SampleAnnotation.xlsx',destfile = 'data/SampleAnnotation.xlsx')
```



``` r
# If you are on a Mac, this data needed to be cleaned up a bit with gsub if you download directly from PRIDE
samples_stool <- as.data.frame(readxl::read_excel('data/SampleAnnotation.xlsx'))
samples_stool$fixed_file_name <- gsub("\\\\", "/", samples_stool$file_name)
samples_stool$fixed_file_name <- basename(samples_stool$fixed_file_name)
samples_stool$sample_name <- gsub(".mzML", "",samples_stool$fixed_file_name)
```





The contaminant information is on GitHub and can be downloaded from their GitHub page.

More info available here: 
Frankenfield AM, Ni J, Ahmed M, Hao L. Protein contaminants matter: building universal protein contaminant libraries for DDA and DIA proteomics. Journal of proteome research. 2022 Jul 6;21(9):2104-13. https://pubs.acs.org/doi/full/10.1021/acs.jproteome.2c00145




``` r
# Code for downloading the contaminant libraries

# curl::curl_download('https://github.com/HaoGroup-ProtContLib/Protein-Contaminant-Libraries-for-DDA-and-DIA-Proteomics/raw/refs/heads/main/Universal%20protein%20contaminant%20FASTA/Contaminant%20protein%20information%20in%20the%20FASTA%20library%20and%20the%20potential%20source%20of%20contaminations.xlsx',destfile = 'data/Contaminants.xlsx')
```





``` r
contaminants <- as.data.frame(readxl::read_excel('data/Contaminants.xlsx'))
```

``` output
New names:
• `` -> `...2`
• `` -> `...3`
• `` -> `...4`
• `` -> `...5`
• `` -> `...6`
• `` -> `...7`
• `` -> `...8`
```

``` r
contaminants_colnames <- as.character(contaminants[1,])
contaminants <- contaminants[2:nrow(contaminants),]
colnames(contaminants) <- contaminants_colnames
```


``` r
# Check how many proteins match contaminants
table(y.peptide$genes$Protein.Group %in% contaminants$`Uniprot ID`)
```

``` output

FALSE  TRUE 
 6202   313 
```

``` r
# Filter out contaminants
y.peptide$genes <- dplyr::filter(y.peptide$genes, !(Protein.Group %in% contaminants$`Uniprot ID`))
y.peptide$E <- y.peptide$E[rownames(y.peptide$genes), ]

# Confirm removal
table(y.peptide$genes$Protein.Group %in% contaminants$`Uniprot ID`)
```

``` output

FALSE 
 6202 
```


### Align sample metadata



``` r
sample_info <- samples_stool[samples_stool$sample_name %in% colnames(y.peptide$E),]
rownames(sample_info) <- sample_info$sample_name
y.peptide$E <- y.peptide$E[,rownames(sample_info)]
```




``` r
# Add experimental metadata
y.peptide$targets <- sample_info[, c('Class', 'Batch')]

# Apply peptide-level filters
y.peptide <- filterNonProteotypicPeptides(y.peptide)
y.peptide <- filterCompoundProteins(y.peptide)
y.peptide <- filterSingletonPeptides(y.peptide, min.n.peptides = 2)
```


For more information on the peptide-level filters, please consult the limpa vignette.


### Dimensionality Reduction and QC


``` r
dpcfit <- dpc(y.peptide)
```

``` output
157 peptides are completely missing in all samples.
```

``` r
plotDPC(dpcfit)
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-11-1.png" style="display: block; margin: auto;" />


By protein:


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
Proteins: 406 Peptides: 4035
```

``` r
plotMDSUsingSEs(y.protein)
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-12-1.png" style="display: block; margin: auto;" />

That plot is extremely messy and doesn't tell us much! Let's clean it up.



``` r
# Class visualisation
Class <- factor(y.peptide$targets$Class)
levels(Class) <- c("Ctrl","aCD","aUC","CDr","UCr")
Class.color <- Class
levels(Class.color) <- c("pink","black",'grey','darkorange','darkgreen')
plotMDSUsingSEs(y.protein, pch=16, col=as.character(Class.color))
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-13-1.png" style="display: block; margin: auto;" />

Hmmm, not a lot of differences. Let's see if there are any batch effects.


``` r
# Batch visualisation
Batch <- factor(y.peptide$targets$Batch)
Batch.color <- Batch
levels(Batch.color) <- c("green","purple",'orange')
plotMDSUsingSEs(y.protein, pch=16, col=as.character(Batch.color))
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-14-1.png" style="display: block; margin: auto;" />

It looks like there are some differences between the batches. Let's remove the batch effects and see how it looks.



``` r
# Correct for batch effects
y.protein.rbe <- removeBatchEffect(y.protein,batch = y.protein$targets$Batch)
```

``` output
design matrix of interest not specified. Assuming a one-group experiment.
```

``` r
plotMDS(y.protein.rbe, pch=16, col=as.character(Batch.color))
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-15-1.png" style="display: block; margin: auto;" />

## Subset Comparison: Ctrl vs aUC


``` r
cdr_bool <- y.peptide$targets$Class %in% c('Ctrl','aUC')

plotMDSUsingSEs(y.protein[, cdr_bool], 
                pch=16, col=as.character(Class.color[cdr_bool]))
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-16-1.png" style="display: block; margin: auto;" />

``` r
plotMDS(y.protein[, cdr_bool], 
         pch=16, col=as.character(Class.color[cdr_bool]))
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-16-2.png" style="display: block; margin: auto;" />

``` r
plotMDS(y.protein.rbe[, cdr_bool], 
         pch=16, col=as.character(Class.color[cdr_bool]))
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-16-3.png" style="display: block; margin: auto;" />

We should add batch as a covariate to our analysis.

## Differential Expression Analysis

We now fit a linear model that includes both Class and Batch effects.


``` r
design <- model.matrix(~0 + Class + Batch)
fit <- dpcDE(y.protein, design, plot = TRUE)
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-17-1.png" style="display: block; margin: auto;" />

``` r
fit <- eBayes(fit)
```


``` r
topTable(fit, coef = "ClassCtrl")
```

``` output
                                 Protein.Group      Protein.Names
DPP2_HUMAN                              Q9UHL4         DPP2_HUMAN
Beta-galactosidase contam_sp|P00722|BGAL_ECOLI Beta-galactosidase
PSME2_HUMAN                             Q9UL46        PSME2_HUMAN
IGLC7_HUMAN                             A0M8Q6        IGLC7_HUMAN
GPA33_HUMAN                             Q99795        GPA33_HUMAN
MYO1D_HUMAN                             O94832        MYO1D_HUMAN
PRDX1_HUMAN                             Q06830        PRDX1_HUMAN
KCRB_HUMAN                              P12277         KCRB_HUMAN
HSPB1_HUMAN                             P04792        HSPB1_HUMAN
HMGB2_HUMAN                             P26583        HMGB2_HUMAN
                                Genes NPeptides     PropObs     logFC   AveExpr
DPP2_HUMAN                       DPP7         2 0.003225806 10.220764 10.221984
Beta-galactosidase Beta-galactosidase         2 0.006451613  9.924804  9.923592
PSME2_HUMAN                     PSME2         2 0.006451613 10.226636 10.228126
IGLC7_HUMAN                     IGLC7         2 0.006451613 10.048682 10.047295
GPA33_HUMAN                     GPA33         3 0.004301075  9.660197  9.662022
MYO1D_HUMAN                     MYO1D         2 0.009677419  9.731588  9.728144
PRDX1_HUMAN                     PRDX1         2 0.009677419 10.432375 10.436082
KCRB_HUMAN                        CKB         2 0.012903226 10.242352 10.244198
HSPB1_HUMAN                     HSPB1         2 0.006451613 10.491556 10.495906
HMGB2_HUMAN                     HMGB2         2 0.003225806 10.312915 10.316145
                          t P.Value adj.P.Val        B
DPP2_HUMAN         4877.269       0         0 850.6721
Beta-galactosidase 6395.855       0         0 841.4553
PSME2_HUMAN        6597.852       0         0 838.1056
IGLC7_HUMAN        3960.205       0         0 830.1724
GPA33_HUMAN        3963.989       0         0 827.7626
MYO1D_HUMAN        4831.810       0         0 826.2010
PRDX1_HUMAN        5230.944       0         0 816.8051
KCRB_HUMAN         4244.864       0         0 806.9717
HSPB1_HUMAN        3998.250       0         0 805.8680
HMGB2_HUMAN        2999.040       0         0 800.8352
```

## Create and Explore Contrasts

We can define custom contrasts (e.g., comparing aCD vs Ctrl) and explore differential proteins.


``` r
# Define contrasts
contrasts_fit <- makeContrasts(Ctrl_aCD = ClassaCD - ClassCtrl, 
                               levels = colnames(design))

# Apply contrasts
fit2 <- contrasts.fit(fit, contrasts = contrasts_fit)
fit2 <- eBayes(fit2)

# Summaries
topTable(fit2, coef = 1)
```

``` output
            Protein.Group Protein.Names    Genes NPeptides    PropObs
ALBU_HUMAN         P02768    ALBU_HUMAN      ALB        81 0.02692155
HBB_HUMAN          P68871     HBB_HUMAN      HBB        11 0.02697947
S10A9_HUMAN        P06702   S10A9_HUMAN   S100A9        30 0.05612903
DPP4_HUMAN         P27487    DPP4_HUMAN     DPP4        54 0.05627240
NEP_HUMAN          P08473     NEP_HUMAN      MME        49 0.05345622
TRFE_HUMAN         P02787    TRFE_HUMAN       TF        32 0.01108871
HEMO_HUMAN         P02790    HEMO_HUMAN      HPX         9 0.01505376
MUC2_HUMAN         Q02817    MUC2_HUMAN     MUC2        79 0.04752960
ANT3_HUMAN         P01008    ANT3_HUMAN SERPINC1        34 0.04724858
PRTN3_HUMAN        P24158   PRTN3_HUMAN    PRTN3        10 0.04709677
                  logFC  AveExpr         t      P.Value   adj.P.Val         B
ALBU_HUMAN   2.23879844 10.40131  4.478478 1.480449e-05 0.006010624  1.677531
HBB_HUMAN    0.36899249 11.35499  2.357978 1.966516e-02 0.998947663 -3.357132
S10A9_HUMAN -0.67889289 12.45692 -2.174909 3.120371e-02 0.998947663 -3.530783
DPP4_HUMAN   0.76076855 11.94575  2.030031 4.412074e-02 0.998947663 -3.771985
NEP_HUMAN    0.51280892 11.55182  1.850375 6.622674e-02 0.998947663 -4.036171
TRFE_HUMAN   0.21123633 10.12176  1.963051 5.148991e-02 0.998947663 -4.251168
HEMO_HUMAN   0.08091332 11.54930  1.974498 5.016131e-02 0.998947663 -4.262076
MUC2_HUMAN   0.66751339 11.82464  1.630523 1.050890e-01 0.998947663 -4.262562
ANT3_HUMAN  -0.31282593 10.81556 -1.740333 8.385128e-02 0.998947663 -4.525959
PRTN3_HUMAN -0.29640473 11.81750 -1.691011 9.291138e-02 0.998947663 -4.634509
```


``` r
# Compare multiple coefficients
topTable(fit, coef = 1)
```

``` output
                                 Protein.Group      Protein.Names
DPP2_HUMAN                              Q9UHL4         DPP2_HUMAN
Beta-galactosidase contam_sp|P00722|BGAL_ECOLI Beta-galactosidase
PSME2_HUMAN                             Q9UL46        PSME2_HUMAN
IGLC7_HUMAN                             A0M8Q6        IGLC7_HUMAN
GPA33_HUMAN                             Q99795        GPA33_HUMAN
MYO1D_HUMAN                             O94832        MYO1D_HUMAN
PRDX1_HUMAN                             Q06830        PRDX1_HUMAN
KCRB_HUMAN                              P12277         KCRB_HUMAN
HSPB1_HUMAN                             P04792        HSPB1_HUMAN
HMGB2_HUMAN                             P26583        HMGB2_HUMAN
                                Genes NPeptides     PropObs     logFC   AveExpr
DPP2_HUMAN                       DPP7         2 0.003225806 10.220764 10.221984
Beta-galactosidase Beta-galactosidase         2 0.006451613  9.924804  9.923592
PSME2_HUMAN                     PSME2         2 0.006451613 10.226636 10.228126
IGLC7_HUMAN                     IGLC7         2 0.006451613 10.048682 10.047295
GPA33_HUMAN                     GPA33         3 0.004301075  9.660197  9.662022
MYO1D_HUMAN                     MYO1D         2 0.009677419  9.731588  9.728144
PRDX1_HUMAN                     PRDX1         2 0.009677419 10.432375 10.436082
KCRB_HUMAN                        CKB         2 0.012903226 10.242352 10.244198
HSPB1_HUMAN                     HSPB1         2 0.006451613 10.491556 10.495906
HMGB2_HUMAN                     HMGB2         2 0.003225806 10.312915 10.316145
                          t P.Value adj.P.Val        B
DPP2_HUMAN         4877.269       0         0 850.6721
Beta-galactosidase 6395.855       0         0 841.4553
PSME2_HUMAN        6597.852       0         0 838.1056
IGLC7_HUMAN        3960.205       0         0 830.1724
GPA33_HUMAN        3963.989       0         0 827.7626
MYO1D_HUMAN        4831.810       0         0 826.2010
PRDX1_HUMAN        5230.944       0         0 816.8051
KCRB_HUMAN         4244.864       0         0 806.9717
HSPB1_HUMAN        3998.250       0         0 805.8680
HMGB2_HUMAN        2999.040       0         0 800.8352
```

``` r
topTable(fit, coef = 3)
```

``` output
                                 Protein.Group      Protein.Names
DPP2_HUMAN                              Q9UHL4         DPP2_HUMAN
Beta-galactosidase contam_sp|P00722|BGAL_ECOLI Beta-galactosidase
PSME2_HUMAN                             Q9UL46        PSME2_HUMAN
IGLC7_HUMAN                             A0M8Q6        IGLC7_HUMAN
GPA33_HUMAN                             Q99795        GPA33_HUMAN
MYO1D_HUMAN                             O94832        MYO1D_HUMAN
PRDX1_HUMAN                             Q06830        PRDX1_HUMAN
HSPB1_HUMAN                             P04792        HSPB1_HUMAN
SCMC1_HUMAN                             Q6NUK1        SCMC1_HUMAN
ACTN3_HUMAN                             Q08043        ACTN3_HUMAN
                                Genes NPeptides     PropObs     logFC   AveExpr
DPP2_HUMAN                       DPP7         2 0.003225806 10.220752 10.221984
Beta-galactosidase Beta-galactosidase         2 0.006451613  9.925068  9.923592
PSME2_HUMAN                     PSME2         2 0.006451613 10.226641 10.228126
IGLC7_HUMAN                     IGLC7         2 0.006451613 10.047665 10.047295
GPA33_HUMAN                     GPA33         3 0.004301075  9.661512  9.662022
MYO1D_HUMAN                     MYO1D         2 0.009677419  9.731603  9.728144
PRDX1_HUMAN                     PRDX1         2 0.009677419 10.432355 10.436082
HSPB1_HUMAN                     HSPB1         2 0.006451613 10.491542 10.495906
SCMC1_HUMAN                  SLC25A24         2 0.009677419 10.453710 10.457892
ACTN3_HUMAN                     ACTN3         2 0.006451613 10.381462 10.385495
                          t P.Value adj.P.Val        B
DPP2_HUMAN         6014.446       0         0 882.1088
Beta-galactosidase 7930.433       0         0 873.7109
PSME2_HUMAN        8137.383       0         0 869.5640
IGLC7_HUMAN        4972.912       0         0 864.3338
GPA33_HUMAN        4951.280       0         0 861.1143
MYO1D_HUMAN        5887.336       0         0 855.8380
PRDX1_HUMAN        6446.345       0         0 848.1433
HSPB1_HUMAN        4930.013       0         0 837.2906
SCMC1_HUMAN        5176.668       0         0 835.4328
ACTN3_HUMAN        4544.208       0         0 833.0233
```

``` r
topTable(fit, coef = 4)
```

``` output
                                 Protein.Group      Protein.Names
DPP2_HUMAN                              Q9UHL4         DPP2_HUMAN
Beta-galactosidase contam_sp|P00722|BGAL_ECOLI Beta-galactosidase
PSME2_HUMAN                             Q9UL46        PSME2_HUMAN
IGLC7_HUMAN                             A0M8Q6        IGLC7_HUMAN
GPA33_HUMAN                             Q99795        GPA33_HUMAN
MYO1D_HUMAN                             O94832        MYO1D_HUMAN
PRDX1_HUMAN                             Q06830        PRDX1_HUMAN
KCRB_HUMAN                              P12277         KCRB_HUMAN
HSPB1_HUMAN                             P04792        HSPB1_HUMAN
ACTN3_HUMAN                             Q08043        ACTN3_HUMAN
                                Genes NPeptides     PropObs     logFC   AveExpr
DPP2_HUMAN                       DPP7         2 0.003225806 10.222675 10.221984
Beta-galactosidase Beta-galactosidase         2 0.006451613  9.925129  9.923592
PSME2_HUMAN                     PSME2         2 0.006451613 10.226647 10.228126
IGLC7_HUMAN                     IGLC7         2 0.006451613 10.046443 10.047295
GPA33_HUMAN                     GPA33         3 0.004301075  9.663383  9.662022
MYO1D_HUMAN                     MYO1D         2 0.009677419  9.730722  9.728144
PRDX1_HUMAN                     PRDX1         2 0.009677419 10.433241 10.436082
KCRB_HUMAN                        CKB         2 0.012903226 10.242268 10.244198
HSPB1_HUMAN                     HSPB1         2 0.006451613 10.492542 10.495906
ACTN3_HUMAN                     ACTN3         2 0.006451613 10.381636 10.385495
                          t P.Value adj.P.Val        B
DPP2_HUMAN         4919.450       0         0 851.9517
Beta-galactosidase 6649.977       0         0 847.2963
PSME2_HUMAN        6749.099       0         0 841.5052
IGLC7_HUMAN        4064.202       0         0 834.0719
GPA33_HUMAN        4044.728       0         0 830.7692
MYO1D_HUMAN        4826.051       0         0 826.0299
PRDX1_HUMAN        5213.601       0         0 816.2985
KCRB_HUMAN         4266.358       0         0 807.7300
HSPB1_HUMAN        4006.021       0         0 806.1518
ACTN3_HUMAN        3763.088       0         0 804.7304
```

``` r
topTable(fit, coef = 6)
```

``` output
                                 Protein.Group      Protein.Names    Genes
DIGEST_Label13C15N                    P9999999 DIGEST_Label13C15N DL13C15N
Trypsin              contam_sp|P00761|TRYP_PIG            Trypsin  Trypsin
CTRC_HUMAN                              Q99895         CTRC_HUMAN     CTRC
MGA_HUMAN                               O43451          MGA_HUMAN     MGAM
CEL3A_HUMAN                             P09093        CEL3A_HUMAN   CELA3A
CDHR2_HUMAN                             Q9BYE9        CDHR2_HUMAN    CDHR2
CBPA1_HUMAN                             P15085        CBPA1_HUMAN     CPA1
DPEP1_HUMAN                             P16444        DPEP1_HUMAN    DPEP1
Serum              contam_sp|P02769|ALBU_BOVIN              Serum    Serum
MUC2_HUMAN                              Q02817         MUC2_HUMAN     MUC2
                   NPeptides    PropObs      logFC  AveExpr         t
DIGEST_Label13C15N        14 0.02672811  0.4528090 12.80679  5.756188
Trypsin                   26 0.05707196 -0.5672180 13.66182 -4.054078
CTRC_HUMAN                37 0.07061901 -0.7645132 14.23116 -3.872798
MGA_HUMAN                 92 0.06367461 -0.7376365 12.84200 -3.460817
CEL3A_HUMAN               34 0.06907021 -0.5811737 14.01516 -3.399523
CDHR2_HUMAN               40 0.08016129 -0.5862240 13.28108 -3.170783
CBPA1_HUMAN               49 0.05543120 -0.6468480 12.69218 -2.109360
DPEP1_HUMAN               27 0.06905615 -0.2534077 12.33003 -1.925096
Serum                     53 0.03992696  0.3432754 10.74604  1.803974
MUC2_HUMAN                79 0.04752960  0.4196754 11.82464  1.599716
                        P.Value    adj.P.Val          B
DIGEST_Label13C15N 4.693993e-08 1.905761e-05  8.2068985
Trypsin            8.058659e-05 1.635908e-02  1.2767389
CTRC_HUMAN         1.602282e-04 2.168421e-02  0.6739633
MGA_HUMAN          7.016031e-04 7.023562e-02 -0.5608611
CEL3A_HUMAN        8.649707e-04 7.023562e-02 -0.8696406
CDHR2_HUMAN        1.843621e-03 1.247517e-01 -1.4593639
CBPA1_HUMAN        3.657416e-02 9.999636e-01 -4.1158372
DPEP1_HUMAN        5.610980e-02 9.999636e-01 -4.6781628
Serum              7.324194e-02 9.999636e-01 -4.7376792
MUC2_HUMAN         1.117658e-01 9.999636e-01 -4.7700018
```

``` r
# Example diagnostic plots
plotMD(fit, coef = 1)
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-20-1.png" style="display: block; margin: auto;" />

``` r
plotMD(fit, coef = 2)
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-20-2.png" style="display: block; margin: auto;" />

## Example: aUC vs Ctrl Comparison


``` r
results <- topTable(fit, coef = 3, number = Inf)

# Visualize a specific protein
plotProtein(y.protein, "S10A9_HUMAN", col = as.character(Class.color))
legend('topleft', legend = levels(Class), fill = levels(Class.color))
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-21-1.png" style="display: block; margin: auto;" />

## Visualize Top Significant Proteins

We identify the top 50 most significant and variable proteins, then visualize via a clustered heatmap.


``` r
sig_proteins <- results %>%
  filter(adj.P.Val < 0.05) %>%
  top_n(50, wt = abs(logFC)) %>%
  pull(Protein.Names)

expr_matrix <- y.protein$E[sig_proteins, ]
scaled_expr <- t(scale(t(expr_matrix)))

col_annotation <- data.frame(Class = Class,
                             row.names = colnames(y.protein$E))

pheatmap(scaled_expr,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = col_annotation)
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-22-1.png" style="display: block; margin: auto;" />


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

``` warning
Warning: One or more p-values is 0. Converting to 10^-1 * current lowest
non-zero p-value...
```

``` warning
Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
ℹ Please use `linewidth` instead.
ℹ The deprecated feature was likely used in the EnhancedVolcano package.
  Please report the issue to the authors.
This warning is displayed once every 8 hours.
Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
generated.
```

``` warning
Warning: The `size` argument of `element_line()` is deprecated as of ggplot2 3.4.0.
ℹ Please use the `linewidth` argument instead.
ℹ The deprecated feature was likely used in the EnhancedVolcano package.
  Please report the issue to the authors.
This warning is displayed once every 8 hours.
Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
generated.
```

``` warning
Warning: Removed 2 rows containing missing values or values outside the scale range
(`geom_vline()`).
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-23-1.png" style="display: block; margin: auto;" />


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

# Convert UniProt IDs to Entrez IDs for enrichment
converted <- bitr(de_proteins$Protein.Group,
                  fromType = "UNIPROT",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
```

``` output
'select()' returned 1:many mapping between keys and columns
```

``` warning
Warning in bitr(de_proteins$Protein.Group, fromType = "UNIPROT", toType =
"ENTREZID", : 16.5% of input gene IDs are fail to map...
```

``` r
# GO Biological Process enrichment
ego <- enrichGO(gene = converted$ENTREZID,
                OrgDb = org.Hs.eg.db,
                ont = "BP",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)

dotplot(ego, showCategory = 10)
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-25-1.png" style="display: block; margin: auto;" />

``` r
# KEGG pathway enrichment
ekegg <- enrichKEGG(gene = converted$ENTREZID,
                    organism = 'hsa',
                    pvalueCutoff = 0.05)
```

``` output
Reading KEGG annotation online: "https://rest.kegg.jp/link/hsa/pathway"...
```

``` output
Reading KEGG annotation online: "https://rest.kegg.jp/list/pathway/hsa"...
```

``` r
dotplot(ekegg, showCategory = 10)
```

<img src="fig/datacleaning-2-rendered-unnamed-chunk-25-2.png" style="display: block; margin: auto;" />


## Potential next steps

* Compare coefficients between classes as an exercise.

* Visualize selected proteins of interest across experimental groups.

* Extend enrichment analyses using Reactome or GSEA approaches.

* Integrate sample metadata for clinical covariates or longitudinal modeling.




::::::::::::::::::::::::::::::::::::: keypoints 

- Data!

::::::::::::::::::::::::::::::::::::::::::::::::

