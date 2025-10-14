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
```

``` error
Error in readDIANN(file = "data/report_subset.parquet", format = "parquet"): arrow package required but is not installed (or can't be loaded)
```

``` r
names(y.peptide)
```

``` error
Error: object 'y.peptide' not found
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

``` error
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function '%in%': object 'y.peptide' not found
```

``` r
# Filter out contaminants
y.peptide$genes <- dplyr::filter(y.peptide$genes, !(Protein.Group %in% contaminants$`Uniprot ID`))
```

``` error
Error: object 'y.peptide' not found
```

``` r
y.peptide$E <- y.peptide$E[rownames(y.peptide$genes), ]
```

``` error
Error: object 'y.peptide' not found
```

``` r
# Confirm removal
table(y.peptide$genes$Protein.Group %in% contaminants$`Uniprot ID`)
```

``` error
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function '%in%': object 'y.peptide' not found
```


### Align sample metadata



``` r
sample_info <- samples_stool[samples_stool$sample_name %in% colnames(y.peptide$E),]
```

``` error
Error in h(simpleError(msg, call)): error in evaluating the argument 'table' in selecting a method for function '%in%': error in evaluating the argument 'x' in selecting a method for function 'colnames': object 'y.peptide' not found
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
# Add experimental metadata
y.peptide$targets <- sample_info[, c('Class', 'Batch')]
```

``` error
Error: object 'sample_info' not found
```

``` r
# Apply peptide-level filters
y.peptide <- filterNonProteotypicPeptides(y.peptide)
```

``` error
Error: object 'y.peptide' not found
```

``` r
y.peptide <- filterCompoundProteins(y.peptide)
```

``` error
Error: object 'y.peptide' not found
```

``` r
y.peptide <- filterSingletonPeptides(y.peptide, min.n.peptides = 2)
```

``` error
Error: object 'y.peptide' not found
```


For more information on the peptide-level filters, please consult the limpa vignette.


### Dimensionality Reduction and QC


``` r
dpcfit <- dpc(y.peptide)
```

``` error
Error: object 'y.peptide' not found
```

``` r
plotDPC(dpcfit)
```

``` error
Error: object 'dpcfit' not found
```


By protein:


``` r
y.protein <- dpcQuant(y.peptide, "Protein.Names", dpc=dpcfit)
```

``` error
Error: object 'y.peptide' not found
```

``` r
plotMDSUsingSEs(y.protein)
```

``` error
Error: object 'y.protein' not found
```

That plot is extremely messy and doesn't tell us much! Let's clean it up.



``` r
# Class visualisation
Class <- factor(y.peptide$targets$Class)
```

``` error
Error: object 'y.peptide' not found
```

``` r
levels(Class) <- c("Ctrl","aCD","aUC","CDr","UCr")
```

``` error
Error: object 'Class' not found
```

``` r
Class.color <- Class
```

``` error
Error: object 'Class' not found
```

``` r
levels(Class.color) <- c("pink","black",'grey','darkorange','darkgreen')
```

``` error
Error: object 'Class.color' not found
```

``` r
plotMDSUsingSEs(y.protein, pch=16, col=as.character(Class.color))
```

``` error
Error: object 'y.protein' not found
```

Hmmm, not a lot of differences. Let's see if there are any batch effects.


``` r
# Batch visualisation
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
levels(Batch.color) <- c("green","purple",'orange')
```

``` error
Error: object 'Batch.color' not found
```

``` r
plotMDSUsingSEs(y.protein, pch=16, col=as.character(Batch.color))
```

``` error
Error: object 'y.protein' not found
```

It looks like there are some differences between the batches. Let's remove the batch effects and see how it looks.



``` r
# Correct for batch effects
y.protein.rbe <- removeBatchEffect(y.protein,batch = y.protein$targets$Batch)
```

``` error
Error: object 'y.protein' not found
```

``` r
plotMDS(y.protein.rbe, pch=16, col=as.character(Batch.color))
```

``` error
Error: object 'y.protein.rbe' not found
```

## Subset Comparison: Ctrl vs aUC


``` r
cdr_bool <- y.peptide$targets$Class %in% c('Ctrl','aUC')
```

``` error
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function '%in%': object 'y.peptide' not found
```

``` r
plotMDSUsingSEs(y.protein[, cdr_bool], 
                pch=16, col=as.character(Class.color[cdr_bool]))
```

``` error
Error: object 'y.protein' not found
```

``` r
plotMDS(y.protein[, cdr_bool], 
         pch=16, col=as.character(Class.color[cdr_bool]))
```

``` error
Error: object 'y.protein' not found
```

``` r
plotMDS(y.protein.rbe[, cdr_bool], 
         pch=16, col=as.character(Class.color[cdr_bool]))
```

``` error
Error: object 'y.protein.rbe' not found
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
Error: object 'y.protein' not found
```

``` r
fit <- eBayes(fit)
```

``` error
Error in eBayes(fit): fit is not a valid MArrayLM object
```


``` r
topTable(fit, coef = "ClassCtrl")
```

``` error
Error in topTable(fit, coef = "ClassCtrl"): fit must be an MArrayLM object
```

## Create and Explore Contrasts

We can define custom contrasts (e.g., comparing aCD vs Ctrl) and explore differential proteins.


``` r
# Define contrasts
contrasts_fit <- makeContrasts(Ctrl_aCD = ClassaCD - ClassCtrl, 
                               levels = colnames(design))
```

``` error
Error in if (levels[1] == "(Intercept)") {: argument is of length zero
```

``` r
# Apply contrasts
fit2 <- contrasts.fit(fit, contrasts = contrasts_fit)
```

``` error
Error: object 'contrasts_fit' not found
```

``` r
fit2 <- eBayes(fit2)
```

``` error
Error: object 'fit2' not found
```

``` r
# Summaries
topTable(fit2, coef = 1)
```

``` error
Error: object 'fit2' not found
```


``` r
# Compare multiple coefficients
topTable(fit, coef = 1)
```

``` error
Error in topTable(fit, coef = 1): fit must be an MArrayLM object
```

``` r
topTable(fit, coef = 3)
```

``` error
Error in topTable(fit, coef = 3): fit must be an MArrayLM object
```

``` r
topTable(fit, coef = 4)
```

``` error
Error in topTable(fit, coef = 4): fit must be an MArrayLM object
```

``` r
topTable(fit, coef = 6)
```

``` error
Error in topTable(fit, coef = 6): fit must be an MArrayLM object
```

``` r
# Example diagnostic plots
plotMD(fit, coef = 1)
```

``` error
Error in as.vector(x, mode): cannot coerce type 'closure' to vector of type 'any'
```

``` r
plotMD(fit, coef = 2)
```

``` error
Error in as.vector(x, mode): cannot coerce type 'closure' to vector of type 'any'
```

## Example: aUC vs Ctrl Comparison


``` r
results <- topTable(fit, coef = 3, number = Inf)
```

``` error
Error in topTable(fit, coef = 3, number = Inf): fit must be an MArrayLM object
```

``` r
# Visualize a specific protein
plotProtein(y.protein, "S10A9_HUMAN", col = as.character(Class.color))
```

``` error
Error: object 'y.protein' not found
```

``` r
legend('topleft', legend = levels(Class), fill = levels(Class.color))
```

``` error
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'levels': object 'Class' not found
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
Error: object 'results' not found
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
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 't': error in evaluating the argument 'x' in selecting a method for function 't': object 'expr_matrix' not found
```

``` r
col_annotation <- data.frame(Class = Class,
                             row.names = colnames(y.protein$E))
```

``` error
Error in h(simpleError(msg, call)): error in evaluating the argument 'x' in selecting a method for function 'colnames': object 'y.protein' not found
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
Error: object 'scaled_expr' not found
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
Error: object 'results' not found
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
Error: object 'results' not found
```

``` r
# Convert UniProt IDs to Entrez IDs for enrichment
converted <- bitr(de_proteins$Protein.Group,
                  fromType = "UNIPROT",
                  toType = "ENTREZID",
                  OrgDb = org.Hs.eg.db)
```

``` error
Error: object 'de_proteins' not found
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
Error: object 'converted' not found
```

``` r
dotplot(ego, showCategory = 10)
```

``` error
Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'dotplot': object 'ego' not found
```

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

``` error
Error: object 'converted' not found
```

``` r
dotplot(ekegg, showCategory = 10)
```

``` error
Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'dotplot': object 'ekegg' not found
```


## Potential next steps

* Compare coefficients between classes as an exercise.

* Visualize selected proteins of interest across experimental groups.

* Extend enrichment analyses using Reactome or GSEA approaches.

* Integrate sample metadata for clinical covariates or longitudinal modeling.




::::::::::::::::::::::::::::::::::::: keypoints 

- Data!

::::::::::::::::::::::::::::::::::::::::::::::::

