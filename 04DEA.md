---
title: 'Differential expression analysis (DEA)'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- How can we identify the differentially expressed proteins between experimental groups?

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- Conduct differential expression analysis using *limpa*
- Visualise the top differentially expressed proteins

::::::::::::::::::::::::::::::::::::::::::::::::

## Differential Expression Analysis

Now it is time to analyse our data and investigate proteins that are differentially expressed between our two experimental groups (control vs active Crohn’s disease). We will make use of the *limma* package and closely related *limpa* functions specifically designed for proteomics data. 

First, we will fit a linear model that includes both Class and Batch effects.


``` r
# Set the control group as our reference
Class <- relevel(Class, ref = "Ctrl")
```

``` error
Error: object 'Class' not found
```

``` r
## We also need to reset our colors to match the correct labels
Class.color <- Class
```

``` error
Error: object 'Class' not found
```

``` r
levels(Class.color) <- hcl.colors(nlevels(Class), palette = "cividis")
```

``` error
Error: object 'Class' not found
```

``` r
# Fit a linear model
design <- model.matrix(~Class + Batch)
```

``` error
Error in eval(predvars, data, env): object 'Class' not found
```

``` r
design
```

``` error
Error: object 'design' not found
```

We can see both Class and Batch accounted for in our design matrix. 

Now we run the *limpa* differential expression analysis function `dpcDE()` to identify the top differentially expressed proteins in our dataset.


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
topTable(fit, coef = 2)
```

``` error
Error in topTable(fit, coef = 2): could not find function "topTable"
```

We can visualise the log2 expression values for differentially expressed proteins, together with standard errors.


``` r
plotProtein(y.protein, "PRTN3_HUMAN", col = as.character(Class.color))
```

``` error
Error in plotProtein(y.protein, "PRTN3_HUMAN", col = as.character(Class.color)): could not find function "plotProtein"
```

``` r
legend('topleft', legend = levels(Class), fill = levels(Class.color))
```

``` error
Error: object 'Class' not found
```

``` r
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

## Visualise differentially expressed proteins

### Heatmap

We can filter our results to the top up to 50 most significant differentially expressed proteins, then visualise via a clustered heatmap.


``` r
# Save our DEA results as a dataframe
results <- topTable(fit, coef = 2, number = Inf)
```

``` error
Error in topTable(fit, coef = 2, number = Inf): could not find function "topTable"
```

``` r
# Filter for the up to 50 most significant results
sig_proteins <- results %>%
  filter(adj.P.Val < 0.05) %>% # filter by p value
  top_n(50, wt = abs(logFC)) %>% # filter by absolute log fold change
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
# Scale the data and visualise via heatmap
scaled_expr <- t(scale(t(expr_matrix)))
```

``` error
Error: object 'expr_matrix' not found
```

``` r
pheatmap(scaled_expr,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         annotation_col = y.protein$targets,
         clustering_method = 'ward.D')
```

``` error
Error in pheatmap(scaled_expr, cluster_rows = TRUE, cluster_cols = TRUE, : could not find function "pheatmap"
```

:::challenge

How would you interpret this heatmap?

:::solution

As expected, our samples cluster into their respective classes. 

The three proteins with decreased expression in the aCD samples are clearly distinguishable from those with increased expression in the aCD samples.

There are no apparent batch effects confounding these results.

:::
:::::

### Volcano plot

We can also visualise our results as a volcano plot using `EnhancedVolcano()`.


``` r
EnhancedVolcano(results,
                lab = results$Genes,
                x = 'logFC',
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2.0,
                labSize = 3.5,
                drawConnectors = TRUE,
                maxoverlapsConnectors = Inf,
                lengthConnectors = unit(0.001, 'npc'),
                typeConnectors ='closed',
                ylim = c(0, 5),
                xlim = c(-5, 5))
```

``` error
Error in EnhancedVolcano(results, lab = results$Genes, x = "logFC", y = "adj.P.Val", : could not find function "EnhancedVolcano"
```


::::::::::::::::::::::::::::::::::::: keypoints 

- We can fit a linear model to our data and conduct differential expression analysis using *limpa* and *limma* functions.
- There are many different ways to visualise differentially expressed proteins.

::::::::::::::::::::::::::::::::::::::::::::::::

