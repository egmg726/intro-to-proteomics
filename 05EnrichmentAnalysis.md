---
title: 'Network and Enrichment Analysis'
teaching: 10
exercises: 2
---

:::::::::::::::::::::::::::::::::::::: questions 

- What is DE analysis

::::::::::::::::::::::::::::::::::::::::::::::::

::::::::::::::::::::::::::::::::::::: objectives

- To do DE analysis

::::::::::::::::::::::::::::::::::::::::::::::::

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

- DEA

::::::::::::::::::::::::::::::::::::::::::::::::

