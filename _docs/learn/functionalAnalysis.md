---
layout: default
title: Functional Analysis
parent: Learn
nav_order: 7
---

## Code Description: 
In the Functional Analysis page, there were three major components : gene set enrichment analysis, KEGG analysis, and gene ontology. The first step is to generate a DEG vector with a given threshold value. Additionally the gene symbols in the data had to be converted to ENTREZIDs using the select() function. (ENTREZIDs are an alternative way to identify genes with databases.) Once the DEG vector is made, the GSEA() function is used for the gene set enrichment analysis. Then to visualize the output of the GSEA() function the gseaplot2() function is optimized. The next part of the functional analysis is the enriched KEGG analysis. To get the KEGG pathways, the enrichKEGG( ) function is used and to visualize the KEGG pathways, the dotplot() function and the barplot() function are optimized. The number of categories that were shown for the KEGG analysis was determined by the user. Then the setReadable() function was used so that the gene concept network for the significant KEGG pathways could be presented as well. Additionally, for gene ontology the user can choose which category of barplot they want to visualize for gene ontology. There are three different aspects the viewer can view: Biological Processes, Cellular Components, and Molecular Functions. Depending on which aspect the user decides to view the gene ontology, the enrichGO() function is used. The results were boxplots with the colors representing the p-value. 

## Significance of Output: 
The plots have significant importance as they represent the top DEGs from the data and with the ENTREZIDs there is additional information to more biological concepts. 

## Important Functions: 
In the functional analysis the select() function is important as it allows the data to be converted to the Entrenzids. The select() function allows the user’s data to identify genes from a NIH database. 
