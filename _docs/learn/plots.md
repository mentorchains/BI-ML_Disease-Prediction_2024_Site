---
layout: default
title: Visualizations
parent: Learn
nav_order: 8
---

## PCA
In the QC visualization, the user can choose to visualize either the boxplot or the PCA plot. The PCA plot was constructed using the prcomp() function and the ggplot() function. The prcomp() function was crucial for calculating the pca for all the genes. Then to visualize the PCA data, the ggplot() function was utilized. In our app, the user can choose the two features that the user wanted to analyze more. 

## Volcano Plot
The Volcano Plot is utilized in the Statistical Analysis for visualizing the top DEGs. The user specifies the LogFC and adjusted P-Values for the Volcano plot. The function used for generating the Volcano Plot is EnhancedVolcano which takes in the user input and the final result for the top table. The Volcano Plot allows the user to visualize how various processes have changed their data. 

## Boxplot
The Box plot was used to visualize the quality control in the aspect of normalization. To generate a BoxPlot the boxplot() function is used to generate the plot with the x axis representing the sample number and the y axis representing the gene expression values. The boxplot visualizes how the data changed from being normalized. 
