---
layout: default
title: Quality Control
parent: Learn
nav_order: 4
---

## QC

### Code Description: 
In the beginning for the User Interface for the Quality Control page, the user chooses the type of QC visualization before the normalization process. The two choices the user gets for batch correction are NUSE and RLE. If the user chooses the NUSE method for quality control then at the server end, depending on the file that was imported, the quality control is completed using the NUSE() function. If the user chooses the RLE method, then in the server end the RLE() function is used for quality control using different methods depending on the file type. Once that is completed, the user can see the visualization of the data before the normalization and batch correction. The user can choose to visualize the data after the normalization with either a Boxplot or a PCA plot. On the server side, if the user chooses Boxplot then the boxplot() function is used to generate the plot with the x axis representing the sample number and the y axis representing the gene expression values. However if the user chooses to visualize the data using a PCA plot the prcomp() function is used to calculate the pca values. Then the user will choose the components for plotting the PCA plot. The user will receive input from the server side if the user enters more or less than 2 components and if the user did not specify a feature. Then on the server side, the ggplot() function outputs the PCA plot.

### Significance of Output: 
 Quality control is important because it helps check if the data is reliable. The QC visualizations before and after normalization helps see the quality control reports, but it also shows how normalization can change the data.

### Advanced Functions: 
The ggplot() is very useful when not only plotting the PCA plot but also for various other types of plots and visualizations. Importation functions for QC were RLE() and NUSE(). The RLE and NUSE functions can produce a boxplot for arrays from a dataframe. 
