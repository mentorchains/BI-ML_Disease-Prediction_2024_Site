---
layout: default
title: Statistical Analysis
parent: Learn
nav_order: 6
---

## Code Description: 
Under the Statistical Analysis tab the user is asked what percent of genes they want filtered.The genes below the indicated percentile are removed from the data. However the annotated data is getting filtered. Prior to the genes being filtered, the data went through annotation in which the duplicates, NA, and duplicate gene symbols are removed. After the user sees that the genes were filtered the user will have the option to group the data based on the metadata feature for limma analysis. Then the user is asked to select which metadata features they want to compare specifically. The code then generates a data frame with the DEGs (differentially expressed genes) using the topTable() function. On the server side the first step was to determine the type of data that was being used, the filtered data or the non filtered data... Once that is complete, it is important to make sure that the expression data only includes the data that the user has selected. Then it is time to create the top table by first creating a design matrix, and then going through limma analysis. Using model.matrix() function a matrix is made, which will then be used for the limma analysis. For the limma analysis, the lmfit() function and the eBays() function are used. Limma analysis is important because it helps analyze gene expression data. For limma analysis, the top table is created with the topTable() function, with the user inputting the cutoff p-value and the percent of DEGs the user wants to view. After the top table is generated, a volcano plot is used to visualize the DEGs. The user gets to distinguish the LogFC and the adjusted p-value cutoff for the volcano plot. 

## Significance of Output: 
The final input is the volcano plot which represents a lot of important components. An important part is the top DEGs which are filtered, annotated, and have undergone the limma analysis. The DEG analysis also represents how the data has changed through the processes of limma analysis and filtering.

## Important Functions: 
The dataTableProxy() function is used to manipulate the existing data table created to select or add rows/columns. This is an important function because it allows the server code to be flexible depending on what the user inputs. The lmfit() function is also crucial for the limma analysis which allowed the data to fit a linear model for the genes and allowed the analysis of the data. The toptable() function is also important because it calculates the DEGs and puts them into a table. The EnhancedVolcano() function is crucial for creating the volcano plot that is used to visualize the top table with the DEGs. 
