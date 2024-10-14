---
layout: default
title: Outliers (Removing Potential Outliers)
parent: Learn
nav_order: 5
---

## Code Description: 
A crucial part for making sure the data the user enters is reliable is to make sure the outliers are identified and removed. The user can choose from three different methods for identifying the outliers; the three methods are KS, sum, and upper quartile. When the user presses the button Find Potential Outliers after selecting a method, the outliers() function is used to find the outliers. Then the user can visualize the potential outliers with the ggplot() function. Then the user is asked what outliers should be removed depending on the result they see. .. After the outliers are removed, the expression matrix is updated on the server end. Then once this process is complete, the user can select the button Show Updated List of Samples and the user will be able to see the data without the outliers. 

## Significance of Output: 
The process of identifying and removing outliers is really important because it allows the data to be more reliable. It also lets the user see how the data changed after the outliers are removed. The user can be sure that extreme samples are not causing the results they see, and that results can truly be attributed to differences in gene function and expression.

## Advanced Functions: 
An important function in finding and removing functions is the outliers(). The outliers() performs outlier detection using three different methods.
