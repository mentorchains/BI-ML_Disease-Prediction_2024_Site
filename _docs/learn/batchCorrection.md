---
layout: default
title: Batch Correction
parent: Learn
nav_order: 3
---

### Code Description: 
Before batch correction occurs there is a drop down menu in the User Interface that asks the user to select which metadata column should provide the batch assignment for each of the samples. Then the user selects the button Perform Batch Correction and the process of batch correction takes place in the Server using the ComBat() function. Once the server code runs and batch correction is complete, the user will see a text output saying “Batch Correction is Complete”. 

### Significance of Output: 
Batch Correction is useful when using various numbers of batches. Batch Correction makes two batches comparable to one another by correcting for differences in equipment and experimentation. As a result, after batch correction, the two batches can be used for the same analysis.

### Advanced Functions: 
We used the ComBat() function to allow the user to adjust their data for batch correction. The imputed object is the normalized data and the returned value is an expression matrix that was corrected for batch effects. 
