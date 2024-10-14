---
layout: default
title: Normalization (RMA GCRMA and MAS 5)
parent: Learn
nav_order: 2
---


### Code Description: 
The user has three options for normalizing their data: RMA, GCRMA, and MAS5. In the user interface , the user will specify between either of these methods using radio buttons and select the button Begin Normalization in order for the data to begin normalizing. The user is warned that normalization only has to be done for CEL files and RMA normalization is the only acceptable type of normalization for Gene ST Arrays and Exon ST Arrays. Based on the method the user chooses, the normalization takes place using the rma() function if the user chooses the RMA method, gcrma() function if the user chooses the GCRMA method, and mas5() function if the user chooses the MAS5 method.

### Significance of Output: 
Normalization of data is very important because it will make feature extraction easier and it will make the organization of the data much more structured. 

### Advanced Functions: 
The rma() function was used if the user chose the rma method. The RMA() function was used to take in an AffyBatch object and turn it into an expression set. The gcrma() function is used if the user chooses the gcrma method. The gcrma() function takes in the affybatch object and turns it into an expression set using the probe sequence. If the user chooses the mas5 method the MAS5() function would be used. The MAS5 function will normalize each array in the data independently while rma and gcrma use a multi-chip model.
