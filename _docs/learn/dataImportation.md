---
layout: default
title: Data Importing
parent: Learn
nav_order: 1
---

# Code Description:
To import the different data files, the fileInput() function is used to allow multiple CSV files to be uploaded. For inputting the raw CEL files, the code includes a file input for tar zipped files that is then unzipped in the server function. The user can specify whether they have a microarray that requires the oligo package using radio buttons, and depending on their input, the code will use the correct package to read in the CEL files. The app can also get data from the GEO database using functions from the GEOQuery package. The user can input an accession number and can specify whether they want the series matrix or raw CEL files. Based on the user input, the code uses different GEOQuery functions to read in these files. Lastly, the code inputs metadata files in a similar manner to the CSV files of normalized and pre-processed expression data.

# Significance of Output: 
The data importation code will serve as the source for all the data importation functions in the final app. It will allow the app to be robust in being able to handle many file types.

# Advanced Functions: 
The GEOQuery package is important for Group B1 to learn how to use, particularly the getGEO() function and the getGEOSuppFiles() function.

