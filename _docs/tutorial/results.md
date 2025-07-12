---
layout: default
title: Output Structure
parent: Tutorial
nav_order: 5
---

## Output Structure

GDAP generates results in the following structure:

```
results/
└── {disease_name}/
    ├── embedding_wheel/
    │   ├── ggvec_wheel_model.bin
    │   ├── node2vec_wheel_model.bin
    │   └── prone_wheel_model.bin
    ├── classifier_models/
    │   ├── Logistic_Regression.pkl
    │   ├── Random_Forest.pkl
    │   └── model_evaluation_plots.png
    └── predictions/
        ├── associated_proteins.csv
        └── non_associated_proteins.csv
```
