---
layout: default
title: Overview
nav_order: 1
has_children: true
---

# Overview

GDAP (Gene-Disease Association Prediction) is a comprehensive machine learning pipeline designed to predict novel gene-disease associations using graph-based approaches and network analysis.

## Purpose

The primary goal of GDAP is to identify potential gene-disease relationships that may not be immediately apparent through traditional experimental methods. This is particularly valuable for:

- **Drug Discovery**: Identifying new therapeutic targets
- **Disease Research**: Understanding genetic mechanisms of diseases
- **Personalized Medicine**: Supporting precision medicine approaches

## How It Works

GDAP operates through a five-step pipeline:

1. **Data Collection**: Gathers disease-gene associations from Open Targets Platform and protein-protein interactions from STRING database
2. **Graph Construction**: Builds a bipartite graph connecting diseases to genes with weighted edges
3. **Embedding Generation**: Creates node representations using graph embedding algorithms
4. **Feature Extraction**: Extracts edge features for machine learning
5. **Model Training**: Trains binary classifiers to predict novel associations

## Key Components

- **Data Sources**: Open Targets Platform (disease-gene data) and STRING database (protein interactions)
- **Embedding Algorithms**: Node2Vec, ProNE, GGVec, and degree-based embeddings
- **Machine Learning Models**: Random Forest, SVM, Logistic Regression, Gradient Boosting
- **Web Interface**: Streamlit app for interactive exploration and visualization
