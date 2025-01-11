---
layout: default
title: Navigating through GDAP
parent: Overview
nav_order: 4

---

## Navigating through GDAP (Gene-Disease Association Prediction)
{: .no_toc }

<details open markdown="block">
  <summary>
    Table of contents
  </summary>
  {: .text-delta }
1. TOC
{:toc}
</details>

---

The GDAP project follows a modular structure where each component handles a specific part of the [gene-disease association prediction pipeline](https://github.com/mentorchains/BI-ML_Disease-Prediction_2024/blob/ahmed/src/main.py). Here's how to navigate through the main components:

### Data Input

  1. Initialize the pipeline using the [`Config`](https://github.com/mentorchains/BI-ML_Disease-Prediction_2024/blob/ahmed/src/config.py#L11-L21) class.
    
  2. Set parameters such as:
    - `DATA_SOURCE`: Choose between "`GraphQLClient`" or "`BigQueryClient`"
    - `DISEASE_ID`: Specify the disease ID (EFO format)
    - `PPI_INTERACTIONS`: Maximum number of PPI interactions to consider
    - `NEGATIVE_TO_POSITIVE_RATIO`: Ratio for negative to positive edges

  3. Load the data using [`load_data()`](https://github.com/mentorchains/BI-ML_Disease-Prediction_2024/blob/ahmed/src/main.py#L19C5-L42C29) function:

  ```python
  config = Config()
  pipeline = GeneDisease(config)
  ppi_df, ot_df = pipeline.load_data()
  ```
---

### Graph Construction

After loading the data, the graph construction step integrates disease-gene associations and PPI data into a bipartite graph using [`create_graph()`](https://github.com/mentorchains/BI-ML_Disease-Prediction_2024/blob/ahmed/src/main.py#L44-L50) function:

  ```python
  pipeline.create_graph(ppi_df, ot_df)
  ```

  - This creates:
    - `pipeline.G`: The main graph object
    - `pipeline.pos_edges`: Positive edge connections (disease-gene associations)
    - `pipeline.neg_edges`: Negative edge connections (PPI network)

---

### Embedding Generation

The embedding generation is a crucial step in GDAP that significantly impacts the quality of predictions. Using [`generate_embeddings()`](https://github.com/mentorchains/BI-ML_Disease-Prediction_2024/blob/ahmed/src/main.py#L52-L75), GDAP supports multiple embedding methods to be configered in the [`Config`](https://github.com/mentorchains/BI-ML_Disease-Prediction_2024/blob/186ac67b9e99b2e2b438d7b7f381ecae5e0da432/src/config.py#L17) class:

1. Simple node embedding based on degree
  - Default Parameters:
    - Dimension: 64


2. Node2Vec (Node2Vec)
  - Default Parameters:
    - Walk length: 10
    - Components: 32
    - Output: Embedding model `n2v_model` and vector embeddings `n2v_wheel_model.bin`


3. ProNE (ProNE)
  - Default Parameters:
    - Components: 32
    - Parameters: step=5, mu=0.2, theta=0.5
    - Output: Embedding model `prone_model` and vector embeddings `prone_wheel_model.bin`


4. GGVec (GGVec)
  - Default Parameters:
    - Components: 64
    - Order: 3
    - Output: Embedding model `ggvec_model` and vector embeddings `ggvec_wheel_model.bin`


The quality of these embeddings directly affects GDAP's ability to capture the complex relationships between genes and diseases, making this step vital for accurate predictions.

---

### Feature Extraction

The feature extraction process converts node embeddings into edge features using [`extract_features_labels()`](https://github.com/mentorchains/BI-ML_Disease-Prediction_2024/blob/ahmed/src/main.py#L77-L90):

  ```python
  X_train, y_train, X_val, y_val, X_test, y_test, edges_train, edges_val, edges_test = pipeline.extract_features_labels(config.TEST_SIZE)
  ```

This step:
  - Maps nodes to indices for embedding lookup
  - Generates feature matrices and labels
  - Splits data into training, validation, and test sets

---

### Model Selection

The model selection and training process is handled by [`train_and_evaluate_model()`](https://github.com/mentorchains/BI-ML_Disease-Prediction_2024/blob/ahmed/src/main.py#L92-L99):

  ```python
  model, models_path = pipeline.train_and_evaluate_model(X_train, y_train, X_test, y_test, X_val, y_val)
  ```

This step:

  1. Trains the selected classification model
  2. Evaluates performance using multiple metrics:
    - Accuracy
    - Precision
    - Recall
    - F1 Score

  3. Saves the trained model to disk

Model evaluation plots can be generated using:

  ```python
  pipeline.plot_model_evaluation(model, X_val, y_val, models_path)
  ```

---

### Model Prediction

Finally, predictions are made using [`predict_and_save_results()`](https://github.com/mentorchains/BI-ML_Disease-Prediction_2024/blob/ahmed/src/main.py#L122-L127):

  ```python
  associated_df, non_associated_df = pipeline.predict_and_save_results(model, X_val, edges_val, models_path)
  ```


This generates:

  - associated_df: Predicted positive associations
  - non_associated_df: Predicted negative associations

The results are saved in the following directory structure:

  ```sh
  output_dir/
  └── disease_name/
      ├── embedding_wheel/
      │   ├── n2v_model
      │   ├── n2v_wheel_model.bin
      │   └── ...
      ├── network/
      │   ├── graph.graphml
      │   ├── negative_edges.npy
      │   ├── positive_edges.npy
      │   └── ...
      ├── predictions/
      │   ├── associated_prediction_results.csv
      │   └── non_associated_results.csv
      └── classifier_model/
          └── model_name.pkl
  ```
