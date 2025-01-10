---
layout: default
title: Session Workflow
parent: Overview
nav_order: 2

---

## GDAP session workflow

<figure style="text-align: center;">
  <img src="/assets/gdap/gdap-workflow.png" alt="gdap-workflow" style="width: 80%;">
  <figcaption>GDAP Workflow</figcaption>
</figure>

The GDAP workflow predicts gene-disease associations by leveraging disease-gene association data and Protein-Protein Interaction (PPI) data.

The GDAP workflow comprises three primary steps:

1. [Data Collection](#1-data-collection): Fetching disease-gene and PPI data.
2. [Graph Construction](#2-graph-construction): Building a graph that integrates disease-gene and PPI data.
3. [Embedding Generation](#3-embedding-generation): Using embedding models to represent graph nodes for prediction tasks.
4. [Feature Extraction](#4-feature-extraction-and-data-splitting): Extract edge features, and split the data. 
5. [Model Selection](#5-model-selection): Training the model and generating predictions based on the features.

## **1. Data Collection**

The first step in the GDAP workflow involves gathering the necessary data from two critical datasets:

### 1.1 Open Targets Dataset

The [Open Targets platform][opentargets-platform] provides a comprehensive resource for identifying therapeutic drug targets, offering data on [disease-gene associations][opentargets-disease]. GDAP utilizes two primary methods for accessing this data: the GraphQL API and Google BigQuery. Both methods offer different approaches for querying disease-target association data which are identified by [Experimental Factor Ontology (EFO)][efo] terms.

<figure style="text-align: center;">
  <img src="/assets/gdap/overview/ot-sample-data-indirect-scores.png" alt="sample of ppi data" style="width: 50%;">
  <figcaption>Sample of disease-target data</figcaption>
</figure>


#### Data Sources:

**(i) GraphQL API:**

- The GraphQL API enables you to retrieve disease-target association data using a flexible [GraphQL querying][graphql-ex] structure. You can interact with the Open Targets [GraphQL schema][graphql-schema], which allows you to request detailed information about diseases and their associated targets, along with scores and other relevant metrics using the disease's experimental variable.

    <details>
        <summary><b>GraphQLClient Endpoints and Queries</b></summary>
        <ul>
            <li>
                <code>GraphQLClient</code> interacts with the Open Targets
                <a href="https://api.platform.opentargets.org/api/v4/graphql/schema" target="_blank">GraphQL Schema</a> to fetch disease-target association data.
            </li>

            <li>
                According to the Open Targets <a href="https://community.opentargets.org/t/returning-all-associations-data-using-the-platform-api/324/2" target="_blank">community forum</a>, when using the GraphQL API, you need to iterate
                through each page of the disease to retrieve the overall results, as their API provides 50 results for each query. This is also outlined in their
                <a href="https://youtu.be/_sZR0VxpwqE?si=kJ7pmAg_Uh3PLbU4&t=2384" target="_blank">webinar</a>.
            </li>

            <li>
                The <code>GraphQLClient</code> first fetches the total number of targets associated with a specific disease using its EFO ID. This provides the necessary data to calculate how many pages of results need to be retrieved. Then, it
                iterates through each page, retrieving detailed disease-target association data, including target symbols and scores, until all results are collected.
            </li>

            <li>
                In summary, the class performs the following steps:
                <ul>
                    <li>1. Retrieves the total count of targets associated with the disease.</li>
                    <li>2. Iterates through all pages of results, collecting the complete dataset.</li>
                </ul>
            </li>

            <li>
                You can also run GraphQL queries through the
                <a
                    href="https://api.platform.opentargets.org/api/v4/graphql/browser?query=%23%20GraphQL%20Query%20example%20used%20in%20GDAP%0Aquery%20GetDiseaseTargets%20%7B%0A%20%20%23%20Disease%20example%3A%20cardiovascular%20disease%0A%20%20disease%28efoId%3A%20%22EFO_0000319%22%29%20%7B%0A%20%20%20%20id%0A%20%20%20%20name%0A%20%20%20%20associatedTargets%28page%3A%20%7B%20size%3A%2050%2C%20index%3A%201%20%7D%29%20%7B%0A%20%20%20%20%20%20count%0A%20%20%20%20%20%20rows%20%7B%0A%20%20%20%20%20%20%20%20target%20%7B%0A%20%20%20%20%20%20%20%20%20%20id%0A%20%20%20%20%20%20%20%20%20%20approvedSymbol%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20score%0A%20%20%20%20%20%20%20%20datasourceScores%20%7B%0A%20%20%20%20%20%20%20%20%20%20id%0A%20%20%20%20%20%20%20%20%20%20score%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%7D%0A%20%20%7D%0A%7D%0A"
                    target="_blank"
                >
                    GraphQL API Playground
                </a>
                to explore or customize your queries.
            </li>
        </ul>
    </details>


**(ii) Google BigQuery:**

- Google BigQuery offers [asynchronous SQL querying][bigquery-docs-ex] capabilities through [Google Cloud][bigquery-console]'s infrastructure, providing an alternative method for extracting disease-target data.

    <details>
        <summary><b>BigQueryClient Endpoints and Queries</b></summary>
        <ul>
            <li>
                The <code>BigQueryClient</code> class allows you to query disease-target association data using
                <a href="https://console.cloud.google.com/marketplace/product/bigquery-public-data/open-targets-platform" target="_blank">Google BigQuery</a>. It supports two main SQL queries: <code>DIRECT_SCORES</code> and
                <code>INDIRECT_SCORES</code>, which can be executed against the <a href="https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/24.06/output/etl/json/" target="_blank">Open Targets BigQuery database</a>.
            </li>

            <ul>
                <li><code>DIRECT_SCORES</code>: Retrieves the direct associations between diseases and targets, along with relevant scores and evidence counts.</li>
                <li><code>INDIRECT_SCORES</code>: Retrieves the indirect associations between diseases and targets, with corresponding scores and evidence counts.</li>
            </ul>

            <li>
                You can also use these SQL queries directly in the BigQuery console for your specific disease and download the dataset in formats like CSV or JSON. Access the
                <a href="https://console.cloud.google.com/bigquery?sq=119305062578:ae33e769ecec43ceb055565f7d6c74df" target="_blank">console</a>.
            </li>
        </ul>
    </details>

        
### 1.2 String Database

The STRING database provides a comprehensive [network][string-ppi-network] of protein-protein interactions (PPIs) for Homo sapiens. Each interaction is assigned a [combined score][string-combined-score], a probabilistic measure derived from multiple evidence sources, such as experimental data and text mining. This score quantifies the likelihood of interaction, forming the backbone for integrating disease-gene data with protein interactions in the GDAP workflow.

<figure style="text-align: center;">
  <img src="/assets/gdap/overview/ppi-sample-data.png" alt="sample of ppi data" style="width: 50%;">
  <figcaption>Sample of PPI Data</figcaption>
</figure>

## **2. Graph Construction**

The graph construction step forms the core of the GDAP pipeline. It integrates both disease-gene associations and protein-protein interactions to form a bipartite graph.

This stage integrates the disease-gene associations and PPI data into a bipartite graph:

- **Nodes:** Represent the target (disease), and the genes (proteins).

- **Edges:** Represent relationships between nodes.
    
    - Positive edges represent disease-gene associations (edges weighted by their relevance scores from Open Targets).
    - Negative edges represent non-associated gene pairs, weighted by their interaction combined scores from the STRING database.

- The goal here is to construct a network that connects disease to genes and represents their biological relationships, along with possible interactions between genes.

<div style="display: flex; justify-content: space-between; align-items: center;">

  <figure style="text-align: center; width: 48%;">
    <img src="/assets/gdap/overview/bigraph-negative-edge.png" alt="Negative PPI edge" style="width: 100%;">
    <figcaption>Negative PPI edge</figcaption>
  </figure>

  <figure style="text-align: center; width: 48%;">
    <img src="/assets/gdap/overview/bigraph-positive-edge.png" alt="Positive disease-gene edge" style="width: 100%;">
    <figcaption>Positive disease-gene edge</figcaption>
  </figure>

</div>

## **3. Embedding Generation**

Node embeddings are generated using various [nodevectors][nodevectors] models to represent graph nodes in a lower-dimensional space.

These embeddings capture the structural and relational properties of nodes within the graph, ensuring that nodes representing similar biological entities or interactions are closely aligned in the vector space, which facilitates accurate predictions.

<figure style="text-align: center;">
  <img src="/assets/gdap/overview/node-embeddings.png" alt="sample of ppi data" style="width: 100%;">
  <figcaption>Node embeddings generation</figcaption>
</figure>


## **4. Feature Extraction and Data Splitting**

Edge features are generated by combining the embeddings of the two nodes connected by each edge, capturing the relationship between them.

- For each edge, the embeddings of the connected nodes (disease and gene) are retrieved.
- The edge feature is created by multiplying these embeddings, resulting in a feature vector representing the relationship.
- The data is then split into training, validation, and testing sets: the training set is used for model fitting, the validation set for hyperparameter tuning, and the testing set to evaluate the model's generalizability.

These edge features, derived from node embeddings, are essential for capturing relationships like disease-gene or protein-protein interactions. They are used to train a machine learning model for tasks such as predicting novel disease-gene associations or interactions.


## **5. Model Selection**

The final step involves training a binary classifier to predict disease-gene associations. Currently, GDAP offers traditional machine learning algorithms and a simple dense sequential model.

**(i) Model Training:**

- **Binary Classification:** The classifier is trained using edge features. Positive edges (true associations) are labeled as 1, while negative edges (random or non-associating gene pairs) are labeled as 0.
- During training, the classifier learns to distinguish between true disease-gene associations and random connections, based on features derived from node embeddings and interaction strengths.

**(ii) Model Prediction:**
- After training, the model can predict new disease-gene associations or protein interactions by applying the learned model to novel edge pairs. The output is a probability score representing the likelihood that an edge (disease-gene or protein-protein interaction) is real.
- A threshold is applied to this score to decide whether the predicted edge is a positive or negative association (default threshold = 0.5).

**(iii) Model Evaluation:**
- The model’s performance is evaluated using metrics such as accuracy, precision, recall, and F1 score. These metrics assess how well the model generalizes to unseen data, ensuring that the predictions are meaningful and reliable.


<!--

===========
REFERENCES
===========
-->

[opentargets-platform]: https://platform.opentargets.org/
[opentargets-disease]: https://platform-docs.opentargets.org/disease-or-phenotype
[efo]: https://www.ebi.ac.uk/efo/

<!-- graphql -->
[GraphQLClient]: https://github.com/mentorchains/BI-ML_Disease-Prediction_2024/blob/ahmed/src/gene_disease/datasets/open_targets/graphql_client.py#L21-L89
[graphql]: https://platform-docs.opentargets.org/data-access/graphql-api
[graphql-schema]: https://api.platform.opentargets.org/api/v4/graphql/schema
[graphql-ex]: https://platform-docs.opentargets.org/data-access/graphql-api#example-graphql-query
[graphql-playground-ex]: https://api.platform.opentargets.org/api/v4/graphql/browser?query=%23%20GraphQL%20Query%20example%20used%20in%20GDAP%0Aquery%20GetDiseaseTargets%20%7B%0A%20%20%23%20Disease%20example%3A%20cardiovascular%20disease%0A%20%20disease%28efoId%3A%20%22EFO_0000319%22%29%20%7B%0A%20%20%20%20id%0A%20%20%20%20name%0A%20%20%20%20associatedTargets%28page%3A%20%7B%20size%3A%2050%2C%20index%3A%201%20%7D%29%20%7B%0A%20%20%20%20%20%20count%0A%20%20%20%20%20%20rows%20%7B%0A%20%20%20%20%20%20%20%20target%20%7B%0A%20%20%20%20%20%20%20%20%20%20id%0A%20%20%20%20%20%20%20%20%20%20approvedSymbol%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%20%20score%0A%20%20%20%20%20%20%20%20datasourceScores%20%7B%0A%20%20%20%20%20%20%20%20%20%20id%0A%20%20%20%20%20%20%20%20%20%20score%0A%20%20%20%20%20%20%20%20%7D%0A%20%20%20%20%20%20%7D%0A%20%20%20%20%7D%0A%20%20%7D%0A%7D%0A

<!-- bigquery -->
[bigquery]: https://platform-docs.opentargets.org/data-access/google-bigquery
[bigquery-console]: https://console.cloud.google.com/marketplace/product/bigquery-public-data/open-targets-platform
[bigquery-blog]: https://www.bing.com/ck/a?!&&p=75ca412274f363d765e58ca8d6552f9c99150d093a89ead8f1a81b64ab97cbffJmltdHM9MTczMjc1MjAwMA&ptn=3&ver=2&hsh=4&fclid=3e17fb8c-d68a-68d3-3203-ee9dd72369bb&psq=open+Targets+Platform+bigquery+google&u=a1aHR0cHM6Ly9ibG9nLm9wZW50YXJnZXRzLm9yZy9uZXh0LWdlbi1wbGF0Zm9ybS1yZWxlYXNlZC8jOn46dGV4dD1Hb29nbGUlMjBCaWdRdWVyeSUyMGFsbG93cyUyMHlvdSUyMHRvJTIwcnVuJTIwYXN5bmNocm9ub3VzJTIwU1FMLEVUTCUyMHBpcGVsaW5lcyUyQyUyMHdpdGglMjBtb3JlJTIwdGhhbiUyMDIwJTIwZGF0YXNldHMlMjBhdmFpbGFibGUu&ntb=1
[bigquery-docs-ex]: https://platform-docs.opentargets.org/data-access/google-bigquery#example-bigquery-sql-queries
[bigquery-db]: https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/24.06/output/etl/json/
[bigquery-console-ex]: https://console.cloud.google.com/bigquery?sq=119305062578:ae33e769ecec43ceb055565f7d6c74df

<!-- string database -->
[string-db]: https://string-db.org/
[string-ppi-network]: https://string-db.org/cgi/covid.pl
[string-combined-score]: https://string-db.org/help/faq/#how-are-the-scores-computed

<!-- nodevectors -->
[nodevectors]: https://github.com/VHRanger/nodevectors