---
layout: default
title: Configuration
parent: Tutorial
nav_order: 2
---

# Configuration

GDAP's behavior is controlled through the [`src/config.py`](https://github.com/ahmedsalim3/GDAP/blob/main/src/config.py) file. This guide explains how to customize the pipeline parameters.

## Configuration File

```python
class Config:
    DISEASE_ID = "EFO_0000319"  # Cardiovascular disease
    DATA_SOURCE = "GraphQLClient"  # GraphQLClient or BigQueryClient
    QUERY = INDIRECT_SCORES  # DIRECT_SCORES or INDIRECT_SCORES
    PPI_INTERACTIONS = 5000000  # Max PPI interactions
    NEGATIVE_TO_POSITIVE_RATIO = 10  # Negative/positive sample ratio
    EMBEDDING_MODE = "GGVec"  # Node2Vec, ProNE, GGVec, degree_avg
    MODEL_NAME = "Logistic_Regression"  # ML model to use
    TEST_SIZE = 0.2  # Train/test split ratio
    OUTPUT_DIR = PROJ_ROOT / "results/"
    CREDENTIALS = PROJ_ROOT / "configs/stemaway-d5a18133ff83.json"
```

## Key Parameters

### Disease Selection

**DISEASE_ID**: EFO (Experimental Factor Ontology) identifier for your target disease.

```python
# Examples
DISEASE_ID = "EFO_0000319"  # Cardiovascular disease
DISEASE_ID = "MONDO_0004979"  # Asthma
```

Find disease EFOs at the [Open Targets Platform](https://platform.opentargets.org/), You can search by name or by an Ontology ID.

### Data Source Configuration

**DATA_SOURCE**: Choose between GraphQL API and BigQuery.

```python
# GraphQL API (default, no credentials needed)
DATA_SOURCE = "GraphQLClient"

# BigQuery (requires Google Cloud credentials)
DATA_SOURCE = "BigQueryClient"
```

**QUERY**: For BigQuery, choose query type:
- `DIRECT_SCORES`: Direct disease-target associations
- `INDIRECT_SCORES`: Indirect associations with evidence

### Embedding Configuration

**EMBEDDING_MODE**: Choose embedding algorithm:

```python
EMBEDDING_MODE = "GGVec"        # Graph-based global vectors (default)
EMBEDDING_MODE = "Node2Vec"     # Random walk-based embeddings
EMBEDDING_MODE = "ProNE"        # Proximity-preserving embeddings
EMBEDDING_MODE = "degree_avg"   # Simple degree-based features
```

### Model Configuration

**MODEL_NAME**: Choose machine learning model:

```python
MODEL_NAME = "Logistic_Regression"  # Linear classifier (default)
MODEL_NAME = "Random_Forest"        # Ensemble tree model
MODEL_NAME = "SVM"                  # Support Vector Machine
MODEL_NAME = "Gradient_Boosting"    # Gradient boosting trees
```

### Data Processing

**PPI_INTERACTIONS**: Maximum number of protein-protein interactions to load:
```python
PPI_INTERACTIONS = 5000000  # 5M interactions (default)
PPI_INTERACTIONS = 1000000  # 1M interactions (faster, less memory)
```

**NEGATIVE_TO_POSITIVE_RATIO**: Balance between positive and negative samples:
```python
NEGATIVE_TO_POSITIVE_RATIO = 10  # 10:1 negative to positive (default)
NEGATIVE_TO_POSITIVE_RATIO = 5   # 5:1 ratio
```

**TEST_SIZE**: Train/test split ratio:
```python
TEST_SIZE = 0.2  # 20% for testing (default)
TEST_SIZE = 0.3  # 30% for testing
```

## Example Configurations

### Quick Test Configuration

```python
class Config:
    DISEASE_ID = "EFO_0000319"
    DATA_SOURCE = "GraphQLClient"
    PPI_INTERACTIONS = 1000000  # Reduced for faster processing
    EMBEDDING_MODE = "degree_avg"  # Simple embeddings
    MODEL_NAME = "Logistic_Regression"
    TEST_SIZE = 0.2
```

### High-Performance Configuration

```python
class Config:
    DISEASE_ID = "EFO_0000319"
    DATA_SOURCE = "BigQueryClient"  # For large datasets
    QUERY = INDIRECT_SCORES
    PPI_INTERACTIONS = 10000000  # More interactions
    EMBEDDING_MODE = "GGVec"
    MODEL_NAME = "Random_Forest"
    NEGATIVE_TO_POSITIVE_RATIO = 15
```

## BigQuery Setup

If using BigQuery, you need Google Cloud credentials:

1. Create a Google Cloud project
2. Enable BigQuery API
3. Create a service account and download JSON key
4. Place the JSON file in `configs/` directory
5. Update the `CREDENTIALS` path in config.py

## Environment Variables

You can also use environment variables for sensitive configuration:

```bash
export GOOGLE_APPLICATION_CREDENTIALS="path/to/credentials.json"
export DISEASE_ID="EFO_0000319"
```

## Google Cloud Setup
1. Create a service account in Google Cloud with the necessary roles (BigQuery User)
    - Go to the [Google Cloud Console](https://console.cloud.google.com/welcome?project=stemaway).
    - Select your project or create a new one if needed.
    - In the navigation menu, go to IAM & Admin > Service Accounts.
    - Create a new Service Account or use an existing one.
    - Grant the necessary roles (e.g., BigQuery User) to this service account ID.
    - Create a key for this service account in JSON format.

2. Download the service account JSON key.
    - After creating the service account and generating the key, download the JSON key file to your local machine. This file will be used for authentication.

3. Set Up Google Cloud Credentials Locally (using `.env` or `direct path`)

- **Option 1: Using `.env` File**

    - Create a `.env` file in your project's root directory (if it doesn't already exist).
    - Add the **full path** to your service account JSON key in the .env file:

        ```sh
        GOOGLE_APPLICATION_CREDENTIALS=/path/to/json/file
        ```
    - Install the `python-dotenv` package if it isn't already installed. You can do this by running:

        ```sh
        pip install python-dotenv
        ```
- **Option 2: Using the environment variable GOOGLE_APPLICATION_CREDENTIALS**

- You can set this in your terminal or provide the [`credentials_path`](./bigquery_client.py#L26) directly when initializing the client
    - To set the GOOGLE_APPLICATION_CREDENTIALS environment variable through terminal, run these commands:

        ```sh
        export GOOGLE_APPLICATION_CREDENTIALS='/path/to/json/file'
        ```

        ```sh
        echo $GOOGLE_APPLICATION_CREDENTIALS  # Verify it
        unset GOOGLE_APPLICATION_CREDENTIALS  # Unset it
        ```

    - Provide the credentials path directly in your code when initializing the BigQuery client

4. Run the [main file](../main.py) to fetch the desired disease dataset

## Deployment

- Add the key file to your local app secrets `.streamlit/secrets.toml`

```toml
# .streamlit/secrets.toml

[gcp_service_account]
type = "service_account"
project_id = "xxx"
private_key_id = "xxx"
private_key = "xxx"
client_email = "xxx"
client_id = "xxx"
auth_uri = "https://accounts.google.com/o/oauth2/auth"
token_uri = "https://oauth2.googleapis.com/token"
auth_provider_x509_cert_url = "https://www.googleapis.com/oauth2/v1/certs"
client_x509_cert_url = "xxx"
```

- Copy your app secrets to the cloud
    Go to the [app dashboard](https://share.streamlit.io/) and in the app's dropdown menu, click on Edit Secrets. Copy the content of secrets.toml into the text area.

## References
- [Connect Streamlit to Google BigQuery](https://docs.streamlit.io/knowledge-base/tutorials/databases/bigquery)