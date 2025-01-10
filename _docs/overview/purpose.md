---
layout: default
title: Purpose
parent: Overview
nav_order: 1

---

## GDAP purpose

Gene-Disease Association Prediction (GDAP) is a tool designed to predict gene-disease associations using machine learning and node embeddings. By integrating data from sources like the  [STRING Database][string-db] (for protein-protein interactions) and [Open Targets][opentargets-platform] (for disease-gene associations), GDAP constructs a bipartite graph where two distinct node types are represented: `disease-gene` and `protein-protein`. This graph captures relationships between genes and diseases, with edges representing known disease associations or protein-protein interactions. Node embeddings are generated using various algorithms to capture these relationships, which are then used for feature extraction to predict both known and novel gene-disease associations.

[string-db]: https://string-db.org/
[opentargets-platform]: https://platform.opentargets.org/
