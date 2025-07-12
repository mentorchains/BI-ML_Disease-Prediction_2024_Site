---
layout: default
title: File Formats
parent: Overview
nav_order: 3

---

## Input and output format

The input files for GDAP are in tables text format (`.CSV`) described in the [Data Collection](../workflow/#1-data-collection) section. Loaded with `pandas` as a dataframe.

<table style="background-color: white; border: 1px solid #ccc; border-collapse: collapse;">
  <thead>
    <tr>
      <th style="background-color: #f4f4f4;">File</th>
      <th style="background-color: #f4f4f4;">Description</th>
      <th style="background-color: #f4f4f4;">File Format</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="background-color: white;">Target Data</td>
      <td style="background-color: white; height: 100px;">The target disease data table consists of symbol Gene Name (Protein) and score type (direct, indirect, global). (<a href="../workflow/#11-open-targets-dataset" target="_blank">Read More</a>)</td>
      <td style="background-color: white;"> pd.Dataframe | .csv</td>
    </tr>
    <tr>
      <td style="background-color: white;">String Database</td>
      <td style="background-color: white; height: 100px;">A protein-protein interactions dataset consisting of GeneName1 (Protein), GeneName2 (Protein), and the combined score. (<a href="../workflow/#12-string-database" target="_blank">Read More</a>)</td>
      <td style="background-color: white;">pd.Dataframe | .csv</td>
    </tr>
    <tr>
      <td style="background-color: white;">GraphML Format</td>
      <td style="background-color: white; height: 100px;">An XML-based file format for storing graphs with nodes and edges. It supports metadata and complex attributes. (<a href="../workflow/#2-graph-construction" target="_blank">Read More</a>)</td>
      <td style="background-color: white;">networkx.Graph object | .graphml (XML)</td>
    </tr>
    <tr>
      <td style="background-color: white;">Edge List</td>
      <td style="background-color: white; height: 100px;">The edge list format consists of pairs of nodes representing the edges in the graph. This is a simple text-based format.</td>
      <td style="background-color: white;">.edgelist (Plain Text)</td>
    </tr>
    <tr>
      <td style="background-color: white;">Edges CSV</td>
      <td style="background-color: white; height: 100px;">A CSV file representing edges in a graph. It contains two columns: source and destination nodes.</td>
      <td style="background-color: white;">.csv (Comma-Separated Values)</td>
    </tr>
    <tr>
      <td style="background-color: white;">Positive Edges</td>
      <td style="background-color: white; height: 100px;">File storing positive edges in a graph. This is typically used in graph learning for positive relationships.</td>
      <td style="background-color: white;">.npy (NumPy Array)</td>
    </tr>
    <tr>
      <td style="background-color: white;">Negative Edges</td>
      <td style="background-color: white; height: 100px;">File storing negative edges in a graph, representing negative relationships or absence of a connection.</td>
      <td style="background-color: white;">.npy (NumPy Array)</td>
    </tr>
  </tbody>
</table>