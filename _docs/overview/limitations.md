---
layout: default
title: Limitations
parent: Overview
nav_order: 3
---

# Limitations

While GDAP provides valuable predictions, it has several important limitations that users should be aware of.

## Data Limitations

### Coverage
- **Disease-specific**: Predictions are limited to diseases with sufficient data in Open Targets Platform
- **Gene coverage**: Only includes genes with protein interaction data in STRING database
- **Species-specific**: Currently limited to Homo sapiens

### Data Quality
- **Experimental bias**: Training data reflects historical research priorities
- **Publication bias**: Well-studied diseases have more data
- **Annotation quality**: Depends on accuracy of source databases

## Model Limitations

### Prediction Accuracy
- **False positives**: Some predicted associations may not be biologically relevant
- **False negatives**: May miss associations due to data limitations
- **Confidence scores**: Should not be interpreted as absolute probabilities

### Generalization
- **Training bias**: Models trained on known associations may not generalize to novel relationships
- **Disease specificity**: Models trained on one disease may not apply to others
- **Temporal changes**: New discoveries may not be reflected in training data

## Technical Limitations

### Computational Resources
- **Memory requirements**: Large datasets require significant RAM
- **Processing time**: Complex embeddings can take hours to generate
- **Storage**: Results can be several GB for large analyses

### Scalability
- **Graph size**: Very large networks may exceed memory limits
- **Embedding quality**: Performance may degrade with extremely large graphs
- **Model complexity**: Some algorithms don't scale well to massive datasets

## Biological Limitations

### Context Dependence
- **Tissue specificity**: Predictions don't account for tissue-specific expression
- **Temporal dynamics**: Don't capture developmental or disease progression changes
- **Environmental factors**: Don't consider external influences on gene expression

### Validation Requirements
- **Experimental validation**: All predictions require laboratory confirmation
- **Biological context**: Predictions need interpretation in disease context
- **Clinical relevance**: Not all associations are clinically meaningful

## Best Practices

### Mitigating Limitations
1. **Use multiple models**: Compare predictions across different algorithms
2. **Validate predictions**: Always confirm with experimental data
3. **Consider context**: Interpret results in biological and clinical context
4. **Update regularly**: Retrain models with new data

### Quality Control
- **Check data sources**: Verify disease and gene data quality
- **Monitor performance**: Track model accuracy over time
- **Document assumptions**: Clearly state limitations in reports