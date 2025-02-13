# Spatial Transcriptomics Iterative Hierarchical Clustering (stIHC): A Novel Method for Identifying Spatial Gene Co-Expression Modules

Recent advancements in spatial transcriptomics technologies allow researchers to 
simultaneously measure RNA expression levels for hundreds to thousands of genes
while preserving spatial information within tissues, providing critical insights
into spatial gene expression patterns, tissue organization, and gene functionality.
However, existing methods for clustering spatially variable genes (SVGs) into
co-expression modules often fail to detect rare or unique spatial expression patterns.
To address this, we present spatial transcriptomics iterative hierarchical
clustering (stIHC), a novel method for clustering spatially variable genes (SVGs)
into co-expression modules, representing groups of genes with shared spatial
expression patterns. Through three simulations and applications to spatial transcriptomics
datasets from technologies such as 10x Visium, 10x Xenium, and
Spatial Transcriptomics, stIHC outperforms clustering approaches used by popular
SVG detection methods, including SPARK, SPARK-X, MERINGUE, and
SpatialDE. Gene Ontology enrichment analysis confirms that genes within each
module share consistent biological functions, supporting the functional relevance
of spatial co-expression. Robust across technologies with varying gene numbers
and spatial resolution, stIHC provides a powerful tool for decoding the spatial
organization of gene expression and the functional structure of complex tissues.

Higgins, C., Li, J.J., Carey, M. Spatial Transcriptomics Iterative Hierarchical Clustering (stIHC): A Novel Method for Identifying Spatial Gene Co-Expression Modules (2025) https://arxiv.org/abs/2502.09574

# Overview 
This repository contains the R implementation of stIHC.

An example demonstrating how to use stIHC is provided in simulations.R, along with the data required to reproduce the simulations described in our paper.

# Data
The Data.zip file contains the datasets required to reproduce the simulations descirbed in Section 3.1, 3.2 and 3.3 of our paper.

# Code
The stIHC.R script provides the necessary functions for spatial transcriptomics data analysis. 

The input to stIHC is a dataset where:
1. The first two columns contain the xy-coordinates of each spatial spot.
2. The remaining columns contain the expression measurements for each gene at the corresponding locations.

The output of stIHC includes:
1. The cluster assignment for each gene.
2. The mean expression level of each cluster.
