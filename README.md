# Protein Structure Generation & Analysis — CAS AML Final Project
**University of Bern — CAS in Advanced Machine Learning**

An advanced computational biology project exploring deep learning approaches to protein structure generation, property prediction, and family classification. Built as a significant extension of prior work on Graph VAEs, this project implements multiple parallel generative and predictive tracks using the SCOP protein database.

---

## Overview

Proteins are represented as molecular graphs where nodes are C-alpha atoms and edges encode spatial proximity. The project investigates how well deep learning models can learn, generate, and reason about protein structure across four interconnected research tracks.

---

## Research Tracks

### 1. Asymmetric VAE — Contact Map Generation
- Custom **Asymmetric Variational Autoencoder (ASYM-VAE)** operating on protein contact maps (C-alpha distance matrices)
- Contact maps encoded as continuous representations using fuzzy banding and soft binning
- Gradient-based optimisation in latent space to refine generated structures
- Recovery of 3D coordinates from generated contact maps

### 2. Conditional Generation — pH-Optimised Proteins
- **Conditional VAE** conditioned on biochemical properties (pH optimum, isoelectric point)
- Generates protein representations optimised for target pH conditions
- Feature engineering incorporating hydrophobicity profiles, BLOSUM62 matrices, and charge calculations

### 3. Graph Neural Network Classification
- **GCN and GAT (Graph Attention Network)** models for SCOP protein family classification
- 7-class classification from structural graph representations
- Demonstrates graph representation learning applied to protein taxonomy

### 4. GRAN — Dual Generative Model
- **Graph Recurrent Attention Network (GRAN)** with dual outputs:
  - Protein contact maps (adjacency matrices)
  - Amino acid sequences
- Multi-head graph attention mechanisms for probabilistic structure generation
- Sampling-based generation with explicit sequence/structure co-design

---

## Key Techniques

| Area | Methods |
|------|---------|
| Protein representation | C-alpha extraction, contact maps, KNN graphs, distance thresholds |
| Feature engineering | Physicochemical properties, BLOSUM62, hydropathy indices, one-hot AA encoding |
| Generative models | VAE, Conditional VAE, GRAN |
| Discriminative models | GCNConv, GATConv, batch normalisation, dropout |
| Optimisation | KL-divergence, reconstruction loss, gradient-based latent optimisation |
| Data | SCOP protein database (thousands of PDB structures) |

---

## Repository Structure

```
├── notebooks/          # Experimental and analysis notebooks
├── models/             # VAE, GRAN, GNN model definitions
├── utils/              # Protein parsing, graph construction, feature extraction
├── protein_analyzer.py # C-alpha extraction, physicochemical feature computation
├── graph_creator.py    # Graph construction with one-hot encoding and padding
└── requirements.txt    # Dependencies
```

---

## Background

This project extends [CAS_Project_2_protein_autoencoders](https://github.com/alexchilton/CAS_Project_2_protein_autoencoders), which established the baseline Graph VAE architecture. The final project significantly expands scope to include conditional generation, dual generative modelling, and supervised classification — moving from unsupervised embedding towards controlled protein design.

---

## Course Context

**Programme:** CAS in Advanced Machine Learning, University of Bern  
**Focus:** Graph neural networks, generative modelling, computational biology  
**Key Libraries:** PyTorch, PyTorch Geometric, torch_geometric, scikit-learn, Biopython
