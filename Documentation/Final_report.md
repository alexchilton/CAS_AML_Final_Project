**Alex Chilton**  
<span style="color: gray; font-size: 0.8em;">(addresse)</span>   
<span style="font-size: 0.8em;">alex_chilton@gmx.co.uk</span>   

**Lauro Foletti**  
<span style="color: gray; font-size: 0.8em;">(addresse)</span>   
<span style="font-size: 0.8em;">lauro.foletti@gmail.com</span> 

**Lara Nonis**  
<span style="color: gray; font-size: 0.8em;">(addresse)</span>   
<span style="font-size: 0.8em;">lara.nonis01@gmail.com</span>  



<div style="height: 15vh"></div>

<div align="left">

<p style="color: #666; font-size: 18px; margin-bottom: 10px;">Data Science Project</p>

<h1 style="font-size: 35px; line-height: 1.2; margin: 20px 0; border-bottom: none; text-decoration: none;">
Methodological Exploration of Protein Structure Generation from PDB Representations
</h1>

<p style="color: #666; font-size: 16px; ">15 June 2025</p>

</div>


<div style="page-break-after: always;"></div>

## Abstract
<div style="text-align: justify;">

We present a comprehensive methodological exploration of computational approaches for generating plausible three-dimensional protein structures from PDB-derived representations. While originally conceived for Physics-Informed Neural Network applications targeting environmental perturbation prediction, our research evolved to address fundamental challenges in protein structure generation methodologies. Our work encompassed systematic data preprocessing using computational biology libraries, followed by extensive exploration of multiple generative approaches including contact map-based Variational Autoencoders with gradient optimization, Graph Attention Network architectures, synthetic graph validation frameworks, and Graph Recurrent Attention Network-inspired dual output systems. Our investigation demonstrates the complexity of computational protein generation while establishing methodological foundations for future Physics-Informed Neural Network development. The most successful approach achieved 99.51% contact map accuracy with visually plausible 3D structures, providing strong foundations for advancing toward environmental perturbation modeling applications.

</div>

<div style="page-break-after: always;"></div>

## Table of contents  

- [Abstract](#abstract)
- [Table of contents](#table-of-contents)
- [1. Introduction and project scope](#1-introduction-and-project-scope)
  - [1.1 Biological and computational context](#11-biological-and-computational-context)
  - [1.2 Project objectives and evolution](#12-project-objectives-and-evolution)
  - [1.3 Methodological framework](#13-methodological-framework)
- [2. Infrastructure and development environment](#2-infrastructure-and-development-environment)
- [3. Data](#3-data)
  - [3.1 Primary data source](#31-primary-data-source)
  - [3.2 Preprocessing and feature engineering](#32-preprocessing-and-feature-engineering)
  - [3.3 Geometric feature engineering and structural analysis](#33-geometric-feature-engineering-and-structural-analysis)
  - [3.4 Physicochemical property integration:](#34-physicochemical-property-integration)
  - [3.5 Graph structure and feature integration](#35-graph-structure-and-feature-integration)
  - [3.6 Additional datasets for model validation and testing](#36-additional-datasets-for-model-validation-and-testing)
  - [3.7 Data quality](#37-data-quality)
- [4. Methodological exploration](#4-methodological-exploration)
  - [4.1 Simple graph attention network VAE architecture](#41-simple-graph-attention-network-vae-architecture)
  - [4.3 Graph neural network exploration and validation framework](#43-graph-neural-network-exploration-and-validation-framework)
- [References](#references)

<div style="page-break-after: always;"></div>

## 1. Introduction and project scope

### 1.1 Biological and computational context

<div style="text-align: justify;">

Proteins are the most abundant biological macromolecules, present in all living creatures. They exhibit substantial diversity in both size and function, ranging from small peptides to large polymeric complexes. Functionally, proteins are central to virtually all biological processes and represent the primary end products of genetic information pathways. Structurally, proteins are linear polymers of amino acids linked by peptide bonds. Twenty standard α-amino acids constitute protein sequences, each containing a central (α) carbon bonded to an amino group, a carboxyl group, a hydrogen atom, and a variable side chain (R group). These side chains differ in size, structure, polarity, and charge, influencing the solubility and conformational properties of the protein. Except for glycine, all amino acids have four distinct substituents on the α-carbon (CA), imparting chirality to the molecule. The repeating sequence of α-carbons forms the protein backbone.  

Protein structure is conventionally described at four hierarchical levels. The primary structure is the linear sequence of the amino acids. The secondary structure refers to the local conformations such as α-helices and β-sheets. The tertiary structure is the full three-dimensional configuration of a single polypeptide chain, while the quaternary structure represents the assembly of multiple folded subunits into a functional complex. Unlike the primary structure, which is genetically encoded, the secondary and tertiary structures are influenced also by environmental factors. Predicting these higher-order structures remains a complex task due to the intricate balance of biochemical constraints and conformational variability.  

### 1.2 Project objectives and evolution

We initially conceptualized this project to investigate and compare computational methodologies for generating plausible three-dimensional protein structures from PDB-derived representations, with particular focus on Physics-Informed Neural Network applications for environmental perturbation prediction. Our original scope envisioned developing systems capable of predicting protein structural modifications under environmental changes such as pH variations, temperature fluctuations, and ionic strength modifications.

However, during implementation, our research evolved to address preliminary challenges in protein structure generation methodologies. The complexity of basic protein generation necessitated comprehensive exploration of underlying computational approaches before advancing to environmental perturbation modeling. This evolution led to systematic investigation of multiple generative architectures, with each approach addressing different aspects of the protein structure generation problem.

### 1.3 Methodological framework 

Our investigation was structured with clearly defined phases:

1. **Systematic Data Preprocessing**: Development of comprehensive pipelines for converting protein structural data into computationally tractable representations suitable for machine learning applications

2. **Multiple Methodological Exploration**: Investigation of various neural network architectures and approaches, building upon the shared preprocessing foundation

3. **Validation and Assessment**: Integration of findings from different approaches to identify optimal methodologies and establish foundations for future development

Our objective was to characterize the respective modeling capabilities, limitations, and potential integration strategies of different approaches for future protein structure prediction frameworks, while acknowledging that environmental modification prediction remains a target for future research.

</div>

## 2. Infrastructure and development environment

<div style="text-align: justify;">

**Core framework**: We used PyTorch as the primary deep learning framework, with PyTorch Geometric providing specialized graph neural network operations for protein graph processing. A Linux-based environment was required due to dependencies on external tools used by DSSP [[1]](#ref1) and Graphein [[2]](#ref2). For implementation, we leveraged Google Colab [[3]](#ref3) with GPU acceleration (Tesla T4, 15GB memory) and UBELIX [[4]](#ref4) for processing datasets across multiple nodes. We used different Python versions for different implementation phases: Python 3.9 for preprocessing and Python 3.11 or later for all other implementations, with required libraries detailed in requirements.txt.  

**Experiment tracking**: we used Weights & Biases [[5]](#ref5) integration for experiment logging and visualization, tracking training metrics, loss components, and model performance across different architectural approaches.   

**Data and Model Distribution**: We made the trained best model and associated metadata publicly available via Hugging Face [[6]](#ref6) to ensure reproducibility and facilitate further research. 

</div>

## 3. Data 

<div style="text-align: justify;">

### 3.1 Primary data source

Data were obtained from the Protein Data Bank (PDB) [[7]](#ref7) using the official API with carefully defined selection criteria.
The initial working dataset was restricted to nanobodies, selected for their relatively uniform length (typically 100–150 amino acids) and substantial structural diversity. In a secondary phase, we extended the pipeline to a second dataset comprising diverse protein types and subsequently cross-referenced with the BRENDA enzyme database [[8]](#ref8) to retrieve available experimental pH annotations for future analyses involving pH-dependent structural features. 

### 3.2 Preprocessing and feature engineering 

Using BioPython [[9]](#ref9), we parsed the structures at the atomic level, where each standard residue was extracted and stored by chain as a tuple representation of the form: 

$$
a_i = (\text{chain\_id}, r_i, t_i, a_i^n, e_i, \vec{x}_i)
$$

where:
- $r_i$ is the residue number,
- $t_i$ is the residue type (three-letter code),
- $a_i^n$ is the atom name,
- $e_i$ is the element,
- $\vec{x}_i \in \mathbb{R}^3$ is the atomic position.

Residue-level secondary structures were computed using DSSP [[4]](#ref4), resulting in a mapping:

$$
(r_i, \text{chain\_id}) \mapsto s_i \in \{H, E, C, G, I, T, S, ?\}
$$

where $s_i$ denoted the DSSP-assigned secondary structure class.

Subsequently, residue-level protein graphs $G = (V, E)$ were built where $V$ represented nodes features and $E$ represented edge features (either peptide bonds or spatial proximity). We defined Edges $(i, j) \in \mathbb{E}$ by the following:

$$
\|\vec{x}_i^{CA} - \vec{x}_j^{CA}\|_2 < \delta \quad \text{with} \quad \delta = 7.0\,\text{\AA}
$$

### 3.3 Geometric feature engineering and structural analysis

To characterize the 3D conformation of th protein chain we computed bond lengths, bond angles and torsion angles. These calculations were restricted to backbone atoms only (N, CA, C, O) to balance computational efficiency with essential structural information retention, though the underlying implementation was designed to support full-atom calculations for future extension.

Bond lengths were defined as Euclidean distances between bonded atoms:

$$
d_{ij} = \|\vec{x}_i - \vec{x}_j\|_2
$$

calculated between pairs of backbone atoms within the same residue or between sequential residues

For triplets of atoms (i, j, k), bond angles $\theta_{ijk}$ provided measures of angles formed by three consecutive atoms:

$$
\theta_{ijk} = \cos^{-1}\left( \frac{ (\vec{x}_i - \vec{x}_j) \cdot (\vec{x}_k - \vec{x}_j) }{ \|\vec{x}_i - \vec{x}_j\|_2 \cdot \|\vec{x}_k - \vec{x}_j\|_2 } \right)
$$

We considered only triplets involving backbone atoms, including intra-residue angles and inter-residue connections.

Torsion angles (dihedral angles) described rotational relationships between four sequential atoms, critical for capturing 3D folding patterns. We computed all the three canonical backbone torsions:

  - $\phi$: Rotation on the $N-CA$ bond
  - $\psi$: Rotation on the $CA–C$ bond
  - $\omega$: from $CA_{i–1}–C_{i–1}–N_i–CA_i$, reflecting the cis/trans conformation across the peptide bond. 

For atoms $(i, j, k, l)$ , the torsion angle $\phi$ were calculated using the angle between the planes formed by $(i, j, k)$ and $(j, k, l)$:

Let $\vec{b}_1 = \vec{x}_j - \vec{x}_i,  \vec{b}_2 = \vec{x}_k - \vec{x}_j,  \vec{b}_3 = \vec{x}_l - \vec{x}_k$

$$
\phi = \text{atan2}\left( \frac{ \vec{b}_2 \cdot (\vec{b}_1 \times \vec{b}_3) }{ \|\vec{b}_1 \times \vec{b}_2\|_2 \cdot \|\vec{b}_2 \times \vec{b}_3\|_2 }, (\vec{b}_1 \times \vec{b}_2) \cdot (\vec{b}_2 \times \vec{b}_3) \right)
$$

where $\vec{x}_i$, $\vec{x}_j$, $\vec{x}_k$, $\vec{x}_l$ were the 3D coordinates of the respective atoms. This formulation was applied uniformly across all torsion types.

### 3.4 Physicochemical property integration:

Charges were estimated using predefined rules based on atom types:

$$
q_i = 
\begin{cases}
-0.5 & \text{if } e_i = \text{N} \\
0.5 & \text{if } e_i = \text{C} \\
-0.5 & \text{if } e_i = \text{O} \\
0 & \text{otherwise}
\end{cases}
$$

Hydrophobic residues were identified using a fixed set $H$, such that:

$$
\text{hydrophobic}(r_i) = \begin{cases}
1 & \text{if } t_i \in H \\
0 & \text{otherwise}
\end{cases}
$$

### 3.5 Graph structure and feature integration

The resulting NetworkX graph structure incorporated comprehensive protein representations as summarized in the following node and edge attribute framework:

**Node Attributes:**  
- **Node ID**: Unique identifier following the format "chain_id:residue_name:residue_number"  
- **Chain ID**: Protein chain identifier
- **Residue Number**: Sequential position in protein
- **Residue Name**: Three-letter amino acid code
- **Secondary Structure**: DSSP classification (H, E, B, G, T, S, ?)
- **Backbone Presence**: Boolean indicating availability of backbone atoms (N, CA, C, O)
- **Coordinates**: 3D coordinates of CA atom
- **Optional Atomic Coordinates**: Full atomic coordinates when available

**Edge Attributes:**  
- **Connection Type**: Classification as {peptide_bond} or {contact}
- **Distance**: C-alpha distances for contact edges (measured in Ångströms)

### 3.6 Additional datasets for model validation and testing 

In addition to the nanobody dataset, which we ultimately used to develop and test our models, we employed several additional datasets to enhance model comprehension and assess methodological limits.  
**Synthetic Graph Validation Dataset**:To isolate our core methodology from biological complexity, we developed a controlled validation framework using synthetic graphs with known, controllable properties. We generated 2,500 synthetic graphs across five topological structures (circles, stars, grids, crosses, and line structures) with 8-20 nodes per graph. Node features included amino acid type encoding (21-dimensional one-hot), color features (5 categories), physicochemical properties (size, charge, hydrophobicity), and structural features (coordinates, degree, clustering coefficient).  
**Biological Validation Datasets**: We also utilized the SCOP dataset [[10]](#ref10) for diagnostic classification experiments to assess our graph representation capabilities, and an additional dataset of fluorescent proteins (which we termed "fluobodies") from the protein data bank.

### 3.7 Data quality

Our nanobody-focused approach provided several methodological advantages including size consistency (~15kDa, 120-130 residues), high-resolution structural data from X-ray/cryo-EM sources, and consistent immunoglobulin fold architecture while maintaining functional diversity. However, the shared scaffold architecture and potential engineering bias toward stable conformations may limit generalization to diverse protein families.
To support our methodological development, we employed additional datasets for different purposes: synthetic graphs with controlled topological properties to test our algorithms on simpler, known structures before tackling complex biological data, SCOP dataset entries for diagnostic classification experiments, and fluorescent proteins (fluobodies) as an alternative biological dataset. These limitations were considered acceptable for our methodological exploration phase focused on investigating generation techniques rather than comprehensive biological coverage.

</div>

## 4. Methodological exploration 

Building upon our systematic preprocessing foundation, we explored multiple computational approaches to protein structure generation. These explorations encompassed both theoretical innovations and practical implementations, ranging from contact map-based representations to advanced graph neural network architectures. Each approach was designed to address specific challenges in protein structure generation while building upon insights gained from previous attempts.

### 4.1 Simple graph attention network VAE architecture




### 4.3 Graph neural network exploration and validation framework

Building on our preprocessed graph representations, and in parallel to the feature map exploration we 
 
VAEs have shown success in molecular generation tasks, and the 
graph structure seemed well-suited to capture protein spatial relationships. 







## References

<a id="ref1"></a>
[1] Kabsch, W., & Sander, C. (1983). Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22(12), 2577–2637. https://doi.org/10.1002/bip.360221211

<a id="ref2"></a>
[2] Graphein - a Python Library for Geometric Deep Learning and Network Analysis on Protein Structures
Arian R. Jamasb, Pietro Lió, Tom L. Blundell
bioRxiv 2020.07.15.204701; doi: https://doi.org/10.1101/2020.07.15.204701

<a id="ref3"></a>
[3] Google Colab. (2017). Colaboratory: A Google research project. Google Research. https://colab.research.google.com/

<a id="ref4"></a>
[4] UBELIX (https://www.id.unibe.ch/hpc), the HPC cluster at the University of Bern.

<a id="ref5"></a>
[5] Biewald, L. (2020). Experiment tracking with Weights and Biases. Software available from wandb.com. URL: https://www.wandb.com/

<a id="ref6"></a>
[6] Casual Labs. (2025). GRAN-Nanobody-Proteins Dataset. Hugging Face. https://huggingface.co/Casual-Labs/gran-nanobody-proteins

<a id="ref7"></a>
[7] Berman, H. M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T. N., Weissig, H., Shindyalov, I. N., & Bourne, P. E. (2000). The Protein Data Bank. Nucleic Acids Research, 28(1), 235-242. https://doi.org/10.1093/nar/28.1.235

<a id="ref8"></a>
[8] Chang A., Jeske L., Ulbrich S., Hofmann J., Koblitz J., Schomburg I., Neumann-Schaal M., Jahn D., Schomburg D. BRENDA, the ELIXIR core data resource in 2021: new developments and updates. (2021), Nucleic Acids Res., 49:D498-D508. DOI: 10.1093/nar/gkaa1025 PubMed: 33211880

<a id="ref9"></a>
[9] Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3 https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878

<a id="ref10"></a>
[10] Antonina Andreeva, Dave Howorth, Cyrus Chothia, Eugene Kulesha, Alexey Murzin, SCOP2 prototype: a new approach to protein structure mining. (2014) Nucl. Acid Res., 42 (D1): D310-D314 and Antonina Andreeva, Eugene Kulesha, Julian Gough, Alexey Murzin, The SCOP database in 2020: expanded classification of representative family and superfamily domains of known protein structures. (2020) Nucl. Acid Res., 48 (D1): D376-D382

