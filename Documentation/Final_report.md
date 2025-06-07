**Alex Chilton**  
<span style="color: gray; font-size: 0.8em;">(addresse)</span>   
<span style="font-size: 0.8em;">alex_chilton@gmx.co.uk</span>   

**Lauro Foletti**  
<span style="color: gray; font-size: 0.8em;">Rue de la Balance 1</br>CH &mdash; 2000 Neuchâtel</span>   
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

We present a comprehensive methodological exploration of computational approaches for generating plausible three-dimensional protein structures from PDB-derived representations. While originally conceived for Physics-Informed Neural Network applications targeting environmental perturbation prediction, our research evolved to address fundamental challenges in protein structure generation methodologies. Our work encompassed systematic data preprocessing using Linux-based computational biology libraries, followed by extensive exploration of multiple generative approaches including contact map-based Variational Autoencoders with gradient optimization, Graph Attention Network architectures, synthetic graph validation frameworks, and Graph Recurrent Attention Network-inspired dual output systems. Our investigation demonstrates the complexity of computational protein generation while establishing methodological foundations for future Physics-Informed Neural Network development. The most successful approach achieved 99.51% contact map accuracy with visually plausible 3D structures, providing strong foundations for advancing toward environmental perturbation modeling applications.

</div>

<div style="page-break-after: always;"></div>

## Table of contents  

- [Abstract](#abstract)
- [Table of contents](#table-of-contents)
- [1. Introduction and project scope](#1-introduction-and-project-scope)
  - [1.1 Biological and computational context](#11-biological-and-computational-context)
  - [1.2 Project objectives and evolution](#12-project-objectives-and-evolution)
  - [1.3 Methodological framework](#13-methodological-framework)
- [System requirements](#system-requirements)
- [Data](#data)
  - [Webscrapring](#webscrapring)
  - [Synthetic datasets](#synthetic-datasets)
  - [Protein Structure Preprocessing](#protein-structure-preprocessing)
- [Data quality](#data-quality)
- [Protein structure prediction](#protein-structure-prediction)
  - [General architecture](#general-architecture)
  - [Optimization Pipeline](#optimization-pipeline)
    - [MDS pipeline](#mds-pipeline)
    - [Gradient-Based Optimization](#gradient-based-optimization)
  - [Empirical Parameter Optimization](#empirical-parameter-optimization)
  - [Datasets](#datasets)
  - [Model](#model)
  - [Results](#results)
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

</div>

## System requirements

<div style="text-align: justify;">

The preprocessing pipeline was developed and executed under the following system configuration and software dependencies:

- Operating System: Linux-based environment

- Python Version: 3.9

- Required Python Packages: numpy ≥ 1.20, BioPython ≥ 1.79 (including Bio.PDB and DSSP interface), Graphein (for protein graph construction), networkx ≥ 2.6

- Hardware: CPU-based system with minimum 30 GB of available RAM to support memory-efficient feature extraction for large PDB files

- Parallel Execution: A parallelized version of the pipeline is available and was deployed on a high-performance computing (HPC) cluster to process datasets in parallel across multiple nodes.

**Note**: A Linux-based system is required due to dependencies on external tools used by DSSP and Graphein, which may rely on Unix-specific system calls, executable paths, or subprocess handling that are not reliably supported on Windows or macOS environments.

</div>

## Data

<div style="text-align: justify;">

### Webscrapring

Data were obtained from the Protein Data Bank (PDB) [[1]](#ref1) using the official API. The initial working dataset was restricted to nanobodies, selected for their relatively uniform length (typically 100–150 amino acids) and substantial structural diversity. In a secondary phase, the pipeline was extended to a second dataset comprising diverse protein types and subsequently cross-referenced with the BRENDA enzyme database [[2]](#ref2) to retrieve available experimental pH annotations for future analyses involving pH-dependent structural features.

### Synthetic datasets

Additionally to the webscraped data two sets of synthetic data were created to compare the generative capabilities of the methods using data with known properties. 

### Protein Structure Preprocessing

A memory-efficient preprocessing pipeline was proposed for the systematic extraction of structural, geometric, and physicochemical features from Protein Data Bank (PDB) files. The implementation was designed as a modular system and built to handle large-scale datasets with minimal memory overhead. All files with the `.pdb` extension were identified and validated by size $S$. Only those satisfying the condition:

$$S \leq S_{\text{max}} \quad \text{with} \quad S_{\\text{max}} = 25\,\text{MB}$$
were processed.

Using BioPython, the structures were parsed at the atomic level, where each standard residue was extracted and stored by chain as a tuple representation of the form:

$$
a_i = (\text{chain\_id}, r_i, t_i, a_i^n, e_i, \vec{x}_i)
$$

where:
- $r_i$ is the residue number,
- $t_i$ is the residue type (three-letter code),
- $a_i^n$ is the atom name,
- $e_i$ is the element,
- $\vec{x}_i \in \mathbb{R}^3$ is the atomic position.

Residue-level secondary structures were computed using DSSP [[3]](#ref3), resulting in a mapping:

$$
(r_i, \text{chain\_id}) \mapsto s_i \in \{H, E, C, G, I, T, S, ?\}
$$

where $s_i$ denoted the DSSP-assigned secondary structure class.

Subsequently, residue-level protein graphs $G = (V, E)$ were built where $V$ represented nodes features and $E$ represented edge features (either peptide bonds or spatial proximity). Edges $(i, j) \in \mathbb{E}$ were defined by the following:

$$
\|\vec{x}_i^{CA} - \vec{x}_j^{CA}\|_2 < \delta \quad \text{with} \quad \delta = 7.0\,\text{\AA}
$$

To characterize the 3D conformation of protein chains, bond lengths, bond angles, and torsion angles were computed. These factors enables to capture local spatial relationships between atoms and provide fundamental features for protein structure and stability. Due to the computational cost associated with calculating all geometric interactions at the atomic level—and in line with the specific objectives of the present study— the calculations were restricted to the backbone atoms only. These include the nitrogen (N), alpha carbon (CA), carbon (C), and oxygen (O) atoms that form the peptide backbone of each amino acid residue. This reduction significantly improved computational efficiency while retaining essential structural information relevant for most downstream analyses. The underlying implementation, however, was been designed to support full-atom calculations and could be extended in future work. Each chain was processed independently.  
**Bond lengths** are defined as the Euclidean distance between two bonded atoms. In this context, they were calculated between pairs of backbone atoms within the same residue or between sequential residues. The bond length $d_{ij}$ between atoms $i$ and $j$ was computed by: 

  $$
  d_{ij} = \|\vec{x}_i - \vec{x}_j\|_2
  $$

  where $\vec{x}_i$ and $\vec{x}_j$ were the 3D coordinates of the respective atoms.  **Bond Angles** (for triplet  $i - j - k$) provide a measure of the angle formed by three consecutive atoms and are useful for assessing local bending in the protein chain. Given a triplet of atoms $(i, j, k)$ , the  bond angle $\theta_{ijk}$  were computed as:
  $$
  \theta_{ijk} = \cos^{-1}\left( \frac{ (\vec{x}_i - \vec{x}_j) \cdot (\vec{x}_k - \vec{x}_j) }{ \|\vec{x}_i - \vec{x}_j\|_2 \cdot \|\vec{x}_k - \vec{x}_j\|_2 } \right)
  $$

  Only triplets involving backbone atoms were considered, including intra-residue angles (e.g., N–CA–C) and inter-residue connections (e.g., C–N–CA across peptide bonds). **Torsion Angles**, also called dihedral angles, describe the rotational relationship between four sequential atoms and are critical for capturing the 3D folding of the protein. Also these angles were computed exclusively on the backbone atoms and included the three canonical torsions: phi $(\phi)$, psi $(\psi)$ and omega $(\omega)$, each corresponding to a specific atom configuration across sequential residues:  

  - $\phi$: from atoms $C_{i-1} - N_i -  CA_i -  C_i$, for the rotation on the $N-CA$ bond
  - $\psi$: from $N_i - CA_i - C_{i} – N_{i-1}$, for the the rotation on the $CA–C$ bond
  - $\omega$: from $CA_{i–1}–C_{i–1}–N_i–CA_i$, reflecting the cis/trans conformation across the peptide bond.    

For atoms $(i, j, k, l)$ , the torsion angle $\phi$ were calculated using the angle between the planes formed by $(i, j, k)$ and $(j, k, l)$:

  Let $\vec{b}_1 = \vec{x}_j - \vec{x}_i,  \vec{b}_2 = \vec{x}_k - \vec{x}_j,  \vec{b}_3 = \vec{x}_l - \vec{x}_k$

  $$
  \phi = \text{atan2}\left( \frac{ \vec{b}_2 \cdot (\vec{b}_1 \times \vec{b}_3) }{ \|\vec{b}_1 \times \vec{b}_2\|_2 \cdot \|\vec{b}_2 \times \vec{b}_3\|_2 }, (\vec{b}_1 \times \vec{b}_2) \cdot (\vec{b}_2 \times \vec{b}_3) \right)
  $$

  where $\vec{x}_i$, $\vec{x}_j$, $\vec{x}_k$, $\vec{x}_l$ were the 3D coordinates of the respective atoms. This formulation was applied uniformly across all torsion types.

Lastly the physiochemical features were computed. Charges were estimated using predefined rules:

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

The graph structure of the NetworkX graph obtained is summarized in [Table 1](#table-1)

<a id="table-1"></a>
**Table 1:** Protein Graph Representation

| Component | Attribute        | Description                                                                 |
|-----------|------------------|-----------------------------------------------------------------------------|
| Node      | `Node ID`        | Unique node identifier: `"chain_id:residue_name:residue_number"`            |
| Node      | `chain_id`       | Protein chain identifier                                                    |
| Node      | `residue_number` | Sequential position of the residue in the protein                           |
| Node      | `residue_name`   | Three-letter amino acid code                                                |
| Node      | `ss`             | Secondary structure class (DSSP code: H, E, B, G, T, S, ?)                   |
| Node      | `has_backbone`   | Boolean indicating presence of backbone atoms (N, CA, C, O)                 |
| Node      | `coords`         | 3D coordinates of the CA atom                                               |
| Node      | `CA`, `N`, `C`, `O` | Optional atomic coordinates if available                                    |
| Edge      | `kind`           | Type of connection: `{peptide_bond}` or `{contact}`                        |
| Edge      | `distance`       | Distance between CA atoms for contact edges (Ångströms)                    |


## Data quality

For the first stage opf the project a datsset composed of nanobodies only was used. While this gave substantial advantages it carried also some limitations. Nanobodies are small in size ( ~15kDa, ~120-130 residues) but they maintain a consistenmt immunoglobulin fold architecture. They are naturally evolved, stable proteins which have been selected for proper folding and function, allowing the models to learn from biologically viable conformation. PDB nanobodiy structures are high-resolution X-ray or cryo-EM structures of high quality, this provides stable ground trith data for trainingm which is crucial for physics-informed approaches of common usew in the field. Despite their small size, nanobodies exhibit wide sequence and binding diversity while maintaining the same overall fold, exposing the models to functional variations without the complexity of completely different protein architecture. 
While though nanobodies are diverse in sequence, they share the same immunoglobulin scaffold which could limit the generalizatiopn ability of the trained models, additionally many PDB nanobodies come from camelids or have been engineered for stability, potentially biasing the training data toward unusually stable conformations which may not represent typical protein behaviour. This could be a potential drwback for a physics informed neural network but it's not an applicable risk for this stage of the project focused primarily on the investigation of the different generation techniques from data preprocessed for a physics informed network. 

## Protein Structure Prediction

### General architecture

### Optimization Pipeline

The preliminary encoding of distance matrices into contact matrices involves a significant loss of information. Continuous values are eliminated in favor of indicative binary values, reducing them to very rough approximations. Under certain conditions, however, these approximations can be combined in such a way as to recover close original values, making the encoding process partially - and surprisingly, reversible. The binary and categorical aspect of contact matrices can be directly assimilated to probability mass distributions, enabling the use of sigmoids for optimization.

Interesting results were obtained with combining two operations: 

1.	Logit initialization. A first matrix of distances is produced by taking a uniform random value within the boundaries of the respective domains expressed by the contacts; the matrix is then reduced to three dimensions via a classical MDS (Multidimensional Scaling); these three values constitute a first approximation of the 3D coordinates and are used to initialize a table of logits;

2.	Optimization. Logits are fed into an optimization loop, which translates them into relative distances and compares them with contact information. This optimization process relies entirely on a double sigmoid which “validates” the proposed values when they fall within the expected range, and penalizes them when they fall outside it, all in a continuous and differentiable manner (soft ranges). Proposed values are compared with target values by BCE.

3.	The process returns a set of three-dimensional coordinates whose relative distances correspond as closely as possible to the contact maps.

Sample preliminary tests were carried out with three distinct spatial structures: a homogeneous point cloud, a heterogeneous cloud and a nanoprotein. These structures are identical in length (120 points), and comparable in scale (approx. 35 distance units). Their spatial distribution differs: uniform, clustered and “organic”.

The specific distribution of each structure leads us to consider several approaches to domain partitioning : 
<ol type="a">
  <li><u>Percentiles</u>. A quantitative boundary ensures the homogeneous distribution of contacts between domains;</li>
  <li><u>Sectors</u>. Qualitative demarcation, distributing distances between sectors of equal width;</li>
  <li><u>Structural</u>. An “organic” demarcation based on the “ripples” visible on certain distribution curves, sensitive to the intrinsic structure of the data.</li>
</ol>

#### MDS pipeline

Our analysis revealed that coordinate recovery was straightforward when complete distance matrices were available, with Multidimensional Scaling achieving extremely low error (MAE = 0.0000 for uniform structures) in very short process time. However, when using contact matrices with information loss, MDS relied on prototype distance matrices based on random sampling from contact domains. We found that:
- Partitioning into sectors of identical width gave optimal results
- Percentile partitioning followed closely in performance  
- Structural partitioning yielded least accurate results
  
![MDS reconstruction accuracy](figures/mds_accuracy_grid.png)

#### Gradient-Based Optimization

Our systematic analysis revealed significant improvements over MDS prototypes using Gradient-Based Optimization:
- **Error Reduction**: Mean absolute error greatly reduced with significant improvement between 0-20% coverage
- **Coverage Optimization**: For protein structures, performance plateaued around 60% coverage
- **Domain Number Impact**: More contact matrix domains consistently improved reconstruction quality
- **Partitioning Strategy Reversal**: Percentile partitioning proved most efficient for optimization, while structural partitioning was least effective

![Optimized reconstruction accuracy](figures/gbo_accuracy_grid.png)

### Empirical Parameter Optimization

Our comparative study established an empirical framework for critical process parameters:

*Chunk Length Optimization:*
- Optimal results achieved above 60% coverage
- Consequently, subsequence length should be 30% at least of maximum studied length
- Training dataset constraint: shortest chain length should correspond to 30% at least of longest chain

*Contact Domain Configuration:*
- 6-domain percentile partitioning identified as optimal
- Balances reconstruction accuracy with computational efficiency
- Provides consistent performance across different structure types

### Model

### Results

<img src="./figures/FLUOPROTEIN_result.png" alt="Fluoprotein reconstruction" style="width: 500px;">

## References

<a id="ref1"></a>
[1] Berman, H. M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T. N., Weissig, H., Shindyalov, I. N., & Bourne, P. E. (2000). The Protein Data Bank. Nucleic Acids Research, 28(1), 235-242. https://doi.org/10.1093/nar/28.1.235

<a id="ref2"></a>
[2] Chang A., Jeske L., Ulbrich S., Hofmann J., Koblitz J., Schomburg I., Neumann-Schaal M., Jahn D., Schomburg D. BRENDA, the ELIXIR core data resource in 2021: new developments and updates. (2021), Nucleic Acids Res., 49:D498-D508. DOI: 10.1093/nar/gkaa1025 PubMed: 33211880

<a id="ref3"></a>
[3] Kabsch, W., & Sander, C. (1983). Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22(12), 2577–2637. https://doi.org/10.1002/bip.360221211
