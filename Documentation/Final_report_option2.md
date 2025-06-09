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
(Title TBD)
</h1>

<p style="color: #666; font-size: 16px; ">15 June 2025</p>

</div>


<div style="page-break-after: always;"></div>

## Abstract
<div style="text-align: justify;">

Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum

</div>

<div style="page-break-after: always;"></div>

## Table of contents  

- [Abstract](#abstract)
- [Table of contents](#table-of-contents)
- [1. Introduction and project scope](#1-introduction-and-project-scope)
  - [1.1 Biological and computational context](#11-biological-and-computational-context)
  - [1.2 Project objectives and evolution](#12-project-objectives-and-evolution)
- [2. Infrastructure and development environment](#2-infrastructure-and-development-environment)
- [3. Data](#3-data)
- [4. Metadata](#4-metadata)
- [5. Data quality](#5-data-quality)
- [6. Data flow](#6-data-flow)
  - [6.1 Parsing](#61-parsing)
  - [6.2 Preprocessing and feature engineering](#62-preprocessing-and-feature-engineering)
    - [6.2.1 Geometric feature engineering and structural analysis](#621-geometric-feature-engineering-and-structural-analysis)
    - [6.2.2 Physicochemical property integration:](#622-physicochemical-property-integration)
    - [6.2.3 Graph structure and feature integration](#623-graph-structure-and-feature-integration)
  - [6.3 Additional datasets for model validation and testing](#63-additional-datasets-for-model-validation-and-testing)
- [7. Model flow](#7-model-flow)
  - [7.1 Graph Variational Autoencoder (GraphVAE)](#71-graph-variational-autoencoder-graphvae)
  - [7.2 Simple Variational AutoEncoder within a Hybrid Learning–Inversion Framework (VAE)](#72-simple-variational-autoencoder-within-a-hybrid-learning-inversion-framework-vae)
- [References](#references)
- [List of Figures](#list-of-figures)
- [List of Tables](#list-of-tables)
- [Appendix](#appendix)

<div style="page-break-after: always;"></div>

## 1. Introduction and project scope

### 1.1 Biological and computational context

<div style="text-align: justify;">

Proteins are the most abundant biological macromolecules, present in all living creatures. They exhibit substantial diversity in both size and function, ranging from small peptides to large polymeric complexes. Functionally, proteins are central to virtually all biological processes and represent the primary end products of genetic information pathways. Structurally, proteins are linear polymers of amino acids linked by peptide bonds. Twenty standard α-amino acids constitute protein sequences, each containing a central (α) carbon bonded to an amino group, a carboxyl group, a hydrogen atom, and a variable side chain (R group). These side chains differ in size, structure, polarity, and charge, influencing the solubility and conformational properties of the protein. Except for glycine, all amino acids have four distinct substituents on the α-carbon (CA), imparting chirality to the molecule. The repeating sequence of α-carbons forms the protein backbone[[1]](#ref1).  

Protein structure is conventionally described at four hierarchical levels. The primary structure is the linear sequence of the amino acids. The secondary structure refers to the local conformations such as α-helices and β-sheets. The tertiary structure is the full three-dimensional configuration of a single polypeptide chain, while the quaternary structure represents the assembly of multiple folded subunits into a functional complex [[1]](#ref1). Unlike the primary structure, which is genetically encoded, the secondary and tertiary structures are influenced also by environmental factors. Predicting these higher-order structures remains a complex task due to the intricate balance of biochemical constraints and conformational variability.  

### 1.2 Project objectives and evolution

We initally conceptualized this project to develop a phisically informed neural network (PINN) able to predict protein tertiary structure changes upon environmental perturbation. The development of this ultimate goal imposed though preliminary studies necessary to understand the complex structure of proteins and the challenges intertwined with the translation, modification and generation of these elaborated biological structures.   
More specifically the preliminary studies involved:   
- **Parsing of cristallographic / X-ray experimental data (PDB files)**: to extract meaningful information for the subsequent phases of the project;
- **Preprocessing and pre-calculation**: of all the properties potentially required to describe the 3D structure of the protein and could play a role in the interaction with the environment;
- **Models exploration**: with state of the art generative models commonly used in the chemical field or in the image processing environment. 
- **Development of a loss function**: to constrain physically the neural network and drive the generation through biologically plausible structures

The present report will focus on the pre-studies leaving the ultimate goal for future development on the project. 

</div>

## 2. Infrastructure and development environment

<div style="text-align: justify;">

We used PyTorch [[2]](#ref2) as the primary deep learning framework, with PyTorch Geometric [[3]](#ref3) providing specialized graph neural network operations for protein graph processing and NetworkX for complex graphs generation and handling [[4]](#ref4). A Linux-based environment was required due to dependencies on external tools used by DSSP [[5]](#ref5) and Graphein [[6]](#ref6), necessary for the preprocessing and feature engineering phase. For implementation, we leveraged Google Colab [[7]](#ref7) with GPU acceleration (Tesla T4, 15GB memory) and UBELIX [[8]](#ref8) for processing datasets across multiple nodes. We used different Python versions for different implementation phases: Python 3.9 for preprocessing and Python 3.11 or later for all other implementations, with required libraries detailed in requirements.txt.  

We used Weights & Biases [[9]](#ref9) for integration for experiment logging and visualization, tracking training metrics, loss components, and model performance across different architectural approaches. 

The trained best model wand associated metadata will be made publicly available via Hugging Face [[10]](#ref10) to ensure reproducibility and facilitate further research.

</div>

## 3. Data 

<div style="text-align: justify;">

Data were obtained from the RCSB Protein Data Bank (PDB) [[11]](#ref11) using the official API with carefully defined selection criteria. The initial working dataset was restricted to nanobodies, selected for their relatively uniform length (typically 100–150 amino acids) and substantial structural diversity. Search included multiple criteria to catch as many nanobodies as possible including among others the family (Camelidae), mentions (VHH) or labels. The script was designed to retrieve all matching strucures, handling pagination and downloading each PDB file (containing 3D protein structure data) individually.   
In a secondary phase, we extended the pipeline to a second dataset comprising diverse protein types and subsequently cross-referenced with the BRENDA enzyme database [[12]](#ref12) to retrieve available experimental pH annotations for future analyses involving pH-dependent structural features, envisioning subsequent phases of the project and the requirements of a PINN within a supervised learning framework. 

</div>

## 4. Metadata

<div style="text-align: justify;">

The preliminary studies were handled in indipendent notebook or python files stored in the respective GitHub folder. Datasets were stored in a Drive shared folder due to dimension cotriction upon GitHub upload. The final preprocessed dataset with the respective metadata of the best model will be stored in the Hugging Face repositry as explained in section 2. 


</div>

## 5. Data quality 

<div style="text-align: justify;">

Our nanobody-focused approach provided several methodological advantages including size consistency (~15kDa, 120-130 residues), high-resolution structural data from X-ray/cryo-EM sources, and consistent immunoglobulin fold architecture while maintaining functional diversity. However, the shared scaffold architecture and potential engineering bias toward stable conformations may limit generalization to diverse protein families.
To support our methodological development, we employed additional datasets for different purposes: synthetic graphs with controlled topological properties to test our algorithms on simpler, known structures before tackling complex biological data, SCOP dataset entries for diagnostic classification experiments, and fluorescent proteins (fluobodies) as an alternative biological dataset. These limitations were considered acceptable for our preliminary phase focused on investigating generation techniques, suitable preprocessing and a meaningful loss function. .

</div>

## 6. Data flow

<div style="text-align: justify;">

### 6.1 Parsing

PDB files were parsed using BioPython [[13]](#ref13). We evaluated several parsing strategies considering atoms only and their position, aminoacid types and the neighborhood informations (number of surrounding aminoacids with average and max distance). Additional small molecules and external groups were removed, to focus on the protein structure only. Distances and positions in the 3D space were computed through the alpha carbon mapped with the amino acid mass to have additional information on the residue in a specific position, with relatively small dataframes. The so obtained matrices were converted in graphs and compared to the original molecular representation by PyMOL [[14]](#ref14).

<div style="text-align: center; margin: 20px 0;">
  <img src="figures/figure_a.png" alt="First parsing" width="600">
  <p style="margin: 10px 40px; font-style: italic;">
    <strong>Figure 1:</strong> Comparison between molecular and graph-based representations of the same protein structure. Left: PyMOL  molecular visualization showing secondary structure elements. Right: NetworkX graph representation displaying nodes (residues) connected by spatial proximity edges. The graph representation captures overall protein connectivity but does not clearly preserve the distinct secondary structure organization visible in the molecular model, illustrating the information loss during graph conversion..
  </p>
</div>

Upon the first visualization we could observe among the issues the existance of several chains, not belonging to the same proteins but being clusers instead, which were treated by the preprocessing as one molecule. We therefore preprocessed the pdb files to have indipendent chains (each representing one unique nanobody) and converted in heterogeneous graphs.   

<div style="text-align: center; margin: 20px 0;">
  <img src="figures/figure_b.png" alt="Second parsing" width="600">
  <p style="margin: 10px 40px; font-style: italic;">
    <strong>Figure 2:</strong> Comparison of protein representation formats: NetworkX graph visualization (left) and PyMOL molecular structure (right) of the same nanobody. The graph format preserves connectivity and secondary structure patterns while simplifying atomic detail..
  </p>
</div>


Ultimately the proteins were preprocessed in different steps, identifying the alpha carbon first, mapping the amino-acid types and the respective coordinates and computing the non-covalent intereactions. The reconstruction can be visualized in the following 3D representation:

<div style="text-align: center; margin: 20px 0;">
  <img src="figures/figure_c.png" alt="3D protein structure" width="600">
  <p style="margin: 10px 40px; font-style: italic;">
    <strong>Figure 3:</strong> 3D visualization of a protein structure showing amino acid residues as colored spheres according to type, connected by backbone bonds (gray lines) and non-covalent interactions (pink dashed lines). The compact, folded structure demonstrates typical protein architecture with diverse amino acid types distributed throughout the 3D space.
  </p>
</div>

Upon success of the single chain strategy the same logic was implemented into the final preprocessing to graphs and associated dictionaries. 

### 6.2 Preprocessing and feature engineering 

Using BioPython, we parsed the structures at the atomic level, where each standard residue was extracted and stored by chain as a tuple representation of the form: 

$$
a_i = (\text{chain\_id}, r_i, t_i, a_i^n, e_i, \vec{x}_i)
$$

where:
- $r_i$ is the residue number,
- $t_i$ is the residue type (three-letter code),
- $a_i^n$ is the atom name,
- $e_i$ is the element,
- $\vec{x}_i \in \mathbb{R}^3$ is the atomic position.

Residue-level secondary structures were computed using DSSP, resulting in a mapping:

$$
(r_i, \text{chain\_id}) \mapsto s_i \in \{H, E, C, G, I, T, S, ?\}
$$

where $s_i$ denoted the DSSP-assigned secondary structure class.

Subsequently, residue-level protein graphs $G = (V, E)$ were built where $V$ represented nodes features and $E$ represented edge features (either peptide bonds or spatial proximity). We defined Edges $(i, j) \in \mathbb{E}$ by the following:

$$
\|\vec{x}_i^{CA} - \vec{x}_j^{CA}\|_2 < \delta \quad \text{with} \quad \delta = 7.0\,\text{\AA}
$$

#### 6.2.1 Geometric feature engineering and structural analysis

To characterize the 3D conformation of the protein chain we computed bond lengths, bond angles and torsion angles. These calculations were restricted to backbone atoms only (N, CA, C, O) to balance computational efficiency with essential structural information retention, though the underlying implementation was designed to support full-atom calculations for future extension.

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

#### 6.2.2 Physicochemical property integration:

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

#### 6.2.3 Graph structure and feature integration

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

### 6.3 Additional datasets for model validation and testing 

In addition to the nanobody dataset, which we ultimately used to develop and test our models, we employed several additional datasets to enhance model comprehension and assess methodological limits.  

**Synthetic Graph Validation Dataset**:To isolate our core methodology from biological complexity, we developed a controlled validation framework using synthetic graphs with known, controllable properties. We generated 2,500 synthetic graphs across five topological structures (circles, stars, grids, crosses, and line structures) with 8-20 nodes per graph. Node features included amino acid type encoding (21-dimensional one-hot), color features (5 categories), physicochemical properties (size, charge, hydrophobicity), and structural features (coordinates, degree, clustering coefficient).  

**Biological Validation Datasets**: We also utilized the SCOP dataset [[15]](#ref15) for diagnostic classification experiments to assess our graph representation capabilities, and an additional dataset of fluorescent proteins (which we termed "fluobodies") from the protein data bank.

</div>

## 7. Model flow

<div style="text-align: justify;">

Building upon our systematic preprocessing foundation, we explored multiple computational approaches to protein structure generation. These explorations encompassed both theoretical innovations and practical implementations, ranging from contact map-based representations to advanced graph neural network architectures. Each approach was designed to address specific challenges in protein structure generation while building upon insights gained from previous attempts.

### 7.1 Graph Variational Autoencoder (GraphVAE)

VAEs have shown success in molecular generation tasks [[16]](#ref16) [[17]](#ref17), and the graph structure seemed well-suited to capture protein spatial relationships. 

We built a custom GraphVAE following the standard variational auto-encoder structure with the following parameters: 
- input: 8-dimensional node features representing physicochemical properties
- hidden channels: 64-dimensional intermediate representations
- latent space: 16-dimensional bottleneck latent space (z)
- 4 attention heads for multi-head attention mechanism
- 1 dimensional edge attribute (distances)

The model was developed to handle variable-sized protein graphs without padding artifacts. The attention mechanism with 4 heads allowed the model to capture multiple types of relationships simultaneously within the protein structure.   

We built the loss function as a composite element including reconstruction loss ($\mathcal{L}$ _recon), KL divergence loss ($\mathcal{L}$ _KL ) and rthogonal regularization ($\mathcal{L}$_ortho ). The reconstruction loss employed task-specific weighting to balance off different protein properties encoded in the graph structure. The KL parameter included a $\beta$-VAE with cyclical annealing to balance reconstruction quality and latent space organization with a warmup period t of 25 epochs: 

$$
\beta(t) = 
\begin{cases}
\beta_{\text{min}} + (\beta_{\text{max}} - \beta_{\text{min}}) \cdot \dfrac{t}{t_{\text{warmup}}} & \text{if } t < t_{\text{warmup}} \\
\beta_{\text{max}} & \text{if } t \geq t_{\text{warmup}}
\end{cases}
$$

the orthogonal regularization enforced normalization in the latent space representation: 

$$
\mu_{\text{normalized}} = \mathrm{F.normalize}(\mu, p=2, \text{dim}=1)
$$

$$
C = \frac{\mu_{\text{normalized}}^\top \mu_{\text{normalized}}}{\text{batch\_size}}
$$

$$
\mathcal{L}_{\text{ortho}} = \lambda_{\text{ortho}} \cdot \left\| C - I \right\|_F^2
$$

Where:

- \( \lambda_{\text{ortho}} = 0.1 \): Orthogonality strength parameter  
- \( C \): Correlation matrix of normalized latent means  
- \( I \): Identity matrix  
- \( \left\| \cdot \right\|_F \): Frobenius norm 

The training involved a learning rate schedule with a reduction factor of 0.5 and 3-epochs patience. Gradient clipping was applied with max_norm=0.1 to prevent the exploding gradient issue. 

For the generator part we used latent sampling $\mathbf{z} \sim \mathcal{N}(0, \sigma^2 I)$ with temperature scaling (NOTE FOR LARA: continue from here)

**Main benefits**:  
Unlike traditional VAEs operating on fixed size inputs, our architecture handled variable-sized graphs through attention based pooling for size invariant encoding, dynamic batching with graph-level indices and flexible decoding supporting arbitrary output sizes. The model  was able to simultaneously learn node-level chemical properties and structural features with graph-level batch indices through the multi-attention heads.  

### 7.2 Simple Variational AutoEncoder within a Hybrid Learning–Inversion Framework (VAE)

**General architecture**:

**Optimization Pipeline**:

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

**MDS pipeline**:

Our analysis revealed that coordinate recovery was straightforward when complete distance matrices were available, with Multidimensional Scaling achieving extremely low error (MAE = 0.0000 for uniform structures) in very short process time. However, when using contact matrices with information loss, MDS relied on prototype distance matrices based on random sampling from contact domains. We found that:
- Partitioning into sectors of identical width gave optimal results
- Percentile partitioning followed closely in performance  
- Structural partitioning yielded least accurate results
  
![MDS reconstruction accuracy](figures/mds_accuracy_grid.png)

**Gradient-Based Optimization**:

Our systematic analysis revealed significant improvements over MDS prototypes using Gradient-Based Optimization:
- **Error Reduction**: Mean absolute error greatly reduced with significant improvement between 0-20% coverage
- **Coverage Optimization**: For protein structures, performance plateaued around 60% coverage
- **Domain Number Impact**: More contact matrix domains consistently improved reconstruction quality
- **Partitioning Strategy Reversal**: Percentile partitioning proved most efficient for optimization, while structural partitioning was least effective

![Optimized reconstruction accuracy](figures/gbo_accuracy_grid.png)

**Empirical Parameter Optimization**:

Our comparative study established an empirical framework for critical process parameters:

*Chunk Length Optimization:*
- Optimal results achieved above 60% coverage
- Consequently, subsequence length should be 30% at least of maximum studied length
- Training dataset constraint: shortest chain length should correspond to 30% at least of longest chain

*Contact Domain Configuration:*
- 6-domain percentile partitioning identified as optimal
- Balances reconstruction accuracy with computational efficiency
- Provides consistent performance across different structure types

**Results**:

<img src="./figures/FLUOPROTEINS_lengths.png" alt="Fluoprotein lengths" style="width: 700px;">
<img src="./figures/FLUOPROTEINS_distances.png" alt="Fluoprotein distances" style="width: 700px;">

<img src="./figures/FLUOPROTEINS%20-%20Original%20contact%20maps%20(ground%20truth).png" alt="Ground truth" style="width: 1200px;">
<img src="./figures/FLUOPROTEINS%20-%20Reconstructed%20contact%20maps.png" alt="Ground truth" style="width: 1200px;">
<img src="./figures/FLUOPROTEINS%20-%20Errors.png" alt="Ground truth" style="width: 1200px;">

<img src="./figures/FLUOPROTEIN_result.png" alt="Fluoprotein reconstruction" style="width: 700px;">

</div>

## References

<a id="ref1"></a>
[1] Nelson, David L., and Michael M. Cox. Lehninger Principles of Biochemistry. 6th ed. W. H. Freeman, 2012.

<a id="ref2"></a>
[2] Paszke, A., Gross, S., Massa, F., Lerer, A., Bradbury, J., Chanan, G., ... & Chintala, S. (2019). PyTorch: An imperative style, high-performance deep learning library. Advances in neural information processing systems, 32.

<a id="ref3"></a>
[3] Fey, M., & Lenssen, J. E. (2019). Fast graph representation learning with PyTorch Geometric. arXiv preprint arXiv:1903.02428.

<a id="ref4"></a>
[4] Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart, “Exploring network structure, dynamics, and function using NetworkX”, in Proceedings of the 7th Python in Science Conference (SciPy2008), Gäel Varoquaux, Travis Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15, Aug 2008

<a id="ref5"></a>
[5] Kabsch, W., & Sander, C. (1983). Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Biopolymers, 22(12), 2577–2637. https://doi.org/10.1002/bip.360221211

<a id="ref6"></a>
[6] Graphein - a Python Library for Geometric Deep Learning and Network Analysis on Protein Structures
Arian R. Jamasb, Pietro Lió, Tom L. Blundell
bioRxiv 2020.07.15.204701; doi: https://doi.org/10.1101/2020.07.15.204701

<a id="ref7"></a>
[7] Google Colab. (2017). Colaboratory: A Google research project. Google Research. https://colab.research.google.com/

<a id="ref8"></a>
[8] UBELIX (https://www.id.unibe.ch/hpc), the HPC cluster at the University of Bern.  

<a id="ref9"></a>
[9] Biewald, L. (2020). Experiment tracking with Weights and Biases. Software available from wandb.com. URL: https://www.wandb.com/

<a id="ref10"></a>
[10] Casual Labs. (2025). GRAN-Nanobody-Proteins Dataset. Hugging Face. https://huggingface.co/Casual-Labs/gran-nanobody-proteins

<a id="ref11"></a>
[11] Berman, H. M., Westbrook, J., Feng, Z., Gilliland, G., Bhat, T. N., Weissig, H., Shindyalov, I. N., & Bourne, P. E. (2000). The Protein Data Bank. Nucleic Acids Research, 28(1), 235-242. https://doi.org/10.1093/nar/28.1.235

<a id="ref12"></a>
[12] Chang A., Jeske L., Ulbrich S., Hofmann J., Koblitz J., Schomburg I., Neumann-Schaal M., Jahn D., Schomburg D. BRENDA, the ELIXIR core data resource in 2021: new developments and updates. (2021), Nucleic Acids Res., 49:D498-D508. DOI: 10.1093/nar/gkaa1025 PubMed: 33211880

<a id="ref13"></a>
[13] Cock, P.J.A. et al. Biopython: freely available Python tools for computational molecular biology and bioinformatics. Bioinformatics 2009 Jun 1; 25(11) 1422-3 https://doi.org/10.1093/bioinformatics/btp163 pmid:19304878

<a id="ref14"></a>
[14] PyMOL, The PyMOL Molecular Graphics System, Version 3.0 Schrödinger, LLC.

<a id="ref15"></a>
[15] Antonina Andreeva, Dave Howorth, Cyrus Chothia, Eugene Kulesha, Alexey Murzin, SCOP2 prototype: a new approach to protein structure mining. (2014) Nucl. Acid Res., 42 (D1): D310-D314 and Antonina Andreeva, Eugene Kulesha, Julian Gough, Alexey Murzin, The SCOP database in 2020: expanded classification of representative family and superfamily domains of known protein structures. (2020) Nucl. Acid Res., 48 (D1): D376-D382

<a id="ref16"></a>
[16] De Cao, Nicola & Kipf, Thomas. (2018). MolGAN: An implicit generative model for small molecular graphs. 10.48550/arXiv.1805.11973. 

<a id="ref17"></a>
[17] Basu, V. (2024, December 2017). Drug molecule generation with VAE. Keras. https://mng.bz/rKve

<div style="page-break-after: always;"></div>

## List of Figures

<div style="page-break-after: always;"></div>


## List of Tables

<div style="page-break-after: always;"></div>


## Appendix
