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

10-20 lines

Lorem ipsum dolor sit amet, consectetuer adipiscing elit, sed diam nonummy nibh euismod tincidunt ut laoreet dolore magna aliquam erat volutpat. Ut wisi enim ad minim veniam, quis nostrud exerci tation ullamcorper.

<div style="page-break-after: always;"></div>

## Table of contents  

- [Abstract](#abstract)
- [Table of contents](#table-of-contents)
- [Project objectives](#project-objectives)
- [Methods](#methods)
- [Data](#data)
  - [Protein Structure Preprocessing](#protein-structure-preprocessing)
    - [Overview](#overview)
    - [Input Filtering](#input-filtering)
    - [Structural Parsing](#structural-parsing)
    - [Secondary Structure Assignment](#secondary-structure-assignment)
    - [Graph Construction](#graph-construction)
    - [Geometric Feature Extraction](#geometric-feature-extraction)
    - [Physicochemical Features](#physicochemical-features)
    - [Output](#output)

<div style="page-break-after: always;"></div>

## Project objectives

<div style="text-align: justify;">

Proteins are the most abundant biological macromolecules, present in all cells and cellular compartments. They exhibit substantial diversity in both size and function, ranging from small peptides to large polymeric complexes. Functionally, proteins are central to virtually all biological processes and represent the primary end products of genetic information pathways. 

Structurally, proteins are linear polymers of amino acids linked by peptide bonds. Twenty standard α-amino acids constitute protein sequences, each containing a central (α) carbon bonded to an amino group, a carboxyl group, a hydrogen atom, and a variable side chain (R group). These side chains differ in size, structure, polarity, and charge, influencing the solubility and conformational properties of the protein. Except for glycine, all amino acids have four distinct substituents on the α-carbon, imparting chirality to the molecule. The repeating sequence of α-carbons forms the protein backbone.  

Protein structure is conventionally described at four hierarchical levels. The primary structure is the linear sequence of amino acids. The secondary structure refers to local conformations such as α-helices and β-sheets. The tertiary structure is the full three-dimensional configuration of a single polypeptide chain, while the quaternary structure represents the assembly of multiple folded subunits into a functional complex. Unlike the primary structure, which is genetically encoded, the secondary and tertiary structures are influenced by both the amino acid sequence and environmental factors. Predicting these higher-order structures remains a complex task due to the intricate balance of biochemical constraints and conformational variability.  

In this work, we aim to investigate and compare computational methodologies for generating plausible three-dimensional protein structures from PDB-derived representations. Focusing on the generative potential of machine learning techniques, we explore a spectrum of approaches including gradient-based optimization, graph-based representations, and deep generative models. The objective is to characterize their respective modeling capabilities, limitations, and potential integration strategies for future protein structure prediction frameworks.

</div>

## Methods 

<div style="text-align: justify;">

0.5-1.0 page

Which infrastructure, tools, software libraries, statistical methods etc do you intend to use. It is clear that you may not know all this at this stage, but try to make yourself some thoughts, even if it is going to change during the CAS.

</div>


## Data

<div style="text-align: justify;">

### Protein Structure Preprocessing

#### Overview

We propose a memory-efficient preprocessing pipeline for the systematic extraction of structural, geometric, and physicochemical features from Protein Data Bank (PDB) files. The workflow is implemented as a modular system and designed to handle large-scale datasets with minimal memory overhead.

#### Input Filtering

All files with the `.pdb` extension are recursively identified. Each file is validated by size $S$, and only those satisfying the condition:

$$S \leq S_{\text{max}} \quad \text{with} \quad S_{\\text{max}} = 25\,\text{MB}$$
are processed.


#### Structural Parsing

Using BioPython, the structure is parsed at the atomic level. Each standard residue is extracted and stored by chain. A tuple representation is maintained:

$$
a_i = (\text{chain\_id}, r_i, t_i, a_i^n, e_i, \vec{x}_i)
$$

where:
- $r_i$ is the residue number,
- $t_i$ is the residue type (three-letter code),
- $a_i^n$ is the atom name,
- $e_i$ is the element,
- $\vec{x}_i \in \mathbb{R}^3$ is the atomic position.


#### Secondary Structure Assignment

Residue-level secondary structure is computed using DSSP, resulting in a mapping:

$$
(r_i, \text{chain\_id}) \mapsto s_i \in \{H, E, C, G, I, T, S, ?\}
$$

where $s_i$ denotes the DSSP-assigned secondary structure class.


#### Graph Construction

A residue-level protein graph $G = (V, E)$ is built where:
- $V$ are nodes representing residues,
- $E$ are edges representing either peptide bonds or spatial proximity.

Edges $(i, j) \in \mathbb{E}$ are defined by:

$$
\|\vec{x}_i^{CA} - \vec{x}_j^{CA}\|_2 < \delta \quad \text{with} \quad \delta = 7.0\,\text{\AA}
$$


#### Geometric Feature Extraction

To characterize the 3D conformation of protein chains, we compute several geometric features: bond lengths, bond angles, and torsion angles. These features capture local spatial relationships between atoms and are fundamental for understanding protein structure and stability.

Due to the computational cost associated with calculating all geometric interactions at the atomic level—and in line with the specific objectives of the present study—we restrict these calculations to the backbone atoms only. These include the nitrogen (N), alpha carbon (CA), carbon (C), and oxygen (O) atoms that form the peptide backbone of each amino acid residue. This reduction significantly improves computational efficiency while retaining essential structural information relevant for most downstream analyses. The underlying implementation, however, has been designed to support full-atom calculations and can be extended in future work.

Each chain is processed independently. For each residue: 

**Bond Lengths**: Bond lengths are defined as the Euclidean distance between two bonded atoms. In this context, they are computed between pairs of backbone atoms within the same residue or between sequential residues. The bond length $d_{ij}$ between atoms $i$ and $j$ is given by: 

  $$
  d_{ij} = \|\vec{x}_i - \vec{x}_j\|_2
  $$

  where $\vec{x}_i$ and $\vec{x}_j$ are the 3D coordinates of the respective atoms.


**Bond Angles** (for triplet  $i - j - k$): Bond angles provide a measure of the angle formed by three consecutive atoms and are useful for assessing local bending in the protein chain. Given a triplet of atoms $(i, j, k)$ , the bond angle $\theta_{ijk}$  is computed as:
  $$
  \theta_{ijk} = \cos^{-1}\left( \frac{ (\vec{x}_i - \vec{x}_j) \cdot (\vec{x}_k - \vec{x}_j) }{ \|\vec{x}_i - \vec{x}_j\|_2 \cdot \|\vec{x}_k - \vec{x}_j\|_2 } \right)
  $$

  Only triplets involving backbone atoms are considered, including intra-residue angles (e.g., N–CA–C) and inter-residue connections (e.g., C–N–CA across peptide bonds).

**Torsion Angles** : also called dihedral angles, they describe the rotational relationship between four sequential atoms and are critical for capturing the 3D folding of the protein. These angles are computed exclusively on the backbone atoms and include the three canonical torsions: phi $(\phi)$, psi $(\psi)$ and omega $(\omega)$. 
Each angle corresponds to a specific actom configuration across sequential residues:  

  - $\phi$: calculated from atoms $C(i-1 - N(i) - CA(i) - C(i)$, measuring rotation on the $N-CA$ bond
  - $\psi$: from $N(i)-CA(i)-C(i)–N(i+1)$, measuring the rotation on the $CA–C$ bond
  - $\omega$: from $CA(i–1)–C(i–1)–N(i)–CA(i)$, reflecting th ecis/trans conformation across the peptide bond.    

For atoms $(i, j, k, l)$ , the torsion angle $\phi$ is calculated using the angle between the planes formed by $(i, j, k)$ and $(j, k, l)$:

  Let $\vec{b}_1 = \vec{x}_j - \vec{x}_i,  \vec{b}_2 = \vec{x}_k - \vec{x}_j,  \vec{b}_3 = \vec{x}_l - \vec{x}_k$

  $$
  \phi = \text{atan2}\left( \frac{ \vec{b}_2 \cdot (\vec{b}_1 \times \vec{b}_3) }{ \|\vec{b}_1 \times \vec{b}_2\|_2 \cdot \|\vec{b}_2 \times \vec{b}_3\|_2 }, (\vec{b}_1 \times \vec{b}_2) \cdot (\vec{b}_2 \times \vec{b}_3) \right)
  $$

  where $\vec{x}_i$, $\vec{x}_j$, $\vec{x}_k$, $\vec{x}_l$ are the 3D coordinates of the respective atoms. 

  This formulation is applied uniformly across all torsion types. While current calculations are limited to backbone atoms, the implementation is generalizable and can be extended to side-chain or full-atom torsions in future work.

#### Physicochemical Features

Charges are heuristically assigned:

$$
q_i = 
\begin{cases}
-0.5 & \text{if } e_i = \text{N} \\
0.5 & \text{if } e_i = \text{C} \\
-0.5 & \text{if } e_i = \text{O} \\
0 & \text{otherwise}
\end{cases}
$$

Hydrophobic residues are identified using a fixed set $H = \{\text{ALA}, \text{VAL}, \text{LEU}, \ldots\}$, such that:

$$
\text{hydrophobic}(r_i) = \begin{cases}
1 & \text{if } t_i \in H \\
0 & \text{otherwise}
\end{cases}
$$


#### Output

For each PDB entry, the following artifacts are saved:
- Serialized graph object with node/edge attributes.
- Pickled feature dictionaries.
- JSON summary reporting file size, number of residues, edges, and status.