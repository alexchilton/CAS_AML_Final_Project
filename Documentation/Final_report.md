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

All files with the `.pdb` extension are recursively identified. Each file is validated by size $( S )$, and only those satisfying the condition:

$$
(S \leq S_{\text{max}} \quad \text{with} \quad S_{\\text{max}} = 25\,\text{MB}
)$$
are processed.


#### Structural Parsing

Using BioPython, the structure is parsed at the atomic level. Each standard residue is extracted and stored by chain. A tuple representation is maintained:

$$(
a_i = (\text{chain\_id}, r_i, t_i, a_i^n, e_i, \vec{x}_i)
)$$

where:
- $( r_i )$ is the residue number,
- $( t_i )$ is the residue type (three-letter code),
- $( a_i^n )$ is the atom name,
- $( e_i )$ is the element,
- $( \vec{x}_i \in \mathbb{R}^3 )$ is the atomic position.


#### Secondary Structure Assignment

Residue-level secondary structure is computed using DSSP, resulting in a mapping:

$$(
(r_i, \text{chain\_id}) \mapsto s_i \in \{H, E, C, G, I, T, S, ?\}
)$$

where $( s_i )$ denotes the DSSP-assigned secondary structure class.


#### Graph Construction

A residue-level protein graph $( G = (V, E) )$ is built where:
- $( V )$ are nodes representing residues,
- $( E )$ are edges representing either peptide bonds or spatial proximity.

Edges $( (i, j) \in \mathbb{E} )$ are defined by:

$$(
\|\vec{x}_i^{CA} - \vec{x}_j^{CA}\|_2 < \delta \quad \text{with} \quad \delta = 7.0\,\text{\AA}
)$$


#### Geometric Feature Extraction

Each chain is processed independently. For each residue:
- **Bond Lengths**:
  $$(
  d_{ij} = \|\vec{x}_i - \vec{x}_j\|_2
  )$$

- **Bond Angles** (for triplet  $i - j - k$):
  $$(
  \theta_{ijk} = \cos^{-1}\left( \frac{ (\vec{x}_i - \vec{x}_j) \cdot (\vec{x}_k - \vec{x}_j) }{ \|\vec{x}_i - \vec{x}_j\|_2 \cdot \|\vec{x}_k - \vec{x}_j\|_2 } \right)
  )$$

- **Torsion Angles** (for quadruplet $i - j - k - l$):

  Let $\vec{b}_1 = \vec{x}_j - \vec{x}_i,  \vec{b}_2 = \vec{x}_k - \vec{x}_j,  \vec{b}_3 = \vec{x}_l - \vec{x}_k$

  $$
  \phi = \text{atan2}\left( \frac{ \vec{b}_2 \cdot (\vec{b}_1 \times \vec{b}_3) }{ \|\vec{b}_1 \times \vec{b}_2\|_2 \cdot \|\vec{b}_2 \times \vec{b}_3\|_2 }, (\vec{b}_1 \times \vec{b}_2) \cdot (\vec{b}_2 \times \vec{b}_3) \right)
  $$

#### Physicochemical Features

Charges are heuristically assigned:

$$q_i = 
\begin{cases}
-0.5 & \text{if } e_i = \text{N} \\
0.5 & \text{if } e_i = \text{C} \\
-0.5 & \text{if } e_i = \text{O} \\
0 & \text{otherwise}
\end{cases}$$

Hydrophobic residues are identified using a fixed set $H = \{\text{ALA}, \text{VAL}, \text{LEU}, \ldots\}$, such that:

$$\text{hydrophobic}(r_i) = \begin{cases}
1 & \text{if } t_i \in H \\
0 & \text{otherwise}
\end{cases}$$


#### Output

For each PDB entry, the following artifacts are saved:
- Serialized graph object with node/edge attributes.
- Pickled feature dictionaries.
- JSON summary reporting file size, number of residues, edges, and status.