# Title TC

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
