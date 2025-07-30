## TiGeRN : Teams in Gene Regulatory Network 

This repository provides an implementation of a spectral clustering pipeline tailored for signed graphs, specifically gene regulatory networks (GRNs) with activation and repression edges. The algorithm computes a truncated affinity matrix, forms a normalized Laplacian, solves the eigenvalue problem, and applies k-means clustering on the eigenvectors to partition the graph into coherent modules.

---

### Table of Contents

1. [Installation](#installation)
2. [Usage](#usage)
3. [Algorithm Overview](#algorithm-overview)
4. [Function Reference](#function-reference)

   * [Affinity Computation](#affinity-computation)
   * [Degree Matrix](#degree-matrix)
   * [Normalized Laplacian](#normalized-laplacian)
   * [Eigen Decomposition](#eigen-decomposition)
   * [Clustering and Evaluation](#clustering-and-evaluation)
   * [Graph Drawing](#graph-drawing)
5. [Full Pipeline](#full-pipeline)
6. [Example](#example)
7. [Dependencies](#dependencies)
8. [License](#license)

---

## Installation

1. Clone the repository:

   ```bash
   git clone https://github.com/yourusername/spectral-grn-clustering.git
   cd spectral-grn-clustering
   ```

2. Install dependencies:

   ```bash
   pip install -r requirements.txt
   ```

   Requirements include:

   * `numpy`
   * `scikit-learn`
   * `networkx`
   * `matplotlib`

---

## Usage

Import the `full_pipeline_signed_laplacian` function and run on your adjacency matrix `A` and gene labels:

```python
from spectral_clustering import full_pipeline_signed_laplacian

# A: signed adjacency matrix (activation=1, repression=-1)
# labels: list of gene names corresponding to A's rows/cols
A_aff, A_sign, D, L, eig_data, clusters = full_pipeline_signed_laplacian(
    A,
    n=4,
    weights=[1,0.5,0.25,0.125],
    affinity_mode=True,
    k_min=2,
    k_max=10,
    labels=gene_list,
    plot=True,
    sd=42
)
```

This will compute and print clustering results, and optionally display the original and clustered GRN side by side.

---

## Algorithm Overview

The pipeline proceeds through the following steps:

1. **Affinity Computation**: Builds a truncated matrix exponential or weighted sum of powers of the signed adjacency to capture multiscale connectivity.
2. **Sign Extraction**: Converts affinities to signed interactions and re-applies the power-sum to reinforce sign structure.
3. **Degree Matrix Construction**: Computes a normalized degree matrix that handles nodes with zero degree gracefully.
4. **Normalized Laplacian**: Constructs \$L = I - D A D\$ for spectral embedding.
5. **Eigen Decomposition**: Solves for eigenvalues and eigenvectors of \$L\$.
6. **Optimal Cluster Count**: Evaluates silhouette scores across candidate \$k\$ to select the ideal number of clusters.
7. **Spectral Embedding & Clustering**: Runs k-means on the top \$k\$ eigenvectors to assign each node to a cluster.
8. **Visualization** *(optional)*: Draws the input and clustered networks with distinct node coloring and edge styling.

---

## Function Reference

### Affinity Computation

* **`compute_affinity(A, alpha=1.0, Nmax=10)`**

  * Approximates the matrix exponential via a truncated series.
  * Returns absolute affinity `K_abs` and raw `K`.

* **`calc_m(A, n=4, weights=[1,0.5,0.25,0.125])`**

  * Computes weighted sum of powers: \$\sum\_{i=1}^n w\_i A^i\$.
  * Returns absolute `M_abs` and raw `M`.

### Degree Matrix

* **`calc_degree_matrix(A)`**

  * Calculates \$D\$ with entries \$1/\sqrt{\sum\_j|A\_{ij}|}\$.
  * Replaces infinities (from zero-degree) with zero.

### Normalized Laplacian

* **`calc_laplacian(A, D)`**

  * Constructs \$L = I - D A D\$.

### Eigen Decomposition

* **`solve_eigen_problem(L)`**

  * Computes eigenvalues and eigenvectors of the symmetric matrix \$L\$.

### Clustering and Evaluation

* **`cluster_embeddings(U, k)`**: Applies k-means on rows of embedding matrix `U`.
* **`find_optimal_k(U, k_min, k_max)`**: Searches \$k\$ by silhouette score, returning the best \$k\$ and score.
* **`cluster_to_dict(labels, k)`**: Converts cluster label array to dictionary mapping clusters to node indices.

### Graph Drawing

* **`make_graph(adj_matrix, genes)`**: Builds a directed Graph with `activation`/`repression` edge attributes.
* **`draw(G, seed, node_colors, ax)`**: Visualizes nodes and styled edges with NetworkX and Matplotlib.

---

## Full Pipeline Function

```python
def full_pipeline_signed_laplacian(
    A,
    n=4,
    weights=[1,0.5,0.25,0.125],
    affinity_mode=True,
    k_min=2,
    k_max=10,
    labels=None,
    plot=False,
    sd=42
):
    """
    Executes the complete spectral clustering pipeline:
     1. Affinity computation (weighted multi-step powers)
     2. Sign extraction and re-affinity
     3. Degree matrix and normalized Laplacian
     4. Spectral embedding (eigenvectors)
     5. Optimal k selection via silhouette score
     6. k-means clustering on eigenvectors
     7. Optional plotting of original vs clustered networks
    Returns intermediate and final results.
    """
```

---

## Example

```python
import numpy as np
from spectral_clustering import full_pipeline_signed_laplacian

# Example adjacency matrix for 6 genes
A = np.array([
    [ 0,  1,  0, -1,  0,  0],
    [ 1,  0,  1,  0, -1,  0],
    [ 0,  1,  0,  0,  0, -1],
    [-1, 0,  0,  0,  1,  0],
    [ 0, -1, 0,  1,  0,  1],
    [ 0,  0, -1, 0,  1,  0]
])
labels = ["G1","G2","G3","G4","G5","G6"]

results = full_pipeline_signed_laplacian(
    A, labels=labels, plot=True
)
```

---

## Dependencies

```text
numpy
scikit-learn
networkx
matplotlib
```

---

## License

Distributed under the MIT License.
