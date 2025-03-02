## BISON: Bi-clustering of spatial omics data with feature selection

### Overview

BISON is a Bayesian bi-clustering framework with feature selection for spatially resolved transcriptomics data. It detects informative genes with heterogeneity among spots and clusters spots and informative genes simultaneously. 

### Dependencies

The required R packages:

* Rcpp
* RcppArmadillo
* SingleCellExperiment
* scran
* scater

### Quick Start

The main function ``BISON()`` can be loaded via

```R
source("R/main.R")
Result <- BISON(count_mat, sg, neighbors, K, R, f=1, n_iters=1000, seed=1)
```

Here is the explanation for the available parameters:

* ``count_mat``: The $p$-by-$n$ count matrix for spatially resolved transcriptomics data
* ``sg``: The $p$-by-$n$ matrix for the production of size factor and gene factor, where the entry $[sg]_{ij} = s_{i}g_{j}$. 
* ``neighbors``: The $n$-by-$Q$ matrix with Q equal to the maximum number of neighbors for each spot. For SRT data from ST platform and $10$x Visium, $Q = 4$ and $Q = 6$, respectively. The entry in ``neighbors`` is the index of neighboring spots. If the number of neighbors is less than $Q$, the other entries are filled with $0$. 
* ``K``: The number of column (spot) clusters.
* ``R``: The number of informative feature (gene) clusters. Since there is a NonDG gene group, the total number of gene groups is $R +1$. 
* ``f``: The hyperparameter in the MRF model controlling the spatial smoothness of spot clusters. 
* ``n_iters``: The max number of MCMC iterations. 
* ``seed``: The seed to reproduce the results.  

### Tutorials

For a tutorial on a real application, please refer to [Tutorial.html](https://new-zbc.github.io/BISON/MOB_tutorial.html)

If there is any bugs without being fixed timely, please contact author directly [zhubencong@gmail.com](mailto:zhubencong@gmail.com)
