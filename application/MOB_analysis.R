rm(list=ls())
load("application/data/MOB_sce_filter.RData")


n_gene = 1000
K = 4
L = 3


library(scran)
sce <- scater::logNormCounts(sce1)
dec <- scran::modelGeneVar(sce, assay.type = "logcounts")
top <- scran::getTopHVGs(dec, n = n_gene)

sce = sce[top, ]



source("R/utils.R")
Adj = find_neighbors(sce, "ST", "lattice")
neighbors = find_neighbor_index(Adj, "ST")


sizefactor = sce$sizefactor
genefactor = rowData(sce)$genefactor


count_mat_sce = assay(sce, "counts")
P = dim(count_mat_sce)[1]
N = dim(count_mat_sce)[2]


s_mat = matrix(rep(sizefactor, each = P), P, N)
g_mat = matrix(rep(genefactor, N), P, N)
sg = s_mat * g_mat



##### initialization
set.seed(seed)
L_init = L
K_init = K
mu_t = matrix(abs(rnorm(L_init*K_init)), L_init, K_init)
mu0_t = runif(P, 1, 2)
rho_t = c(sample(0:L_init, size = P-L_init-1, replace = TRUE), 0:L_init)

group_t = c(sample(K_init, size = N-K_init, replace = TRUE), 1:K_init)
pi_t = 0.1


library(Rcpp)
library(RcppArmadillo)
library(SingleCellExperiment)
sourceCpp("R/co_clustering_mcmc.cpp")
n_iter = 10000

begin = proc.time()
result = runMCMC(group_t-1 , rho_t, mu_t, mu0_t, pi_t, count_mat_sce, sg, 
                 alpha0 = 1, beta0 = 1, a0 = 1, b0 = 1, f = 1, G = neighbors, 
                 BETA = 1, GAMMA = 1, max_iters = n_iter)
end = proc.time()
time = end - begin

library(flexclust)
ARIvalue = randIndex(table(result$group_iter[sce$label != "unknown", n_iter] +1, sce$label[sce$label != "unknown"]))
res_ICL = ICL(K = K_init, L = L_init, data = count_mat_sce, sg = sg, result = result, iter = n_iter)

save(result, file = paste0("MOB/chain/", "seed_", seed, ".RData"))

out = c(seed, ARIvalue, res_ICL, time)
names(out) = c("seed", "ARI", "ICL", "time")
print(out)
