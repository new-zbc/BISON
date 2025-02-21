library(Rcpp)
library(RcppArmadillo)
library(SingleCellExperiment)
sourceCpp("R/co_clustering_mcmc.cpp")


BISON <- function(count_mat, sg, neighbors, K, R, f = 1, n_iters = 1000, seed = 1){
  
  P = dim(count_mat)[1]
  N = dim(count_mat)[2]
  
  # hyperparameter
  alpha0 = 1
  beta0 = 1
  a0 = 1
  b0 = 1
  BETA = 1
  GAMMA = 1
  
  ##### initialization
  set.seed(seed)
  L_init = L
  K_init = K
  mu_t = matrix(1, L_init, K_init)
  mu0_t = rep(1, P)
  rho_t = sample(L_init, size = P, replace = TRUE)
  rho_t[sample(P, size = 10)] = 0
  group_t = sample(K_init, size = N, replace = T)
  pi_t = 0.01
  
  result = runMCMC(group_t-1 , rho_t, mu_t, mu0_t, pi_t, count_mat, sg, 
                   alpha0 = alpha0, beta0 = beta0, a0 = a0, b0 = b0, f = f, G = neighbors, 
                   BETA = BETA, GAMMA = BETA, max_iters = n_iters)
  
  pred_spot_label = result$group_iter[, n_iters]
  pred_gene_label = result$rho_iter[, n_iters]
  output = list()
  output$pred_spot_label = pred_spot_label
  output$pred_gene_label = pred_gene_label
  output$MCMC = result
  
  return(output)
}