library(foreach)
ARGV = commandArgs(trailingOnly = TRUE)
data_folder = ARGV[1]
P = as.numeric(ARGV[2])
prop = as.numeric(ARGV[3])
data_folder = paste0(data_folder, "/P_", P, "/prop_", prop)
n.data_set = 30
n.iters = 50

############################################
### proposed method
############################################
if(!dir.exists(paste0(data_folder, "/model_selection"))){
  dir.create(paste0(data_folder, "/model_selection"))
}

if(!dir.exists(paste0(data_folder, "/summary"))){
  dir.create(paste0(data_folder, "/summary"))
}

mydata1 <- foreach(i=1:n.data_set, .combine = "rbind") %do% 
{
  library(SingleCellExperiment)
  library(flexclust)
  data_file = paste0(data_folder, "/data/", i, ".RData")
  load(data_file)
  count_mat = assay(sce, "counts")
  P = dim(count_mat)[1]
  N = dim(count_mat)[2]
  cluster_true = sce$label
  rho_true = rowData(sce)$rho_true
  K_true = length(unique(cluster_true))
  L_true = length(unique(rho_true))

  store = data.frame()
  Ks = c(2, 3, 4, 5, 6)
  Ls = c(1, 2, 3, 4, 5)
  for(i in 1:length(Ks)){
    for(j in 1:length(Ls)){

  ##### initialization
  set.seed(1)
  L_init = Ls[j]
  K_init = Ks[i]

  mu_t = matrix(1, L_init, K_init)
  mu0_t = rep(1, P)
  rho_t = sample(L_init, size = P, replace = TRUE)
  rho_t[sample(P, size = 10)] = 0
  group_t = sample(K_init, size = N, replace = T)
  pi_t = 0.01
  
  
  library(Rcpp)
  library(RcppArmadillo)
  library(SingleCellExperiment)
  sourceCpp("R/co_clustering_mcmc.cpp")
  source("R/utils.R")
  
  begin = proc.time()
  result = runMCMC(group_t-1 , rho_t, mu_t, mu0_t, pi_t, count_mat, sg, 
                   alpha0 = 1, beta0 = 1, a0 = 1, b0 = 1, f = 1, G = neighbors, 
                   BETA = 1, GAMMA = 1, max_iters = n.iters)
  end = proc.time()
  time = end - begin

  res_ICL1 = ICL(K = K_init, L = L_init, data = count_mat, sg = sg, result = result, iter = n.iters)
  res_ICL2 = ICL4(K = K_init, L = L_init, data = count_mat, sg = sg, result = result, iter = n.iters)

  ARI_col = randIndex(table(result$group_iter[, n.iters], cluster_true))
  ARI_row = randIndex(table(result$rho_iter[, n.iters], rho_true))

  output = c(K_init, L_init, ARI_col, ARI_row, res_ICL1, res_ICL2)
  
  save(output, file = paste0(data_folder, "/model_selection/", i, ".RData"))
  
    }
  }

  index = which.max(output[, 6])
#   ARI_col = randIndex(table(result$group_iter[, n.iters], cluster_true))
#   ARI_row = randIndex(table(result$rho_iter[, n.iters], rho_true))
#   K_col = length(unique(result$group_iter[, n.iters]))
#   L_row = length(unique(result$rho_iter[, n.iters]))
  out = output[index, ]
  print(out)
  out
}

colnames(mydata1) = c("K_init", "L_init", "ARI_col", "ARI_row", "res_ICL1", "res_ICL2")
write.table(mydata1, file = paste0(data_folder, "/model_selection/", "model_selection", ".txt"))
