library(foreach)
ARGV = commandArgs(trailingOnly = TRUE)
data_folder = ARGV[1]
P = as.numeric(ARGV[2])
prop = as.numeric(ARGV[3])
data_folder = paste0(data_folder, "/P_", P, "/prop_", prop)
n.data_set = 30
n.iters = 200

############################################
### proposed method
############################################
if(!dir.exists(paste0(data_folder, "/results"))){
  dir.create(paste0(data_folder, "/results"))
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

  ##### initialization
  set.seed(1)
  L_init = 3
  K_init = K_true
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
  
  begin = proc.time()
  result = runMCMC(group_t-1 , rho_t, mu_t, mu0_t, pi_t, count_mat, sg, 
                   alpha0 = 1, beta0 = 1, a0 = 1, b0 = 1, f = 1, G = neighbors, 
                   BETA = 1, GAMMA = 1, max_iters = n.iters)
  end = proc.time()
  time = end - begin

  save(result, file = paste0(data_folder, "/results/", i, ".RData"))

  ARI_col = randIndex(table(result$group_iter[, n.iters], cluster_true))
  ARI_row = randIndex(table(result$rho_iter[, n.iters], rho_true))
  K_col = length(unique(result$group_iter[, n.iters]))
  L_row = length(unique(result$rho_iter[, n.iters]))
  out = c(ARI_col, ARI_row, K_col, L_row, time)
  print(out)
}
write.table(mydata1, file = paste0(data_folder, "/summary/", "bison", ".txt"))
