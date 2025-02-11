library(foreach)
ARGV = commandArgs(trailingOnly = TRUE)
data_folder = ARGV[1]
P = as.numeric(ARGV[2])
prop = as.numeric(ARGV[3])
data_folder = paste0(data_folder, "/P_", P, "/prop_", prop)
n.data_set = 50


if(!dir.exists(paste0(data_folder, "/results"))){
  dir.create(paste0(data_folder, "/results"))
}

if(!dir.exists(paste0(data_folder, "/summary"))){
  dir.create(paste0(data_folder, "/summary"))
}
###########################
#Kmeans
############################
mydata2 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar% 
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
    ##### Kmeans 
    set.seed(1)
    log_count = log(count_mat / sg + 1)
    res_kmean = kmeans(x=t(log_count), centers = K_true, nstart = 5)
    ARI_col = randIndex(table(res_kmean$cluster, cluster_true))
    res_kmean = kmeans(x=log_count, centers = L_true, nstart = 5)
    ARI_row = randIndex(table(res_kmean$cluster, rho_true))
    
    c(ARI_col, ARI_row)
  }
write.table(mydata2, file = paste0(data_folder, "/summary/", "Kmeans", ".txt"))


#########################
###  BC
##########################
mydata3 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar% 
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
    ##### BC
    set.seed(1)
    log_count = log(count_mat / sg + 1)
    library(sparseBC)
    res_sparse_BC = sparseBC(t(log_count), K_true, L_true, lambda = 0)
    ARI_col = randIndex(table(res_sparse_BC$Cs, cluster_true))
    ARI_row = randIndex(table(res_sparse_BC$Ds, rho_true))
    
    c(ARI_col, ARI_row)
  }
write.table(mydata3, file = paste0(data_folder, "/summary/", "BC", ".txt"))


#########################
###sparse BC
##########################
mydata4 <- foreach(i=1:n.data_set, .combine = "rbind") %dopar% 
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
    ##### sparse BC
    set.seed(1)
    log_count = log(count_mat / sg + 1)
    library(sparseBC)
    res_sparse_BC = sparseBC(t(log_count), K_true, L_true, lambda = 20)
    ARI_col = randIndex(table(res_sparse_BC$Cs, cluster_true))
    ARI_row = randIndex(table(res_sparse_BC$Ds, rho_true))
    
    c(ARI_col, ARI_row)
  }
write.table(mydata4, file = paste0(data_folder, "/summary/", "sparseBC", ".txt"))