
library(foreach)
ARGV = commandArgs(trailingOnly = TRUE)
data_folder = ARGV[1]
P = as.numeric(ARGV[2])
prop = as.numeric(ARGV[3])
data_folder = paste0(data_folder, "/P_", P, "/prop_", prop)
n.data_set = 30


if(!dir.exists(paste0(data_folder, "/results"))){
  dir.create(paste0(data_folder, "/results"))
}

if(!dir.exists(paste0(data_folder, "/summary"))){
  dir.create(paste0(data_folder, "/summary"))
}
############################
#### Spartaco
#############################
if(!dir.exists(paste0(data_folder, "/spartaco"))){
  dir.create(paste0(data_folder, "/spartaco"))
}


mydata5 <- foreach(i=1:n.data_set, .combine = "rbind") %do% 
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
    ##### spartaco
    set.seed(1)
    log_count = log(count_mat / sg + 1)
    coordinates = as.matrix(cbind(sce$row, sce$col))
    library(spartaco)
    res_spartaco = spartaco(log_count, coordinates = coordinates, K = L_true, R = K_true, 
                            max.iter = 40, verbose = F)
    save(res_spartaco, file = paste0(data_folder, "/spartaco/", i, ".RData"))
    ARI_col = randIndex(table(res_spartaco$Ds, cluster_true))
    ARI_row = randIndex(table(res_spartaco$Cs, rho_true))
    
    c(ARI_col, ARI_row)
    print(c(i, ARI_col, ARI_row))
  }
write.table(mydata5, file = paste0(data_folder, "/summary/", "spartaco", ".txt"))
