library(SingleCellExperiment)
library(flexclust)
ARGV = commandArgs(trailingOnly = TRUE)
data_folder = ARGV[1]
P = as.numeric(ARGV[2])
prop = as.numeric(ARGV[3])
folder = paste0(data_folder, "/P_", P, "/prop_", prop)


out_df = data.frame()

for(i in 1:30){
  if(file.exists(paste0(folder, "/spartaco/", i, ".RData"))){
    
    data_file = paste0(folder, "/data/", i, ".RData")
    load(data_file)
    
    cluster_true = sce$label
    rho_true = rowData(sce)$rho_true
    
    load(paste0(folder, "/spartaco/", i, ".RData"))

    ARI_col = randIndex(table(res_spartaco$Ds, cluster_true))
    ARI_row = randIndex(table(res_spartaco$Cs, rho_true))
    
    out = c(ARI_col, ARI_row)
    
    out_df = rbind(out_df, out)
  }
}

write.table(out_df, file = paste0(folder, "/summary/", "spartaco", ".txt"))