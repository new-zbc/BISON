library(SingleCellExperiment)
library(flexclust)
ARGV = commandArgs(trailingOnly = TRUE)
data_folder = ARGV[1]
P = as.numeric(ARGV[2])
prop = as.numeric(ARGV[3])
folder = paste0(data_folder, "/P_", P, "/prop_", prop)


out_df = data.frame()

for(i in 1:30){
  if(file.exists(paste0(folder, "/results/", i, ".RData"))){
    
    data_file = paste0(folder, "/data/", i, ".RData")
    load(data_file)
    
    cluster_true = sce$label
    rho_true = rowData(sce)$rho_true
    
    load(paste0(folder, "/results/", i, ".RData"))

    est = result$rho_iter[, 200]
    est[est >= 1] = 1

    rho_true[rho_true >= 1] = 1

    cft <- table(est, rho_true)
    tp <- cft[2, 2]
    tn <- cft[1, 1]
    fp <- cft[2, 1]
    fn <- cft[1, 2]
    sensitivity <- tp / (tp + fn)
    specificity <- tn/(fp + tn)
    accuracy <- (tp + tn)/(tp + tn + fp + fn)
    out = c(sensitivity, specificity, accuracy)
    out_df = rbind(out_df, out)
  }
}
colnames(out_df) = c("Sensitivity", "Specificity", "Accuracy")
write.table(out_df, file = paste0(folder, "/summary/", "accuracy", ".txt"))