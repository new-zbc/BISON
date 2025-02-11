library(SingleCellExperiment)
library(flexclust)
source("R/utils.R")

find_coordinate <- function(height, width){
  result = matrix(0, height*width, 2)
  for(i in 1:(height*width)){
    part1 = i %/% height
    part2 = i %% height
    
    if(part2 != 0){
      result[i, 1] = part1 + 1
      result[i, 2] = part2
    }
    else{
      result[i, 1] = part1
      result[i, 2] = part2 + height
    }
  }
  return(result)
}




ARGV = commandArgs(trailingOnly = TRUE)
delta = as.numeric(ARGV[1])
P = as.numeric(ARGV[2])
prop = as.numeric(ARGV[3])

if(delta == 0.5){
  data_folder = "simulation3/patternMOB1"
  dir.create(data_folder)
}else if(delta == 1){
  data_folder = "simulation3/patternMOB2"
  dir.create(data_folder)
}else if(delta == 1.5){
  data_folder = "simulation3/patternMOB3"
  dir.create(data_folder)
}else{stop("Wrong delta")}

#P = 500 # total number of variables
P0 = round(P * prop)

if(!dir.exists(paste0(data_folder, "/P_", P))){
  dir.create(paste0(data_folder, "/P_", P))
}

if(!dir.exists(paste0(data_folder, "/P_", P, "/prop_", prop))){
  dir.create(paste0(data_folder, "/P_", P, "/prop_", prop))
}

data_folder = paste0(data_folder, "/P_", P, "/prop_", prop)


load("simulation3/mouse_olfactory_bulb_replicate_12.RData")
sizefactor_ture = rowMeans(count) / mean(rowMeans(count)) 


metadata = metadata[!is.na(metadata$Layer), ]
position = cbind(metadata$y, metadata$x)
label = metadata$Layer
cluster_true = rep(0, length(label))
cluster_true[label == "GCL"] = 1
cluster_true[label == "GL"] = 2
cluster_true[label == "ONL"] = 3
cluster_true[label == "MCL"] = 4

L = 4
K = length(unique(cluster_true))

#height = 40
#width = 40
rho_true = c(rep(0, P0), c(1:(L-1)), sample(1:(L-1), size = P - P0 - L +1, replace = TRUE))
rho_true = rho_true[sample(length(rho_true))]
#load("simulation1/3_clusters_pattern.RData")

for(data_index in 1:50){
  
  set.seed(data_index)
  
  #cluster_truth = c(labels) + 1
  #position = find_coordinate(height, width)
  
  ##### generate logcounts
  #N = height * width
  N = length(cluster_true)
  
  s_true = matrix(rep(sample(sizefactor_ture, N, replace = TRUE), each = P), P, N) 
  g_true = matrix(rep(runif(P, 0.5, 1.5), N), P, N)
  sg = s_true * g_true
  
  mu_true = matrix(0, L-1, K)
  # for(l in 1:(L-1)){
  #   for(k in 1:K){
  #     mu_true[l, k] = rgamma(1, 4, 1)
  #   }
  # }
  
  mu_true = matrix(0, L-1, K)
  for(l in 1:(L-1)){
    for(k in 1:K){
      mu_true[l, k] = 4 + delta *(k-1)  + delta * (l-1)
    }
  }
  
  
  count_mat = matrix(0, P, N)
  for(j in 1:P){
    if(rho_true[j] == 0){
      mu_j = runif(1, 2, 6) 
      count_mat[j, ] = rpois(N, sg[j, ] * mu_j)
    }else{
      for(i in 1:N){
        mu_ji = mu_true[rho_true[j], cluster_true[i]] + runif(1, -0.1, 0.1)
        count_mat[j, i] = rpois(1, sg[j, i] * mu_ji)
      }
    }
  }
  
  
  colnames(count_mat) = paste0("spot_", 1:dim(count_mat)[2])
  rownames(count_mat) = paste0("gene_", 1:dim(count_mat)[1])
  
  #print(mean(count_mat == 0))
  
  
  sce <- SingleCellExperiment(list(counts = count_mat), 
                              colData = DataFrame(row = position[, 1], 
                                                  col = position[, 2], 
                                                  label = cluster_true,
                                                  sizefactor = s_true[1, ]),
                              rowData = DataFrame(genefactor = g_true[, 1],
                                                  rho_true = rho_true))
  
  # #sce <- logNormCounts(sce)
  # sce <- scater::runPCA(sce, subset_row=1:2000, ncomponents=15, 
  #                       exprs_values="logcounts")
  
  Adj = find_neighbors(sce, "ST", "lattice")
  neighbors = find_neighbor_index(Adj, "ST")
  
  if(!dir.exists(paste0(data_folder, "/data"))){
    dir.create(paste0(data_folder, "/data"))
  }
  file_name = paste0(data_folder, "/", "data/", data_index, ".RData")
  save(sce, Adj, neighbors, sg,  file = file_name)
  
}

print(paste0(data_folder, ": Done"))
