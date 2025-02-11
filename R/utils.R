library(Rcpp)
library(RcppArmadillo)
library(SingleCellExperiment)


# Summarize the spatial information
# sce: SingleCellExperiment with colData including row and col
# platform: "ST" or "Visium"
# coordinate : "lattice" or "image"
find_neighbors <- function(sce, platform="ST", coordinate="lattice"){
  library(purrr)
  library(Matrix)
  if (platform == "Visium") {
    ## Spots to left and right, two above, two below
    offsets <- data.frame(x.offset=c(-2, 2, -1,  1, -1, 1),
                          y.offset=c( 0, 0, -1, -1,  1, 1))
    cutoff = sqrt(3)
  } else if (platform == "ST") {
    ## L1 radius of 1 (spots above, right, below, and left)
    offsets <- data.frame(x.offset=c( 0, 1, 0, -1),
                          y.offset=c(-1, 0, 1,  0))
    cutoff = sqrt(2)
  } else {
    stop("find_neighbors: Unsupported platform \"", platform, "\".")
  }
  
  n <- ncol(sce)
  D <- sparseMatrix(i = 1:n, j = 1:n, x = 0)
  
  if(coordinate == "lattice"){
    ## Get array coordinates (and label by index of spot in SCE)
    spot.positions <- colData(sce)[, c("col", "row")]
    spot.positions$spot.idx <- seq_len(nrow(spot.positions))
    
    ## Compute coordinates of each possible spot neighbor
    neighbor.positions <- merge(spot.positions, offsets)
    neighbor.positions$x.pos <- neighbor.positions$col + neighbor.positions$x.offset
    neighbor.positions$y.pos <- neighbor.positions$row + neighbor.positions$y.offset
    
    ## Select spots that exist at neighbor coordinates
    neighbors <- merge(as.data.frame(neighbor.positions), 
                       as.data.frame(spot.positions), 
                       by.x=c("x.pos", "y.pos"), by.y=c("col", "row"),
                       suffixes=c(".primary", ".neighbor"))
    
    
    neighbors <- neighbors[order(neighbors$spot.idx.primary, 
                                 neighbors$spot.idx.neighbor), ]
    #df_j <- split(neighbors$spot.idx.neighbor, neighbors$spot.idx.primary)
    #df_j <- unname(df_j)
    
    
    for (i in 1:n) {
      index = which(neighbors$spot.idx.primary == i)
      if(length(index) != 0)
        D[i, neighbors$spot.idx.neighbor[index]] <- 1
    }
  }else if(coordinate == "image"){
    
    spot.positions <- as.matrix(colData(sce)[, c("col", "row")])
    
    temp1 = apply(spot.positions, 1, function(rvec) crossprod(rvec))
    
    
    temp = matrix(rep(temp1, n), nrow = n) + t(matrix(rep(temp1, n), nrow = n)) - 2* tcrossprod(spot.positions, spot.positions)
    
    temp = sqrt(temp)
    
    threshold = min(temp[temp > 0])*(cutoff + 1)*0.5 
    
    for(i in 1:n){
      index = which(temp[i, ] > 0 & temp[i, ] < threshold)
      if(length(index) > 0){
        D[i, index] = 1
      }
      
    }
  }else{
    stop("wrong coordinate.")
  }
  
  return(D)
}


####### Transfer the sparse adjacent matrix to list
# Adj a sparse matrix represent the adjacent matrix
# platform: "ST" or "Visium"
find_neighbor_index <- function(Adj, platform = "Visium"){
  if (platform == "Visium"){
    num_neighbor = 6
  }else if(platform == "ST"){
    num_neighbor = 4
  }else{
    stop("wrong platform")
  }
  
  N = dim(Adj)[1]
  result = matrix(0, N, num_neighbor)
  for(i in 1:N){
    index = which(Adj[i, ] == 1)
    K = length(index)
    if(K > 0 && K <= num_neighbor){
      for(k in 1:K){result[i, k] = index[k]}
    }
  }
  return(result)
}




## Get neighbor information for other platforms
## Uses Voronoi tessellation to find nearest neighbours (share a boundary line) and their respective distances.
makeixy = function(data, formula, scale){
  m = model.frame(formula, data=data)
  if(ncol(m)!=3){
    stop("incorrect adjacency formula: id~x+y needed")
  }
  names(m)=c("id","x","y")
  m[,2]=m[,2]/scale
  m[,3]=m[,3]/scale
  m
}

voronoi_adjacency = function(data, formula, scale=1, PLOT=FALSE){
  data = makeixy(data, formula, scale)
  
  P=dim(data)[1];  # number of rows
  
  dd = deldir::deldir(data$x,data$y,suppressMsge=TRUE,plotit=PLOT);  # find adjacencies
  
  ## create adjacency matrix
  A=matrix(0,P,P);
  A[as.matrix(dd$delsgs[,c("ind1","ind2")])] = 1;
  A[as.matrix(dd$delsgs[,c("ind2","ind1")])] = 1;
  
  ## create distance matrix
  D=matrix(NA,P,P);
  D[as.matrix(dd$delsgs[,c("ind1","ind2")])] = sqrt((dd$delsgs[,c("x1")]-dd$delsgs[,c("x2")])^2+(dd$delsgs[,c("y1")]-dd$delsgs[,c("y2")])^2);
  D[as.matrix(dd$delsgs[,c("ind2","ind1")])] = sqrt((dd$delsgs[,c("x1")]-dd$delsgs[,c("x2")])^2+(dd$delsgs[,c("y1")]-dd$delsgs[,c("y2")])^2);
  
  ## create data frame of results
  N=matrix(colSums(A),P,1); # number of adjacencies for each xy$id
  
  ## create neighbor matrix
  n_neighbor <- max(N)
  nei <- matrix(0, nrow = P, ncol = n_neighbor)
  
  for (i in 1:P){
    k <- 1
    for (j in 1:P){
      if (A[i, j] == 1){
        nei[i, k] <- j
        k <- k + 1
      }
    }
  }
  
  return(list(tessellation = dd, G=A,  P=nei, Distances=D, NumNeighbours=N, ids=data[,"id"], coords=data[,c("x","y")]));
}




## preprocess of spatial data
spatialPreprocess <- function(sce, n.PCs=15, n.HVGs=2000, skip.PCA=FALSE,
                              log.normalize=TRUE, assay.type="logcounts",
                              BSPARAM=ExactParam()) {
  library(scran)
  library(scater)
  
  ## Run PCA on HVGs, log-normalizing if necessary
  if (!skip.PCA) {
    if (log.normalize){sce <- logNormCounts(sce)}
    
    dec <- modelGeneVar(sce, assay.type=assay.type)
    top <- getTopHVGs(dec, n=n.HVGs)
    sce <- scater::runPCA(sce, subset_row=top, ncomponents=n.PCs, 
                  exprs_values=assay.type)
    rowData(sce)[["is.HVG"]] <- (rownames(sce) %in% top)
  }
  
  sce
}


### sample-wise and feature-wise quality control
spatialFilter <- function(sce, cutoff_sample = 100, cutoff_feature = 0.1, cutoff_max = 10){
  library(SingleCellExperiment)
  ## Sample-wise quality control
  count = t(assay(sce, "counts"))
  index1 <- which(rowSums(count) >= cutoff_sample)
  
  
  ## Feature-wise quality control
  index2 <- which(colSums(count != 0) >= dim(count)[1]*cutoff_feature & apply(count, 2, max) >= cutoff_max)
  sce1 <- sce[index2, index1]
  
  sce1
}




ICL <- function(K, L, data, sg, result, iter){
  N = dim(data)[2]
  p = dim(data)[1]
  
  rho_hat = result$rho_iter[, iter]
  group_hat = result$group_iter[, iter]
  mu_hat = result$mu_iter[, , iter]
  mu0_hat = result$mu0_iter[, iter]
  p_hat = sum(rho_hat != 0)
  
  likelihood1 = 0
  likelihood2 = 0
  for(i in 1:N){
    for(j in 1:p){
      if(rho_hat[j]  != 0){
        likelihood1 = likelihood1 + dpois(data[j, i], lambda = sg[j, i]*mu_hat[rho_hat[j], group_hat[i]+1], log = T)
      }else{
        likelihood2 = likelihood2 + dpois(data[j, i], lambda = sg[j, i]*mu0_hat[j], log = T)
      }
    }
  }
  term1 = (-2*likelihood1 + N*log(K) + p_hat*log(L) + (K * L) *log(N*p_hat))
  term2 = (-2*likelihood2 + (p-p_hat)*log(N))
  #likelihood = (likelihood1 / p_hat + likelihood2 / (p - p_hat))*(p)
  #likelihood = likelihood1 + likelihood2
  #result = -2 *likelihood + N*log(K) + p_hat*log(L) + (K * L) *log(N*p_hat) + (p - p_hat)*log(N)
  #term1 = term1/(p_hat)
  #term2 = term2 / (p - p_hat)
  result = (term1 + term2)/2
  return(result)
}



########### scaled data
ICL2 <- function(K, L, data, sg, result, iter){
  N = dim(data)[2]
  p = dim(data)[1]
  
  rho_hat = result$rho_iter[, iter]
  group_hat = result$group_iter[, iter]
  mu_hat = result$mu_iter[, , iter]
  mu0_hat = result$mu0_iter[, iter]
  p_hat = sum(rho_hat != 0)
  
  est_prob_group = table(group_hat)
  posterior_z = sum(log(est_prob_group/N) * est_prob_group)
  
  est_prob_rho = table(rho_hat)[1:L]
  posterior_rho = sum(log(est_prob_rho/p_hat) * est_prob_rho)
  
  
  likelihood1 = 0
  likelihood2 = 0
  for(i in 1:N){
    for(j in 1:p){
      if(rho_hat[j]  != 0){
        likelihood1 = likelihood1 + dpois(data[j, i], lambda = sg[j, i]*mu_hat[rho_hat[j], group_hat[i]+1], log = T)
      }else{
        likelihood2 = likelihood2 + dpois(data[j, i], lambda = sg[j, i]*mu0_hat[j], log = T)
      }
    }
  }
  term1 = (-2*(likelihood1+posterior_rho+posterior_z) + K*log(N) + L*log(p_hat) + (K * L) *log(N*p_hat))
  term2 = (-2*likelihood2 + (p-p_hat)*log(N))
  #likelihood = (likelihood1 / p_hat + likelihood2 / (p - p_hat))*(p)
  #likelihood = likelihood1 + likelihood2
  #result = -2 *likelihood + N*log(K) + p_hat*log(L) + (K * L) *log(N*p_hat) + (p - p_hat)*log(N)
  term1 = term1/(p_hat)
  term2 = term2 / (p - p_hat)
  result = (term1 + term2 )/2
  return(result)
}



################ ICL z and rho
ICL3 <- function(K, L, data, sg, result, iter){
  N = dim(data)[2]
  p = dim(data)[1]
  
  rho_hat = result$rho_iter[, iter]
  group_hat = result$group_iter[, iter]
  mu_hat = result$mu_iter[, , iter]
  mu0_hat = result$mu0_iter[, iter]
  p_hat = sum(rho_hat != 0)
  
  est_prob_group = table(group_hat)
  posterior_z = sum(log(est_prob_group/N) * est_prob_group)
  
  est_prob_rho = table(rho_hat)[1:L]
  posterior_rho = sum(log(est_prob_rho/p_hat) * est_prob_rho)
  
  likelihood1 = 0
  likelihood2 = 0
  for(i in 1:N){
    for(j in 1:p){
      if(rho_hat[j]  != 0){
        likelihood1 = likelihood1 + dpois(data[j, i], lambda = sg[j, i]*mu_hat[rho_hat[j], group_hat[i]+1], log = T)
      }else{
        likelihood2 = likelihood2 + dpois(data[j, i], lambda = sg[j, i]*mu0_hat[j], log = T)
      }
    }
  }
  term1 = (-2*(likelihood1+posterior_rho+posterior_z) + K*log(N) + L*log(p_hat) + (K * L) *log(N*p_hat))
  term2 = (-2*likelihood2 + (p-p_hat)*log(N))
  #likelihood = (likelihood1 / p_hat + likelihood2 / (p - p_hat))*(p)
  #likelihood = likelihood1 + likelihood2
  #result = -2 *likelihood + N*log(K) + p_hat*log(L) + (K * L) *log(N*p_hat) + (p - p_hat)*log(N)
  result = (term1 + term2 )/2
  return(result)
}




################ ICL z and rho
ICL4 <- function(K, L, data, sg, result, iter, weight=1, c = 1){
  N = dim(data)[2]
  p = dim(data)[1]
  
  rho_hat = result$rho_iter[, iter]
  group_hat = result$group_iter[, iter]
  mu_hat = result$mu_iter[, , iter]
  mu0_hat = result$mu0_iter[, iter]
  p_hat = sum(rho_hat != 0)
  
  est_prob_group = table(group_hat)
  posterior_z = sum(log(est_prob_group/N) * est_prob_group)
  
  est_prob_rho = table(rho_hat)[1:L]
  posterior_rho = sum(log(est_prob_rho/p_hat) * est_prob_rho)
  
  likelihood1 = 0
  likelihood2 = 0
  for(i in 1:N){
    for(j in 1:p){
      if(rho_hat[j]  != 0){
        likelihood1 = likelihood1 + dpois(data[j, i], lambda = sg[j, i]*mu_hat[rho_hat[j], group_hat[i]+1], log = T)
      }else{
        likelihood2 = likelihood2 + dpois(data[j, i], lambda = sg[j, i]*mu0_hat[j], log = T)
      }
    }
  }
  term1 = (-2*(likelihood1+posterior_rho+posterior_z) + c*(K*log(N) + L*log(p_hat) + (K * L) *log(N*p_hat)))
  term2 = weight * (-2*likelihood2 + (p-p_hat)*log(N))
  #likelihood = (likelihood1 / p_hat + likelihood2 / (p - p_hat))*(p)
  #likelihood = likelihood1 + likelihood2
  #result = -2 *likelihood + N*log(K) + p_hat*log(L) + (K * L) *log(N*p_hat) + (p - p_hat)*log(N)
  result = (term1 + term2 )/2
  return(result)
}


