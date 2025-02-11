#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <Rmath.h>
#include <iostream>

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends("RcppArmadillo")]]


vec rdirichlet(const vec alpha){
  
  int n = alpha.n_elem;
  vec prob = vec(n, fill::zeros);
  for(int k = 0; k < n; k++){
    prob(k) = as_scalar(randg(1, distr_param(alpha(k), 1.0)));
  }
  
  prob = prob / sum(prob);
  return prob;
}

// calculate the density of Poisson distribution
double dpois(int x, double lambda, bool log_value = true){
  if(lambda <= 0){
    throw std::logic_error("lambda should be positive");
  }
  double r;
  if(x == 0){
    r = - lambda;
  }else{
    vec s1 = arma::linspace(1, x, x);
    r = -sum(log(s1)) + x * log(lambda) - lambda;
  }
  if(log_value){return r;}
  else{return exp(r);}
}


// update parameter \pi_{0}
void update_pi(double &pi_t, const vec &rho_t, const double a0, const double b0){
  int p = rho_t.n_elem;
  uvec index = find(rho_t == 0);
  int q = index.n_elem;
  vec alpha = {a0 + q, b0 + p - q};
  pi_t = rdirichlet(alpha)(0);
}


//update mu
mat update_mu(const vec &rho_t, const vec &group_t, const int &K, const int &L,
               const mat &X, const mat &sg, const double &a, const double b){
  mat mu_t = zeros<mat>(L, K);
  for(int k=0; k < K; k++){
    for(int l = 0; l < L; l++){
      uvec index_col = find(group_t == k);
      uvec index_row = find(rho_t == l+1);
      mat x_k_l = X(index_row, index_col);
      mat sg_k_l = sg(index_row, index_col);
      double term1 = accu(x_k_l) + a ;
      double term2 = accu(sg_k_l) + b;
      mu_t(l, k) = as_scalar(randg(1, distr_param(term1, 1.0/term2)));
    }
  }
  return mu_t;
}



//update mu0
void update_mu0(vec &mu0_t, const vec &rho_t, const mat &X, const mat &sg, const double &a, const double b){
  uvec index_row = find(rho_t == 0);
  if(index_row.n_elem != 0){
    int L0 = index_row.n_elem;
    for(int j=0; j < L0; j++){
      vec x_j = X.row(index_row(j)).t();
      vec sg_j = sg.row(index_row(j)).t();
      double term1 = accu(x_j) + a ;
      double term2 = accu(sg_j) + b;
      mu0_t(index_row(j)) = as_scalar(randg(1, distr_param(term1, 1.0/term2)));
    }
  }
}

// [[Rcpp::export]]
double prob_gene_existing(const int &l, const int &j, const vec &rho_t, const vec &group_t, const mat &X, const mat &sg,
                          const double &pi, const double &beta, const double &alpha0, const double &beta0){
  
  vec clusters_name = unique(group_t);
  int N = group_t.n_elem;
  int K = clusters_name.n_elem;
  int p = X.n_rows;
  uvec index_0 = find(rho_t == 0);
  int q = index_0.n_elem;
  vec x_j = X.row(j).t();
  vec s_j = sg.row(j).t();
  
  double output = 0;
  if(l == 0){
    double term1 = sum(x_j);
    double term2 = sum(s_j);
    output += log(pi) + alpha0 * log(beta0) + lgamma(alpha0 + term1) - 
      lgamma(alpha0) - (alpha0 + term1) * log(beta0 + term2);
  }else{
    uvec index_l = find(rho_t == l);
    int n_l = index_l.n_elem;
    mat X_l = X.rows(index_l);
    mat sg_l = sg.rows(index_l);
    vec x_j = X.row(j).t();
    vec s_j = sg.row(j).t();
    double alpha1;
    double beta1;
    if(rho_t(j) == l){
      n_l = n_l - 1;
      for(int k=0; k < K; k++){
        uvec index_k = find(group_t == k);
        mat X_l_k = X_l.cols(index_k);
        mat sg_l_k = sg_l.cols(index_k);
        vec x_j_k = x_j(index_k);
        vec s_j_k = s_j(index_k);
        double term1 = sum(x_j_k);
        double term2 = sum(s_j_k);
        
        alpha1 = alpha0 + accu(X_l_k) - sum(x_j_k);
        beta1 = beta0 + accu(sg_l_k) - sum(s_j_k);
  
        output +=  alpha1 * log(beta1) + lgamma(alpha1 + term1) - lgamma(alpha1) - (alpha1 + term1) * log(beta1 + term2);
      }
    }else{
      for(int k=0; k < K; k++){
        uvec index_k = find(group_t == k);
        mat X_l_k = X_l.cols(index_k);
        mat sg_l_k = sg_l.cols(index_k);
        vec x_j_k = x_j(index_k);
        vec s_j_k = s_j(index_k);
        
        double term1 = sum(x_j_k);
        double term2 = sum(s_j_k);
        alpha1 = alpha0 + accu(X_l_k);
        beta1 = beta0 + accu(sg_l_k);
        output +=  alpha1 * log(beta1) + lgamma(alpha1 + term1) - lgamma(alpha1) - (alpha1 + term1) * log(beta1 + term2);
      }
    }
    
    output += log(1 - pi) + log((n_l + 0.001) / (p - q + beta));
  }
  return output;
}


double prob_gene_new(const int &L, const int &j, const vec &rho_t, const vec &group_t, const mat &X, const mat &sg,
                     const double &pi, const double &beta, const double &alpha0, const double &beta0){
  vec x_j = X.row(j).t();
  vec s_j = sg.row(j).t();
  vec clusters_name = unique(group_t);
  int K = clusters_name.n_elem;
  int p = X.n_rows;
  uvec index_0 = find(rho_t == 0);
  int q = index_0.n_elem;
  double output = 0;
  for(int k = 0; k < K; k++){
    uvec index_k = find(group_t == k);
    output += alpha0 * log(beta0) + lgamma(alpha0 + sum(x_j(index_k))) - 
      lgamma(alpha0) - (alpha0 + sum(x_j(index_k))) * log(beta0 + sum(s_j(index_k)));
  }
  output += log(1 - pi) + log(beta / (p - q + beta));
  return output-999999;
}


// [[Rcpp::export]]
void update_rho(vec &rho_t, const vec &group_t, const mat &X, const mat &sg, const double &pi, 
                const double &beta, const double &alpha0, const double &beta0, const int &Lmax, const double &temperature = 1.0){
  vec clusters_name = unique(rho_t);
  int L = clusters_name.n_elem;
  int p = X.n_rows;
  
  for(int j = 0; j < p; j++){
    int label_old = rho_t(j);
    int label_new;
    uvec group_j = find(rho_t == label_old);
    int c_size = group_j.n_elem;
    
    //not singleton case
    if(c_size > 1){
      vec prob_j = zeros<vec>(L+1);
      for(int l = 0; l < L; l++){
        prob_j(l) = prob_gene_existing(l, j, rho_t, group_t, X, sg, pi, beta, alpha0, beta0);
      }
      
      prob_j(L) = prob_gene_new(L, j, rho_t, group_t, X, sg, pi, beta, alpha0, beta0);
      
      //cout << prob_j << "\n";
      
      prob_j = (prob_j - max(prob_j))/temperature;
      prob_j = exp(prob_j) / sum(exp(prob_j));
      
      label_new = as_scalar(Rcpp::RcppArmadillo::sample(regspace(0, L), 1, false, prob_j));
      
      // cout << "label_new = " << label_new << "\n";
      
      if(label_new >= L && L < Lmax){
        
        rho_t(j) = label_new;
        clusters_name = unique(rho_t);
        L = clusters_name.n_elem;
      }
      
      else if(label_new >= L && L >= Lmax){
        
        rho_t(j) = label_old;
        clusters_name = unique(rho_t);
        L = clusters_name.n_elem;
      }
      else{
        rho_t(j) = label_new;
        clusters_name = unique(rho_t);
        L = clusters_name.n_elem;
      }
    }
    
    
    // singleton case
    else{
      if(L == 2){continue;}
      vec prob_j = zeros<vec>(L+1);
      for(int l = 0; l < L; l++){
        prob_j(l) = prob_gene_existing(l, j, rho_t, group_t, X, sg, pi, beta, alpha0, beta0);
      }
      
      prob_j(L) = prob_gene_new(L, j, rho_t, group_t, X, sg, pi, beta, alpha0, beta0);
      
      
      prob_j = (prob_j - max(prob_j))/temperature;
      prob_j = exp(prob_j) / sum(exp(prob_j));
      
      label_new = as_scalar(Rcpp::RcppArmadillo::sample(regspace(0, L), 1, false, prob_j));
      
      if(label_new == label_old || label_new == L){
        rho_t(j) = label_old;
        clusters_name = unique(rho_t);
        L = clusters_name.n_elem;
      }
      //else{
      //  uvec index = find(clusters_name != label_old);
      //  rho_t(j) = label_new;
      //  for( int id = 0; id < p; id++){
      //    if(rho_t(id) > label_old){
      //      rho_t(id) += -1;
      //    }
      //  }
      //  clusters_name = unique(rho_t);
      //  L = clusters_name.n_elem;
      //}
    }
  }
}


//function to calculate MRF energy function
double MRF_energy(const int &k, const int &l, const vec &group_t, const umat &G,
                  const double f, const bool log_value = true){
  
  uvec G_l = G.row(l).t();
  uvec neighbor = G_l.elem(find(G_l > 0)) - 1;
  
  uvec neighbor_id = find(group_t(neighbor) == k);
  double result = f * neighbor_id.n_elem;
  
  if(log_value){return result;}
  else{ return exp(result); }
}


//function to calculate factors for a new cluster
vec VN(const int Kmax, const int N, const double gamma, const double lambda=1){
  vec vn = vec(Kmax, fill::zeros);
  int iters = max(N+100, 1000);
  for(int k = 1; k < Kmax +1 ; k++){
    double r = -datum::inf;
    for(int t = k; t < iters + 1; t++){
      double b = 0;
      vec s1 = arma::linspace(t-k+1, t, k);
      b +=  sum(log(s1));
      vec s2 = arma::linspace(t*gamma, t*gamma + N -1, N);
      b += - sum(log(s2));
      double s3 = dpois(t-1, lambda);
      b +=  s3;
      double m = max(b, r);
      r = log(exp(r-m) + exp(b-m)) + m;
    }
    vn(k-1) = r;
  }
  return vn;
}



double prob_spot_existing(const int &k, const int &i, const vec &group_t, const vec &rho_t, 
                          const mat &X, const mat &sg, const double &alpha0, const double &beta0, 
                          const umat &G, const double &f, const double &GAMMA){
  vec gene_groups_name = unique(rho_t);
  int L = gene_groups_name.n_elem;
  double output = 0;
  uvec index_k = find(group_t == k);
  int n_k = index_k.n_elem;
  double alpha1;
  double beta1;
  if(group_t(i) == k){
    n_k = n_k - 1;
    for(int l = 1; l < L; l++){
      uvec index_l = find(rho_t == l);
      mat X_l_k = X(index_l, index_k);
      mat sg_l_k = sg(index_l, index_k);
      vec x_i = X.col(i);
      vec s_i = sg.col(i);
      double term1 = sum(x_i(index_l));
      double term2 = sum(s_i(index_l));
      alpha1 = alpha0 + accu(X_l_k) - term1;
      beta1 = beta0 + accu(sg_l_k) - term2;
      output +=  alpha1 * log(beta1) + lgamma(alpha1 + term1) - 
        lgamma(alpha1) - (alpha1 + term1) * log(beta1 + term2);
    }
  }else{
    for(int l = 1; l < L; l++){
      uvec index_l = find(rho_t == l);
      mat X_l_k = X(index_l, index_k);
      mat sg_l_k = sg(index_l, index_k);
      vec x_i = X.col(i);
      vec s_i = sg.col(i);
      double term1 = sum(x_i(index_l));
      double term2 = sum(s_i(index_l));
      alpha1 = alpha0 + accu(X_l_k);
      beta1 = beta0 + accu(sg_l_k);
      output +=  alpha1 * log(beta1) + lgamma(alpha1 + term1) - 
        lgamma(alpha1) - (alpha1 + term1) * log(beta1 + term2);
    }
  }
  output += log(n_k + GAMMA) + MRF_energy(k, i, group_t, G, f);
  return output;
}


double prob_spot_new(const int &K, const int &i, const vec &rho_t, const vec &vn, const mat &X,
                     const mat &sg, const double &alpha0, const double &beta0, const double &GAMMA){
  vec gene_groups_name = unique(rho_t);
  int L = gene_groups_name.n_elem;
  vec x_i = X.col(i);
  vec s_i = sg.col(i);
  double h = 0;
  for(int l = 1; l < L; l++){
    uvec index_l = find(rho_t == l);
    vec x_i_l = x_i(index_l);
    vec s_i_l = s_i(index_l);
    double term1 = sum(x_i_l);
    double term2 = sum(s_i_l);
    h += alpha0 * log(beta0) + lgamma(alpha0 + term1) - 
      lgamma(alpha0) - (alpha0 + term1) * log(beta0 + term2);
  }
  double output = max((h + log(GAMMA) + vn(K) - vn(K-1)), -999999.0);
  return output-999999;
}

// [[Rcpp::export]]
void update_group(vec &group_t, const vec &rho_t, const mat &X, const mat &sg, const vec &vn,
                  const double &alpha0, const double &beta0, 
                  const umat &G, const double &f, const double &GAMMA, const double &temperature){
  int Kmax = vn.n_elem - 1;
  int N = X.n_cols;
  vec clusters_name = unique(group_t);
  int K = clusters_name.n_elem;
  
  
  for(int i = 0; i < N; i++){
    int label_old = group_t(i);
    int label_new;
    uvec group_i = find(group_t == label_old);
    int c_size = group_i.n_elem;
    
    // cout << "i = " << i << "\n";
    
    if(c_size > 1){
      vec prob_i = zeros<vec>(K+1);
      for(int k = 0; k < K; k++){
        prob_i(k) = prob_spot_existing(k, i, group_t, rho_t, X, sg, alpha0, beta0, G, f, GAMMA);
      }
      
      prob_i(K) = prob_spot_new(K, i, rho_t, vn, X, sg, alpha0, beta0, GAMMA);
      
       
      prob_i = (prob_i - max(prob_i))/temperature;
      prob_i = exp(prob_i) / sum(exp(prob_i));
      
      label_new = as_scalar(Rcpp::RcppArmadillo::sample(regspace(0, K), 1, false, prob_i));
      
      // cout << "label_new = " << label_new << "\n";
      
      if(label_new >= K && K < Kmax){
        
        group_t(i) = label_new;
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
      }
      
      else if(label_new >= K && K >= Kmax){
        
        group_t(i) = label_old;
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
      }
      else{
        group_t(i) = label_new;
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
      }
      
    }
    
    else{
      
      //The cluster is a singleton, only K choices
      vec prob_i = vec(K+1, fill::zeros);
      for(int k=0; k < K; k++){
        prob_i(k) = prob_spot_existing(k, i, group_t, rho_t, X, sg, alpha0, beta0, G, f, GAMMA);
      }
      
      prob_i(K) = prob_spot_new(K, i, rho_t, vn, X, sg, alpha0, beta0, GAMMA);
      
      prob_i = (prob_i - max(prob_i))/temperature;
      prob_i = exp(prob_i) / sum(exp(prob_i));
      
      label_new = as_scalar(Rcpp::RcppArmadillo::sample(regspace(0, K), 1, false, prob_i));
      //cout << "else " << label_new << "\n";
      
      if(label_new == label_old || label_new == K){
        group_t(i) = label_old;
        clusters_name = unique(group_t);
        K = clusters_name.n_elem;
      }
      //else{
      //  uvec index = find(clusters_name != label_old);
      //  group_t(i) = label_new;
      //  for( int id = 0; id < N; id++){
      //    if(group_t(id) > label_old){
      //      group_t(id) += -1;
      //    }
      //   }
      //  clusters_name = unique(group_t);
      //  K = clusters_name.n_elem;
      //  
      //}
    }
  }
}


// [[Rcpp::export]]
Rcpp::List runMCMC(vec &group_t,  vec &rho_t, mat &mu_t, vec &mu0_t, double &pi_t, const mat &X, const mat &sg,
                   const double &alpha0, const double &beta0,  const double &a0, const double &b0, 
                   const double &f, const umat G, const double &BETA, const double &GAMMA, 
                   const int max_iters, const int seed = 12569){

  int N = X.n_cols;
  int p = X.n_rows;
  int Kmax = 30;
  int Lmax = 60;
  vec vn = VN(Kmax + 1, N, GAMMA);
  
  arma_rng::set_seed(seed);
  
  //store the parameters and log prob
  mat group_iter = zeros<mat>(N, max_iters);
  mat rho_iter = zeros<mat>(p, max_iters);
  cube mu_iter = zeros<cube>(Lmax, Kmax, max_iters);
  mat mu0_iter = zeros<mat>(p, max_iters);
  vec K_iter = zeros<vec>(max_iters);
  vec L_iter = zeros<vec>(max_iters);
  vec pi_iter = zeros<vec>(max_iters);
  
  // store initialized value
  vec spot_clusters = unique(group_t);
  int K = spot_clusters.n_elem;
  group_iter.col(0) = group_t;
  vec gene_groups = unique(rho_t);
  int L = gene_groups.n_elem;
  rho_iter.col(0) = rho_t;
  mu_iter.slice(0).submat(0, 0, L-2, K-1) = mu_t;
  mu0_iter.col(0) = mu0_t;
  
  uvec index_0 = find(rho_t == 0);
  //mu0_iter(index_0, 0) = mu0_t;
  K_iter(0) = K;
  L_iter(0) = L;
  pi_iter(0) = pi_t;
  
  for(int t = 1; t < max_iters; t++){
    double temperature = pow(0.9999, t);
    
    update_pi(pi_t, rho_t, a0, b0);
    
    // cout << pi_t << "\n";
    
    update_rho(rho_t, group_t, X, sg, pi_t, BETA, alpha0, beta0, Lmax, temperature = 1.0);
    gene_groups = unique(rho_t);
    L = gene_groups.n_elem;
    
    update_group(group_t,rho_t, X, sg, vn, alpha0, beta0, G, f, GAMMA, temperature);
    spot_clusters = unique(group_t);
    K = spot_clusters.n_elem;
    
    mu_t = update_mu(rho_t, group_t, K, L-1, X, sg, alpha0, beta0);
    
    update_mu0(mu0_t, rho_t, X, sg, alpha0, beta0);
    
    //cout << "L="<< L << "; " <<"K="<< K << '\n';
    
    group_iter.col(t) = group_t;
    rho_iter.col(t) = rho_t;
    mu_iter.slice(t).submat(0, 0, L-2, K-1) = mu_t;
    mu0_iter.col(t) = mu0_t;
    K_iter(t) = K;
    L_iter(t) = L;
    pi_iter(t) = pi_t;
    
  }
  return List::create(Named("K_iter") = K_iter,
                      Named("group_iter") = group_iter, 
                      Named("L_iter") = L_iter,
                      Named("rho_iter") = rho_iter, 
                      Named("mu_iter") = mu_iter,
                      Named("mu0_iter") = mu0_iter,
                      Named("pi_iter") = pi_iter);
}
  
  

  
