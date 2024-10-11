

## Define a function which composes a direct sum from a list of matrices (ie. a block diagonal matrix)
block_diag_function <- function(mat_list){
  s = length(mat_list)
  u1 = mat_list[[1]]
  dims <- dim(u1)
  r = dims[1]
  diagmat <- Matrix::Matrix(nrow = (r*s), ncol = (r*s), data = 0, sparse = TRUE)
  for(i in 1:s){
      diagmat = diagmat + kronecker(E_matrix(i,i,s,s), mat_list[[i]])
  }
  return(diagmat)
}

## Mixing distributions in the 2-sex model (note this is a joint age-stage dist. We must extract the age-dist below)
## Non-parity case
pi_mix <- function(Uf, Um, Ff, Fm, alpha, na, ns){
  n <- length(Uf[1,])
  F_block <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE)
  F_block[1:n, 1:n] <- (1-alpha)*Ff
  F_block[ (1+n):(2*n), 1:n] <- alpha*Ff
  A <- block_diag_function(list(Uf,Um)) + F_block
  stable_dist_vec <- SD(A)
  ### Joint distributions
  pi_f <-  Matrix::t( rep(1, na*ns) %*% Ff )*stable_dist_vec[1:n] 
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  Matrix::t( rep(1, na*ns) %*% Fm )*stable_dist_vec[(1+n):(2*n)] 
  pi_m <- pi_m / abs(sum(pi_m))
  ### Age distributions
  pi_F <- kronecker( diag(na), Matrix::t(rep(1, ns)) ) %*% (pi_f)
  pi_M <- kronecker( diag(na), Matrix::t(rep(1, ns)) ) %*% (pi_m)
  return(list(c(pi_f,pi_m),pi_f,pi_m,pi_F,pi_M))
}
## Time-varying analogue
pi_mix_TV <- function(Ff, Fm, alpha, na, ns, previous_age_stage_dist){
  n <- length(Ff[1,])
  ### Joint distributions
  pi_f <-  Matrix::t( rep(1,na*ns) %*% Ff )*previous_age_stage_dist[1:n] 
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  Matrix::t( rep(1,na*ns) %*% Fm )*previous_age_stage_dist[(1+n):(2*n)] 
  pi_m <- pi_m / abs(sum(pi_m))
  ### Age distributions
  pi_F <- kronecker( diag(na), Matrix::t(rep(1, ns)) ) %*% (pi_f)
  pi_M <- kronecker( diag(na), Matrix::t(rep(1, ns)) ) %*% (pi_m)
  return(list(c(pi_f,pi_m),pi_f,pi_m,pi_F,pi_M))
}


## Parity case, as used in DemoKin
pi_mix_parity <- function(Uf, Um, Ff, Fm, alpha, na, ns){
  n <- length(Uf[1,])
  F_block <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE)
  F_block[1:n, 1:n] <- (1-alpha)*Ff
  F_block[ (1+n):(2*n), 1:n] <- alpha*Ff
  A <- block_diag_function(list(Uf,Um)) + F_block
  stable_dist_vec <- SD(A)
  pi_f <-  Matrix::t( rep(1, na*ns) %*% Ff )*stable_dist_vec[1:n]
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  Matrix::t( rep(1, na*ns) %*% Fm )*stable_dist_vec[(1+n):(2*n)]
  pi_m <- pi_m / abs(sum(pi_m))
  m_mat <- pi_f %*% Matrix::t(rep(1,na))
  d_mat <- pi_m %*% Matrix::t(rep(1,na))
  pi_F <- kronecker( diag(1, na), Matrix::t(rep(1,ns)) ) %*% pi_f
  pi_M <- kronecker( diag(1, na), Matrix::t(rep(1,ns)) ) %*% pi_m
  for(i in 1:na){
    m_mat[,i] <- kronecker( E_matrix(i,i,na,na) , Matrix::diag( c(0, rep(1, ns-1)) ) ) %*% m_mat[,i]
    d_mat[,i] <- kronecker( E_matrix(i,i,na,na) , Matrix::diag( c(0, rep(1, ns-1)) ) ) %*% d_mat[,i]
  }
  out_mum <- m_mat %*% MASS::ginv(Matrix::diag(Matrix::colSums(m_mat)))
  out_dad <- d_mat %*% MASS::ginv(Matrix::diag(Matrix::colSums(d_mat)))
  ### Joint distributions
  pi_f <- out_mum %*% pi_F
  pi_m <- out_dad %*% pi_M
  return(list(c(pi_f,pi_m), pi_f, pi_m, pi_F, pi_M))
}
## Time-variant analogue
pi_mix_TV_parity <- function(Ff, Fm, alpha, na, ns, previous_age_stage_dist){
  n <- length(Ff[1,])
  pi_f <-  Matrix::t( rep(1,na*ns) %*% Ff )*previous_age_stage_dist[1:n]
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  Matrix::t( rep(1,na*ns) %*% Fm )*previous_age_stage_dist[(1+n):(2*n)]
  pi_m <- pi_m / abs(sum(pi_m))
  m_mat <- pi_f %*% Matrix::t(rep(1,na))
  d_mat <- pi_m %*% Matrix::t(rep(1,na))
  pi_F  <- kronecker( Matrix::diag(1, na), Matrix::t(rep(1,ns)) ) %*% pi_f
  pi_M <- kronecker( Matrix::diag(1, na), Matrix::t(rep(1,ns)) ) %*% pi_m
  for(i in 1:na){
    m_mat[,i] <- kronecker( E_matrix(i,i,na,na) , Matrix::diag( c(0, rep(1, ns-1)) ) ) %*% m_mat[,i]
    d_mat[,i] <- kronecker( E_matrix(i,i,na,na) , Matrix::diag( c(0, rep(1, ns-1)) ) ) %*% d_mat[,i]
  }
  out_mum <- m_mat %*% MASS::ginv(Matrix::diag(Matrix::colSums(m_mat)))
  out_dad <- d_mat %*% MASS::ginv(Matrix::diag(Matrix::colSums(d_mat)))
  ### Joint distributions
  pi_f <- out_mum %*% pi_F
  pi_m <- out_dad %*% pi_M
  return(list(c(pi_f,pi_m),pi_f,pi_m,pi_F,pi_M))
}

######################################################### Some useful utility functions
## A matrix which projects Focal over age and stages
get_G <- function(U, na, ns){
  sig <- Matrix::t(rep(1,na*ns)) %*% U
  diag <- Matrix::diag(sig[1,])
  G <- U %*% MASS::ginv(diag)
  return(G)
}

# The growth rate -- the spectral radius of PM
lambda <- function(PM) {
  lead_eig <- (abs(eigen(PM, only.values = TRUE)$values))
  lead_eig <- lead_eig[which.max(lead_eig)]
  return(lead_eig)}

# stable age/stage distribution
SD <- function(PM) {
  spectral_stuff <- eigen(PM)
  spectral_stuff <- Re(spectral_stuff$vectors[, which.max(abs(spectral_stuff$values))])
  # normalise...
  vec_lambda <- spectral_stuff/sum(spectral_stuff)
  return(vec_lambda)}

# reproductive values as the (left) eigenvector -- lambda
RD <- function(PM) {
  spectral_stuff <- eigen(t(PM))
  spectral_stuff <- Re(spectral_stuff$vectors[, which.max(abs(spectral_stuff$values))])
  # normalise...
  vec_lambda <- spectral_stuff/sum(spectral_stuff)
  return(vec_lambda)}

## The marginal stage distribution (i.e., summing over all ages)
marg_stage_dist <- function(na, ns, full_dist){
  return(kronecker( Matrix::t(rep(1, na)) , Matrix::diag(ns) ) %*% full_dist)}

# The marginal age dist (i.e., summing over all stages)
marg_age_dist <- function(na, ns, full_dist){
  return(kronecker( Matrix::diag(na) , Matrix::t(rep(1, ns)) ) %*% full_dist)
}

## Matirx operations -- defining the vec permutation martix 

e_vector <- function(i, n){
  e <- rep(0, n)
  e[i] <- 1
  return(e)
}

E_matrix <- function(i,j,n,m){
  #E <- matrix(0, nrow = n, ncol = m, byrow = TRUE)
  E <- Matrix::Matrix(nrow = (n), ncol = (m), data = 0, sparse = TRUE)
  E[i,j] <- 1
  return(E)
  
}

K_perm_mat <- function(n,m){
  perm <- Matrix::Matrix(nrow = (n*m), ncol = (n*m), data = 0, sparse = TRUE)
  for(i in 1:n){
    for(j in 1:m){
      perm = perm + kronecker( E_matrix(i,j,n,m) , Matrix::t(E_matrix(i,j,n,m)) )
    }
  }
  return(perm)
}



