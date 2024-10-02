## Define a function which composes a direct sum from a list of matrices (ie. a block diagonal matrix)
block_diag_function <- function(mat_list){
  s = length(mat_list)
  u1 = mat_list[[1]]
  dims <- dim(u1)
  r = dims[1]
  diagmat <- as(matrix(0, nrow = r*s, ncol = r*s, byrow = TRUE),"sparseMatrix")
  for(i in 1:s){
    diagmat = diagmat + kronecker(E_matrix(i,i,s,s), mat_list[[i]])
  }
  return(diagmat)
}

## Mixing distributions in the 2-sex model (note this is a joint age-stage dist. We must extract the age-dist below)
pi_mix <- function(Uf, Um, Ff, Fm, alpha, na, ns){
  n <- length(Uf[1,])
  F_block <- as(matrix(0, nrow = n*2, ncol = n*2),"sparseMatrix")
  F_block[1:n, 1:n] <- (1-alpha)*Ff
  F_block[ (1+n):(2*n), 1:n] <- alpha*Ff
  A <- block_diag_function(list(Uf,Um)) + F_block
  stable_dist_vec <- SD(A)
  pi_f <-  t(rep(1, na*ns)%*%Ff)*stable_dist_vec[1:n] 
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  t(rep(1, na*ns)%*%Fm)*stable_dist_vec[(1+n):(2*n)] 
  pi_m <- pi_m / abs(sum(pi_m))
  return(list(c(pi_f,pi_m),pi_f,pi_m))
}

pi_mix_parity <- function(Uf, Um, Ff, Fm, alpha, na, ns){
  n <- length(Uf[1,])
  F_block <- as(matrix(0, nrow = n*2, ncol = n*2),"sparseMatrix")
  F_block[1:n, 1:n] <- (1-alpha)*Ff
  F_block[ (1+n):(2*n), 1:n] <- alpha*Ff
  A <- block_diag_function(list(Uf,Um)) + F_block
  stable_dist_vec <- SD(A)
  pi_f <-  t(rep(1, na*ns)%*%Ff)*stable_dist_vec[1:n]
  pi_f[seq(1, length(pi_f), ns)] <- 0
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  t(rep(1, na*ns)%*%Fm)*stable_dist_vec[(1+n):(2*n)]
  pi_m[seq(1, length(pi_m), ns)] <- 0
  pi_m <- pi_m / abs(sum(pi_m))
  return(list(c(pi_f,pi_m),pi_f,pi_m))
}


pi_mix_parity2 <- function(Uf, Um, Ff, Fm, alpha, na, ns){
  n <- length(Uf[1,])
  F_block <- as(matrix(0, nrow = n*2, ncol = n*2),"sparseMatrix")
  F_block[1:n, 1:n] <- (1-alpha)*Ff
  F_block[ (1+n):(2*n), 1:n] <- alpha*Ff
  A <- block_diag_function(list(Uf,Um)) + F_block
  stable_dist_vec <- SD(A)
  pi_f <-  t(rep(1, na*ns)%*%Ff)*stable_dist_vec[1:n]
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  t(rep(1, na*ns)%*%Fm)*stable_dist_vec[(1+n):(2*n)]
  pi_m <- pi_m / abs(sum(pi_m))
  
  Z <- diag(1, ns)
  Z[1,1] <- 0
  Iom <- diag(1, na)
  onesa  <- t(rep(1,na))
  ones <- t(rep(1,ns))
  momarray <- pi_f %*% onesa
  dadarray <- pi_m %*% onesa
  
  piage_f  <- kronecker(Iom,ones) %*% pi_f
  piage_m <- kronecker(Iom,ones) %*% pi_m
  for(i in 1:na){
    E <- Iom[,i] %*% t(Iom[i,])
    momarray[,i] <- kronecker(E,Z) %*% momarray[,i]
    dadarray[,i] <- kronecker(E,Z) %*% dadarray[,i]
  }
  momarray <- momarray %*% MASS::ginv(diag(colSums(momarray)))
  dadarray <- dadarray %*% MASS::ginv(diag(colSums(dadarray)))
  # no 0 parity mothers: (momarray %*% piage)[seq(1,600,6)]
  pi_f <- momarray %*% piage_f
  pi_m <- dadarray %*% piage_m
  
  
  return(list(c(pi_f,pi_m), pi_f, pi_m, piage_f, piage_m))
}


pi_mix_TV <- function(Ff, Fm, alpha, na, ns, previous_age_stage_dist){
  n <- length(Ff[1,])
  pi_f <-  t(rep(1,na*ns)%*%Ff)*previous_age_stage_dist[1:n] 
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  t(rep(1,na*ns)%*%Fm)*previous_age_stage_dist[(1+n):(2*n)] 
  pi_m <- pi_m / abs(sum(pi_m))
  return(list(c(pi_f,pi_m),pi_f,pi_m))
}

pi_mix_TV_parity <- function(Ff, Fm, alpha, na, ns, previous_age_stage_dist){
  n <- length(previous_age_stage_dist)
  pi_f <-  t(rep(1,na*ns)%*%Ff)*previous_age_stage_dist[1:(n/2)]
  pi_f <- pi_f / abs(sum(pi_f))
  pi_m <-  t(rep(1,na*ns)%*%Fm)*previous_age_stage_dist[(1+n/2):(n)]
  pi_m <- pi_m / abs(sum(pi_m))
  
  Z <- diag(1, ns)
  Z[1,1] <- 0
  Iom <- diag(1, na)
  onesa  <- t(rep(1,na))
  ones <- t(rep(1,ns))
  momarray <- pi_f %*% onesa
  dadarray <- pi_m %*% onesa
  
  piage_f  <- kronecker(Iom,ones) %*% pi_f
  piage_m <- kronecker(Iom,ones) %*% pi_m
  for(i in 1:na){
    E <- Iom[,i] %*% t(Iom[i,])
    momarray[,i] <- kronecker(E,Z) %*% momarray[,i]
    dadarray[,i] <- kronecker(E,Z) %*% dadarray[,i]
  }
  momarray <- momarray %*% MASS::ginv(diag(colSums(momarray)))
  dadarray <- dadarray %*% MASS::ginv(diag(colSums(dadarray)))
  # no 0 parity mothers: (momarray %*% piage)[seq(1,600,6)]
  pi_f <- momarray %*% piage_f
  pi_m <- dadarray %*% piage_m
  
  return(list(c(pi_f,pi_m),pi_f,pi_m,piage_f,piage_m))
}

## Age-distributions as extracted from above
pi_age <- function(Uf, Um, Ff, Fm, alpha, no_ages, no_stages){
  pi_F <- kronecker( diag(no_ages), t(rep(1, no_stages)) )%*%(pi_mix(Uf, Um, Ff, Fm, alpha, no_ages, no_stages)[[2]])
  pi_M <- kronecker( diag(no_ages), t(rep(1, no_stages)) )%*%(pi_mix(Uf, Um, Ff, Fm, alpha, no_ages, no_stages)[[3]])
  return(list(pi_F,pi_M))
}

pi_age_parity <- function(Uf, Um, Ff, Fm, alpha, no_ages, no_stages, females_age_stage, males_age_stage){
  pi_F <- kronecker( diag(no_ages), t(rep(1, no_stages)) )%*%females_age_stage
  pi_M <- kronecker( diag(no_ages), t(rep(1, no_stages)) )%*%males_age_stage
  return(list(pi_F,pi_M))
}

pi_age_TV <- function(pop_fem, pop_male, no_ages, no_stages){
  pi_F <- kronecker( diag(no_ages), t(rep(1, no_stages)) )%*%(pop_fem)
  pi_M <- kronecker( diag(no_ages), t(rep(1, no_stages)) )%*%(pop_male)
  return(list(pi_F,pi_M))
}



## A matrix which projects Focal over age and stages
get_G <- function(U, no_ages, no_stages){
  sig <- t(rep(1,no_ages*no_stages))%*%U
  diag <- diag(sig[1,])
  G <- U%*%ginv(diag)
  return(G)
}

# A file of general functions

# Let PM be the projection matrix

# The growth rate -- the spectral radius of PM

lambda <- function(PM) {
  lead_eig <- (abs(eigen(PM, only.values = TRUE)$values))
  lead_eig <- lead_eig[which.max(lead_eig)]
  return(lead_eig)}


# stable age/stage distribution
# as the (right) eigenvector associated with lambda
SD <- function(PM) {
  spectral_stuff <- eigen(PM)
  spectral_stuff <- Re(spectral_stuff$vectors[,which.max(abs(spectral_stuff$values))])
  # normalise...
  vec_lambda <- spectral_stuff/sum(spectral_stuff)
  return(vec_lambda)}

# reproductive values as the (left) eigenvector -- lambda
RD <- function(PM) {
  spectral_stuff <- eigen(t(PM))
  spectral_stuff <- Re(spectral_stuff$vectors[,which.max(abs(spectral_stuff$values))])
  # normalise...
  vec_lambda <- spectral_stuff/sum(spectral_stuff)
  return(vec_lambda)}


# For age by stage structured models

# The marginal stage distribution (i.e., summing over all ages)
# ags : number of ages considered, number stages, and full dist with stage 
# embedded within age
marg_stage_dist <- function(no_ages, no_stages, full_dist){
  return(kronecker(t(rep(1, no_ages)) , diag(no_stages))%*%full_dist)}

# The marginal age dist (i.e., summing over all stages)
marg_age_dist <- function(no_ages, no_stages, full_dist){
  return(kronecker(diag(no_ages) , t(rep(1, no_stages)) ) %*%full_dist)
}

## Matirx operations -- defining the vec permutation martix 

# Given A of dimension (n,m) the permuation K_nm : K_nm vec A = vec A'
# Here E_ij is a matrix of appropriate dimension (n,m) with 0's everywhere but at (i,j)

#

e_vector <- function(i, n){
  e <- rep(0, n)
  e[i] <- 1
  return(e)
}

E_matrix <- function(i,j,n,m){
  E <- matrix(0, nrow = n, ncol = m, byrow = TRUE)
  E[i,j] <- 1
  return(E)
  
}

K_perm_mat <- function(n,m){
  perm = matrix(0, nrow = n*m, ncol = n*m, byrow = TRUE)
  for(i in 1:n){
    for(j in 1:m){
      perm = perm + kronecker( E_matrix(i,j,n,m) ,t( E_matrix(i,j,n,m) ) )
    }
  }
  return(perm)
}



