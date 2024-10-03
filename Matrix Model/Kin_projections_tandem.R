


source(here::here("Matrix Model", "Functions_required.R" ))

##################### Not time varying (boundary time = 0 kin output)

all_kin_dy_tandem <- function(Uf, 
                              Um, 
                              Ff, 
                              Fm, 
                              alpha, ## alpha = sex ratio
                              na, ## na = number of ages
                              ns, ## ns = number of stages
                              Parity, 
                              sex_Focal, ## binary "F" or "M"
                              stage_Focal){
  
  n <- nrow(Uf) ## number of ages * stages for each sex
  
  ## Projection matrices:
  
  ## Uproj is a block diagonal matrix of block-structured Age*Stage matrices over sex...
  ## ...independently over sex transfers individuals across stage and up age
  Uproj <- as(block_diag_function(list(Uf, Um)),"sparseMatrix") ## ( block_diag_func in "Functions_required.R" )
  
  ## Fproj is a Sex-block-structured matrix of block-structured Age*Stage matrices where males and females BOTH reproduce (by stage)
  Fproj <- as(matrix(0, nrow = 2*n, ncol = 2*n, byrow = T), "sparseMatrix")
  Fproj[1:n, 1:n] <- (1-alpha)*Ff ## Ff is Age*Stage block structured giving rate at which females in age-stage produce individuals in age-stage
  Fproj[(n+1):(2*n), 1:n] <- alpha*Ff
  Fproj[1:n, (n+1):(2*n)] <- (1-alpha)*Fm ## Fm is Age*Stage block structured giving rate at which males in age-stage produce individuals in age-stage
  Fproj[(n+1):(2*n), (n+1):(2*n)] <- alpha*Fm 
  
  ## Fprojstar is a Sex-block-structured matrix of block-structured Age*Stage matrices where ONLY females reproduce
  Fprojstar <- as(matrix(0, nrow = 2*n, ncol = 2*n, byrow = T), "sparseMatrix") ## Block structured F_tilde
  Fprojstar[1:n, 1:n] <- (1-alpha)*Ff
  Fprojstar[(n+1):(2*n), 1:n] <- alpha*Ff
  
  ## The stable population structure is an age*stage*sex vector:
  ##                                                            1:n gives the female age*stage structure
  ##                                                            (1+n):2n gives the male age*stage structure
  population_age_stage_structure <- SD(Uproj + Fprojstar)
  
  ### Stable distribution of mothers needs adjusting if we work with parity
  if(Parity){
    stage_Focal <- 1
    
    population_age_stage_of_parenting <- pi_mix_parity2(Uf, Um, Ff, Fm, alpha, na, ns)
    mothers_age_stage <- population_age_stage_of_parenting[[2]]
    fathers_age_stage <- population_age_stage_of_parenting[[3]]
    
    mothers_age_dist <- population_age_stage_of_parenting[[4]]
    fathers_age_dist <- population_age_stage_of_parenting[[5]]
    
  }
  else{
    population_age_stage_of_parenting <- pi_mix(Uf, Um, Ff, Fm, alpha, na, ns)
    mothers_age_stage <- population_age_stage_of_parenting[[2]]
    fathers_age_stage <- population_age_stage_of_parenting[[3]]
    
    mothers_age_dist <- population_age_stage_of_parenting[[4]]
    fathers_age_dist <- population_age_stage_of_parenting[[5]]
    
  }
  
  ####################################### The dynamics of Kinship, starting with Focal who is no longer a unit vector
  
  ### Focal requires its own dynamic: G_tilde constructed below tracks Focal's age*stage advancement over the time-scale 
  f_t <- get_G(Uf, na, ns) ## get_G function in "Functions_required.R" 
  m_t <- get_G(Um, na, ns)
  G_tilde <- block_diag_function(list(f_t,m_t))
  X_Focal  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_Focal <- rep(0, 2*n)
  if(sex_Focal == "Female"){
    entry <- 1 + (stage_Focal-1)*na
    IC_Focal[entry] <- 1}
  else{
    entry <- n + 1 + (stage_Focal-1)*na
    IC_Focal[entry] <- 1
  }
  
  ### empty kin matrices for all of Focal's kin
  X_children <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_grand_children <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_great_grand_children <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_parents <- as(matrix(0, nrow = 2*n, ncol = na, byrow = T),"sparseMatrix")
  X_grand_parents <- as( matrix(0, nrow = 2*n, ncol = na, byrow = T) ,"sparseMatrix")
  X_older_sibs  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_younger_sibs  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_older_niece_nephew  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_younger_niece_nephew  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_older_aunt_uncle  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_younger_aunts_uncles <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_older_cousins  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_younger_cousins  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  
  
  ### Initial distributions for kin with non-zero deterministic initial conditions: 
  # Focal, parents, children, grand+great children, younger siblings, and younger nieces/nehpews
  X_Focal[,1] <- IC_Focal
  X_parents[, 1] <- mothers_age_stage
  
  ### projection all kin with deterministic initial conditions
  foreach(i = 1 : (na-1))%do%{
    X_Focal[,i+1] <- G_tilde%*%X_Focal[,i] 
    X_parents[, i+1] <- Uproj%*%X_parents[, i] 
    X_younger_sibs[,i+1] <- Uproj%*%X_younger_sibs[,i] + Fprojstar%*%X_parents[,i]
    X_younger_niece_nephew[,i+1] <- Uproj%*%X_younger_niece_nephew[,i] + Fproj%*%X_younger_sibs[,i]
    X_children[,i+1] <- Uproj%*%X_children[,i] + Fproj%*%X_Focal[,i] 
    X_grand_children[,i+1] <- Uproj%*%X_grand_children[,i] + Fproj%*%X_children[,i]
    X_great_grand_children[,i+1] <- Uproj%*%X_great_grand_children[,i] + Fproj%*%X_grand_children[,i]
  }
  
  ### IC for kin which are derived from above kin (Focal, parents, children, grand+great children, younger siblings, and younger nieces/nehpews): 
  # grand parents, older sibs, younger aunts/uncles, older nieces/nephews
  IC_f_grand_pars <- mothers_age_dist
  IC_m_grand_pars <- fathers_age_dist
  IC_older_sibs_f <- mothers_age_dist
  pi_younger_aunts_uncles_f <- mothers_age_dist
  pi_younger_aunts_uncles_m <- fathers_age_dist
  IC_older_niece_nephew_f <- mothers_age_dist
  foreach(ic = 1 : (na))%do%{
    X_grand_parents[, 1] <- X_grand_parents[, 1] + (IC_f_grand_pars[ic] + IC_m_grand_pars[ic])*X_parents[,ic] ## IC the sum of parents of Focal's parents,
    X_older_sibs[,1] <- X_older_sibs[,1] + IC_older_sibs_f[ic]*X_children[,ic]
    X_older_niece_nephew[,1] <- X_older_niece_nephew[,1] + IC_older_niece_nephew_f[ic]*X_grand_children[,ic]
    X_younger_aunts_uncles[,1] <- X_younger_aunts_uncles[,1] + (pi_younger_aunts_uncles_f[ic] + pi_younger_aunts_uncles_m[ic])*X_younger_sibs[,ic]
  }
  
  ### Projections of grand parenst, older sibs, younger aunts/uncles, older nieces/nephews
  foreach(i = 1: (na-1))%do%{
    X_grand_parents[, i+1] <- Uproj%*%X_grand_parents[, i] 
    X_older_sibs[,i+1] <- Uproj%*%X_older_sibs[,i] 
    X_older_niece_nephew[,i+1] <- Uproj%*%X_older_niece_nephew[,i] + Fproj%*%X_older_sibs[,i]
    X_younger_aunts_uncles[,i+1] <- Uproj%*%X_younger_aunts_uncles[,i] + Fprojstar%*%X_grand_parents[,i]
  }
  
  ### IC for kin which are derived from above kin (older sibs, younger aunts/uncles, older nieces/nephews):
  ## older unts/uncles, older cousins, younger cousins
  pi_older_aunt_uncle_f <- mothers_age_dist
  pi_older_aunt_uncle_m <- fathers_age_dist
  IC_older_cousins_f <- mothers_age_dist
  IC_older_cousins_m <- fathers_age_dist
  IC_younger_cousins_f <- mothers_age_dist
  IC_younger_cousins_m <- fathers_age_dist
  foreach(ic = 1 : (na-1))%do%{
    X_older_aunt_uncle[,1] <- X_older_aunt_uncle[,1] + (pi_older_aunt_uncle_f[ic] + pi_older_aunt_uncle_m[ic])*X_older_sibs[,ic]
    X_older_cousins[,1] <- X_older_cousins[,1] + (IC_older_cousins_f[ic] + IC_older_cousins_m[ic])*X_older_niece_nephew[,ic]
    X_younger_cousins[,1] <- X_younger_cousins[,1] + (IC_younger_cousins_f[ic] + IC_younger_cousins_m[ic])*X_younger_niece_nephew[,ic]
  }
  
  ## Projections of older unts/uncles, older cousins, younger cousins
  foreach(i = 1: (na-1))%do%{
    X_older_aunt_uncle[,i+1] <- Uproj%*%X_older_aunt_uncle[,i] 
    X_older_cousins[,i+1] <- Uproj%*%X_older_cousins[,i] + Fproj%*%X_older_aunt_uncle[,i]
    X_younger_cousins[,i+1] <- Uproj%*%X_younger_cousins[,i] + Fproj%*%X_younger_aunts_uncles[,i]
  }

  #### OUTPUT of all kin
  return(list(X_Focal,
              X_children,
              X_grand_children,
              X_great_grand_children,
              X_parents,
              X_grand_parents,
              X_older_sibs,
              X_younger_sibs,
              X_older_niece_nephew,
              X_younger_niece_nephew,
              X_older_aunt_uncle,
              X_younger_aunts_uncles,
              X_older_cousins,
              X_younger_cousins,
              population_age_stage_structure
  ))
}

########################################## Now a time varying kin output -- all annotations the same as above

all_kin_dy_TV_tandem <- function(Uf, 
                          Um, 
                          Ff, 
                          Fm, 
                          alpha, 
                          na, 
                          ns, 
                          Parity, 
                          sex_Focal, 
                          stage_Focal,
                          previous_kin_Focal,
                          prev_kin_children, 
                          prev_kin_grandchildren,
                          prev_kin_greatgrandchildren,
                          prev_kin_parents,
                          prev_kin_grand_parents,
                          prev_kin_older_sibs,
                          prev_kin_younger_sibs,
                          prev_kin_older_niece_nephew,
                          prev_kin_younger_niece_nephew,
                          prev_kin_older_aunts_uncles,
                          prev_kin_younger_aunts_uncles,
                          prev_kin_older_cousins,
                          prev_kin_younger_cousins, 
                          previous_population_age_stage_structure){
  
  n <- nrow(Uf)
  Uproj <- block_diag_function(list(Uf,Um))
  Fproj <- as(matrix(0, nrow = 2*n, ncol = 2*n, byrow = T),"sparseMatrix") 
  Fproj[1:n, 1:n] <- (1-alpha)*Ff
  Fproj[(n+1):(2*n), 1:n] <- alpha*Ff
  Fproj[1:n, (n+1):(2*n)] <- (1-alpha)*Fm
  Fproj[(n+1):(2*n), (n+1):(2*n)] <- alpha*Fm 
  Fprojstar <- as(matrix(0, nrow = 2*n, ncol = 2*n, byrow = T),"sparseMatrix") ## Block structured F_tilde
  Fprojstar[1:n, 1:n] <- (1-alpha)*Ff
  Fprojstar[(n+1):(2*n), 1:n] <- alpha*Ff
  
  population_age_stage_structure <- previous_population_age_stage_structure
  population_age_stage_structure <- population_age_stage_structure/sum(population_age_stage_structure)
  population_age_stage_structure_next <- (Uproj + Fprojstar)%*%population_age_stage_structure
  
  ### Stable distribution of mothers needs adjusting if we work with parity
  if(Parity){
    stage_Focal <- 1
    
    population_age_stage_of_parenting <- pi_mix_TV_parity(Ff, Fm, alpha, na, ns, population_age_stage_structure)
    mothers_age_stage <- population_age_stage_of_parenting[[2]]
    fathers_age_stage <- population_age_stage_of_parenting[[3]]
    
    mothers_age_dist <- population_age_stage_of_parenting[[4]]
    fathers_age_dist <- population_age_stage_of_parenting[[5]]
    
  }
  else{
    
    population_age_stage_of_parenting <- pi_mix_TV(Ff, Fm, alpha, na, ns, population_age_stage_structure)
    mothers_age_stage <- population_age_stage_of_parenting[[2]]
    fathers_age_stage <- population_age_stage_of_parenting[[3]]
    
    mothers_age_dist <- population_age_stage_of_parenting[[4]]
    fathers_age_dist <- population_age_stage_of_parenting[[5]]
    
  }

  ### Focal requires its own dynamic: G_tilde constructed below tracks Focal's age*stage advancement over the time-scale 
  f_t <- get_G(Uf, na, ns) ## get_G function in "Functions_required.R" 
  m_t <- get_G(Um, na, ns)
  G_tilde <- block_diag_function(list(f_t,m_t))
  X_Focal  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_Focal <- rep(0, 2*n)
  if(sex_Focal == "Female"){
    entry <- 1 + (stage_Focal-1)*na
    IC_Focal[entry] <- 1}
  else{
    entry <- n + 1 + (stage_Focal-1)*na
    IC_Focal[entry] <- 1
  }
  
  ### empty kin matrices for all of Focal's kin
  X_children <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_grand_children <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_great_grand_children <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_parents <- as(matrix(0, nrow = 2*n, ncol = na, byrow = T),"sparseMatrix")
  X_grand_parents <- as( matrix(0, nrow = 2*n, ncol = na, byrow = T) ,"sparseMatrix")
  X_older_sibs  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_younger_sibs  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_older_niece_nephew  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_younger_niece_nephew  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_older_aunt_uncle  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_younger_aunts_uncles <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_older_cousins  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  X_younger_cousins  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  
  ### Initial distributions for kin with non-zero deterministic initial conditions:
  ## Focal, parents, children, grand+great children, younger siblings, and younger nieces/nehpews
  X_Focal[,1] <- IC_Focal
  X_parents[, 1] <- mothers_age_stage
  ### projection all above kin with deterministic initial conditions
  foreach(i = 1 : (na-1))%do%{
    X_Focal[,i+1] <- G_tilde%*%previous_kin_Focal[,i] 
    X_parents[, i+1] <- Uproj%*%prev_kin_parents[, i] 
    X_younger_sibs[,i+1] <- Uproj%*%prev_kin_younger_sibs[,i] + Fprojstar%*%prev_kin_parents[,i]
    X_younger_niece_nephew[,i+1] <- Uproj%*%prev_kin_younger_niece_nephew[,i] + Fproj%*%prev_kin_younger_sibs[,i]
    X_children[,i+1] <- Uproj%*%prev_kin_children[,i] + Fproj%*%previous_kin_Focal[,i] 
    X_grand_children[,i+1] <- Uproj%*%prev_kin_grandchildren[,i] + Fproj%*%prev_kin_children[,i]
    X_great_grand_children[,i+1] <- Uproj%*%prev_kin_greatgrandchildren[,i] + Fproj%*%prev_kin_grandchildren[,i]
  }
  
  ### IC for kin which are derived from above kin (Focal, parents, children, grand+great children, younger siblings, and younger nieces/nehpews): 
  # grand parents, older sibs, younger aunts/uncles, older nieces/nephews
  IC_f_grand_pars <- mothers_age_dist
  IC_m_grand_pars <- fathers_age_dist
  pi_younger_aunts_uncles_f <- mothers_age_dist
  pi_younger_aunts_uncles_m <- fathers_age_dist
  IC_older_sibs_f <- mothers_age_dist
  IC_older_niece_nephew_f <- mothers_age_dist
  foreach(ic = 1 : (na))%do%{
    X_grand_parents[, 1] <- X_grand_parents[, 1] + (IC_f_grand_pars[ic] + IC_m_grand_pars[ic])*prev_kin_parents[,ic] ## IC the sum of parents of Focal's parents,
    X_older_sibs[,1] <- X_older_sibs[,1] + IC_older_sibs_f[ic]*prev_kin_children[,ic]
    X_older_niece_nephew[,1] <- X_older_niece_nephew[,1] + IC_older_niece_nephew_f[ic]*prev_kin_grandchildren[,ic]
    X_younger_aunts_uncles[,1] <- X_younger_aunts_uncles[,1] + (pi_younger_aunts_uncles_f[ic] + pi_younger_aunts_uncles_m[ic])*prev_kin_younger_sibs[,ic]
  }
  
  ### Projections of older sibs, younger aunts/uncles, older nieces/nephews
  foreach(i = 1: (na-1))%do%{
    X_grand_parents[, i+1] <- Uproj%*%prev_kin_grand_parents[, i] 
    X_older_sibs[,i+1] <- Uproj%*%prev_kin_older_sibs[,i] 
    X_older_niece_nephew[,i+1] <- Uproj%*%prev_kin_older_niece_nephew[,i] + Fproj%*%prev_kin_older_sibs[,i]
    X_younger_aunts_uncles[,i+1] <- Uproj%*%prev_kin_younger_aunts_uncles[,i] + Fprojstar%*%prev_kin_grand_parents[,i]
  }
  
  ### IC for kin which are derived from above kin (older sibs, younger aunts/uncles, older nieces/nephews):
  ## older unts/uncles, older cousins, younger cousins
  pi_older_aunt_uncle_f <- mothers_age_dist
  pi_older_aunt_uncle_m <- fathers_age_dist
  IC_older_cousins_f <- mothers_age_dist
  IC_older_cousins_m <- fathers_age_dist
  IC_younger_cousins_f <- mothers_age_dist
  IC_younger_cousins_m <- fathers_age_dist
  foreach(ic = 1 : (na-1))%do%{
    X_older_aunt_uncle[,1] <- X_older_aunt_uncle[,1] + (pi_older_aunt_uncle_f[ic] + pi_older_aunt_uncle_m[ic])*prev_kin_older_sibs[,ic]
    X_older_cousins[,1] <- X_older_cousins[,1] + (IC_older_cousins_f[ic] + IC_older_cousins_m[ic])*prev_kin_older_niece_nephew[,ic]
    X_younger_cousins[,1] <- X_younger_cousins[,1] + (IC_younger_cousins_f[ic] + IC_younger_cousins_m[ic])*prev_kin_younger_niece_nephew[,ic]
  }
  
  ## Projections of older unts/uncles, older cousins, younger cousins
  foreach(i = 1: (na-1))%do%{
    X_older_aunt_uncle[,i+1] <- Uproj%*%prev_kin_older_aunts_uncles[,i] 
    X_older_cousins[,i+1] <- Uproj%*%prev_kin_older_cousins[,i] + Fproj%*%prev_kin_older_aunts_uncles[,i]
    X_younger_cousins[,i+1] <- Uproj%*%prev_kin_younger_cousins[,i] + Fproj%*%prev_kin_younger_aunts_uncles[,i]
  }
  
  return(list(X_Focal,
              X_children,
              X_grand_children,
              X_great_grand_children,
              X_parents,
              X_grand_parents,
              X_older_sibs,
              X_younger_sibs,
              X_older_niece_nephew,
              X_younger_niece_nephew,
              X_older_aunt_uncle,
              X_younger_aunts_uncles,
              X_older_cousins,
              X_younger_cousins,
              population_age_stage_structure_next))
}