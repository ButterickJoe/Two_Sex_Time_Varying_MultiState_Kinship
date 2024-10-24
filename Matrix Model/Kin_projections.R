


source(here::here("Matrix Model", "matrix_operations.R" ))



##################### Not time varying (boundary time = 0 kin output)

#' Title time invariant two-sex multi-state kin projections
#'
#' @param Uf matrix (block structured). transfers female individuals across stages and advances their age (conditional on survial)
#' @param Um matrix (block structured). transfers male individuals across stages and advances their age (conditional on survial)
#' @param Ff matrix (block structured). accounts for female reproduction, and assigns newborns into given age*stage
#' @param Fm matrix (block structured). accounts for male reproduction; assigns newborns into age-class, and stage
#' @param alpha scalar. birth ratio (male:female)
#' @param na scalar. number of ages.
#' @param ns scalar. number of stages.
#' @param Parity logical. If true then we omit mothers of parity 0, and re-scale the mother's age*stage of parenting
#' @param sex_Focal logical. Female or Male
#' @param Initial_stage_Focal numeric. Any natural number {1,2,3,4,...}
#'
#' @return a list of matrices. Each list entry represents a particular kin. Each kin is chacacterised by a matrix of dimension:
#' nrow = 2* na * ns (2-sex age-stage structured) and ncol = na (Focal's age)
#' yielding the age*stage distribution of kin for each age of Focal
#'
#' @export
#'
all_kin_dy <- function(Uf,
                       Um,
                       Ff,
                       Fm,
                       alpha, ## alpha = sex ratio male:female (i.e., 1 - birth_female)
                       na, ## na = number of ages
                       ns, ## ns = number of stages
                       Parity,
                       sex_Focal, ## binary "F" or "M"
                       Initial_stage_Focal){
  
  n <- nrow(Uf) ## number of ages * stages for each sex
  
  ## Projection matrices:
  
  ## Uproj is a block diagonal matrix of block-structured Age*Stage matrices over sex...
  ## ...independently over sex transfers individuals across stage and up age
  Uproj <- Matrix::Matrix(block_diag_function(list(Uf, Um)), sparse = TRUE)
  ## Fproj is a Sex-block-structured matrix of block-structured Age*Stage matrices where males and females BOTH reproduce (by stage)
  Fproj <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE)
  Fproj[1:n, 1:n] <- (1-alpha)*Ff ## Ff is Age*Stage block structured giving rate at which females in age-stage produce individuals in age-stage
  Fproj[(n+1):(2*n), 1:n] <- alpha*Ff
  Fproj[1:n, (n+1):(2*n)] <- (1-alpha)*Fm ## Fm is Age*Stage block structured giving rate at which males in age-stage produce individuals in age-stage
  Fproj[(n+1):(2*n), (n+1):(2*n)] <- alpha*Fm
  
  ## Fprojstar is a Sex-block-structured matrix of block-structured Age*Stage matrices where ONLY females reproduce
  Fprojstar <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE) ## Block structured F_tilde
  Fprojstar[1:n, 1:n] <- (1-alpha)*Ff
  Fprojstar[(n+1):(2*n), 1:n] <- alpha*Ff
  
  ## The stable population structure is an age*stage*sex vector:
  ##                                                            1:n gives the female age*stage structure
  ##                                                            (1+n):2n gives the male age*stage structure
  population_age_stage_structure <- SD(Uproj + Fprojstar)
  
  ### Stable distribution of mothers needs adjusting if we work with parity
  if(Parity){
    Initial_stage_Focal <- 1
    
    population_age_stage_of_parenting <- pi_mix_parity(Uf, Um, Ff, Fm, alpha, na, ns)
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
  X_Focal  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  IC_Focal <- rep(0, 2*n)
  if(sex_Focal == "Female"){
    entry <- 1 + (Initial_stage_Focal-1)*na
    IC_Focal[entry] <- 1}
  else{
    entry <- n + 1 + (Initial_stage_Focal-1)*na
    IC_Focal[entry] <- 1
  }
  
  ### empty kin matrices for all of Focal's kin
  X_children <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_grand_children <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_great_grand_children <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_parents <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_grand_parents <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_sibs  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_sibs  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_niece_nephew  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_niece_nephew  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_aunt_uncle  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_aunts_uncles <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_cousins  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_cousins  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  
  
  ### Initial distributions for kin with non-zero deterministic initial conditions:
  # Focal, parents, children, grand+great children, younger siblings, and younger nieces/nehpews
  X_Focal[,1] <- IC_Focal
  X_parents[, 1] <- mothers_age_stage
  
  ### projection all kin with deterministic initial conditions
  for(i in 1 : (na-1)){
    X_Focal[,i+1] <- G_tilde %*% X_Focal[,i]
    X_parents[, i+1] <- Uproj %*% X_parents[, i]
    X_younger_sibs[,i+1] <- Uproj %*% X_younger_sibs[,i] + Fprojstar %*% X_parents[,i]
    X_younger_niece_nephew[,i+1] <- Uproj %*% X_younger_niece_nephew[,i] + Fproj %*% X_younger_sibs[,i]
    X_children[,i+1] <- Uproj %*% X_children[,i] + Fproj %*% X_Focal[,i]
    X_grand_children[,i+1] <- Uproj %*% X_grand_children[,i] + Fproj %*% X_children[,i]
    X_great_grand_children[,i+1] <- Uproj %*% X_great_grand_children[,i] + Fproj %*% X_grand_children[,i]
  }
  
  ### IC for kin which are derived from above kin (Focal, parents, children, grand+great children, younger siblings, and younger nieces/nehpews):
  # grand parents, older sibs, younger aunts/uncles, older nieces/nephews
  IC_f_grand_pars <- mothers_age_dist
  IC_m_grand_pars <- fathers_age_dist
  IC_older_sibs_f <- mothers_age_dist
  IC_younger_aunts_uncles_f <- mothers_age_dist
  IC_younger_aunts_uncles_m <- fathers_age_dist
  IC_older_niece_nephew_f <- mothers_age_dist
  for(ic in 1 : (na)){
    X_grand_parents[, 1] <- X_grand_parents[, 1] + (IC_f_grand_pars[ic] + IC_m_grand_pars[ic])*X_parents[,ic] ## IC the sum of parents of Focal's parents,
    X_older_sibs[,1] <- X_older_sibs[,1] + IC_older_sibs_f[ic]*X_children[,ic]
    X_older_niece_nephew[,1] <- X_older_niece_nephew[,1] + IC_older_niece_nephew_f[ic]*X_grand_children[,ic]
    X_younger_aunts_uncles[,1] <- X_younger_aunts_uncles[,1] + (IC_younger_aunts_uncles_f[ic] + IC_younger_aunts_uncles_m[ic])*X_younger_sibs[,ic]
  }
  
  ### Projections of grand parenst, older sibs, younger aunts/uncles, older nieces/nephews
  for(i in 1: (na-1)){
    X_grand_parents[, i+1] <- Uproj %*% X_grand_parents[, i]
    X_older_sibs[,i+1] <- Uproj %*% X_older_sibs[,i]
    X_older_niece_nephew[,i+1] <- Uproj %*% X_older_niece_nephew[,i] + Fproj %*% X_older_sibs[,i]
    X_younger_aunts_uncles[,i+1] <- Uproj %*% X_younger_aunts_uncles[,i] + Fprojstar %*% X_grand_parents[,i]
  }
  
  ### IC for kin which are derived from above kin (older sibs, younger aunts/uncles, older nieces/nephews):
  ## older unts/uncles, older cousins, younger cousins
  IC_older_aunt_uncle_f <- mothers_age_dist
  IC_older_aunt_uncle_m <- fathers_age_dist
  IC_older_cousins_f <- mothers_age_dist
  IC_older_cousins_m <- fathers_age_dist
  IC_younger_cousins_f <- mothers_age_dist
  IC_younger_cousins_m <- fathers_age_dist
  for(ic in 1 : (na-1)){
    X_older_aunt_uncle[,1] <- X_older_aunt_uncle[,1] + (IC_older_aunt_uncle_f[ic] + IC_older_aunt_uncle_m[ic])*X_older_sibs[,ic]
    X_older_cousins[,1] <- X_older_cousins[,1] + (IC_older_cousins_f[ic] + IC_older_cousins_m[ic])*X_older_niece_nephew[,ic]
    X_younger_cousins[,1] <- X_younger_cousins[,1] + (IC_younger_cousins_f[ic] + IC_younger_cousins_m[ic])*X_younger_niece_nephew[,ic]
  }
  
  ## Projections of older unts/uncles, older cousins, younger cousins
  for(i in 1: (na-1)){
    X_older_aunt_uncle[,i+1] <- Uproj %*% X_older_aunt_uncle[,i]
    X_older_cousins[,i+1] <- Uproj %*% X_older_cousins[,i] + Fproj %*% X_older_aunt_uncle[,i]
    X_younger_cousins[,i+1] <- Uproj %*% X_younger_cousins[,i] + Fproj %*% X_younger_aunts_uncles[,i]
  }
  
  #### OUTPUT of all kin
  return(list("Focal" = X_Focal,
              "d" = X_children,
              "gd" = X_grand_children,
              "ggd" = X_great_grand_children,
              "m" = X_parents,
              "gm" = X_grand_parents,
              "os" = X_older_sibs,
              "ys" = X_younger_sibs,
              "nos" = X_older_niece_nephew,
              "nys" = X_younger_niece_nephew,
              "oa" = X_older_aunt_uncle,
              "ya" = X_younger_aunts_uncles,
              "coa" = X_older_cousins,
              "cya" = X_younger_cousins,
              "ps" = population_age_stage_structure
  ))
}

########################################## Now a time varying kin output -- all annotations the same as above

#' Title time-variant two-sex multi-state kin projections
#'
#' @param Uf matrix (block structured). transfers female individuals across stages and advances their age (conditional on survial)
#' @param Um matrix (block structured). transfers male individuals across stages and advances their age (conditional on survial)
#' @param Ff matrix (block structured). accounts for female reproduction, and assigns newborns into given age*stage
#' @param Fm matrix (block structured). accounts for male reproduction; assigns newborns into age-class, and stage
#' @param alpha scalar. birth ratio (male:female)
#' @param na scalar. number of ages.
#' @param ns scalar. number of stages.
#' @param Parity logical. If true then we omit mothers of parity 0, and re-scale the mother's age*stage of parenting
#' @param sex_Focal logical. Female or Male
#' @param Initial_stage_Focal numeric. Any natural number {1,2,3,4,...}
#' @param previous_kin_Focal matrix. last years kinship output.
#' @param prev_kin_children matrix. last years kinship output.
#' @param prev_kin_grandchildren matrix. last years kinship output.
#' @param prev_kin_greatgrandchildren matrix. last years kinship output.
#' @param prev_kin_parents matrix. last years kinship output.
#' @param prev_kin_grand_parents matrix. last years kinship output.
#' @param prev_kin_older_sibs matrix. last years kinship output.
#' @param prev_kin_younger_sibs matrix. last years kinship output.
#' @param prev_kin_older_niece_nephew matrix. last years kinship output.
#' @param prev_kin_younger_niece_nephew matrix. last years kinship output.
#' @param prev_kin_older_aunts_uncles matrix. last years kinship output.
#' @param prev_kin_younger_aunts_uncles matrix. last years kinship output.
#' @param prev_kin_older_cousins matrix. last years kinship output.
#' @param prev_kin_younger_cousins matrix. last years kinship output.
#' @param previous_population_age_stage_structure vector. The transient "population structure" (age*stage distributed)
#'
#' @return a list of matrices. Each list entry represents a particular kin. Each kin is chacacterised by a matrix of dimension:
#' nrow = 2* na * ns (2-sex age-stage structured) and ncol = na (Focal's age)
#' yielding the age*stage distribution of kin for each age of Focal
#'
#' @export
all_kin_dy_TV <- function(Uf,
                          Um,
                          Ff,
                          Fm,
                          alpha, ## alpha = sex ratio male:female (i.e., 1 - birth_female)
                          na, ## number of ages
                          ns, ## number of stages
                          Parity,
                          sex_Focal,
                          Initial_stage_Focal,
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
  Uproj <- Matrix::Matrix(block_diag_function(list(Uf, Um)), sparse = TRUE)
  Fproj <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE)
  Fproj[1:n, 1:n] <- (1-alpha)*Ff
  Fproj[(n+1):(2*n), 1:n] <- alpha*Ff
  Fproj[1:n, (n+1):(2*n)] <- (1-alpha)*Fm
  Fproj[(n+1):(2*n), (n+1):(2*n)] <- alpha*Fm
  Fprojstar <- Matrix::Matrix(nrow = (2*n), ncol = (2*n), data = 0, sparse = TRUE) ## Block structured F_tilde
  Fprojstar[1:n, 1:n] <- (1-alpha)*Ff
  Fprojstar[(n+1):(2*n), 1:n] <- alpha*Ff
  
  population_age_stage_structure <- previous_population_age_stage_structure
  population_age_stage_structure <- population_age_stage_structure/sum(population_age_stage_structure)
  population_age_stage_structure_next <- (Uproj + Fprojstar)%*%population_age_stage_structure
  
  ### Stable distribution of mothers needs adjusting if we work with parity
  if(Parity){
    Initial_stage_Focal <- 1
    
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
  X_Focal  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  IC_Focal <- rep(0, 2*n)
  if(sex_Focal == "Female"){
    entry <- 1 + (Initial_stage_Focal-1)*na
    IC_Focal[entry] <- 1}
  else{
    entry <- n + 1 + (Initial_stage_Focal-1)*na
    IC_Focal[entry] <- 1
  }
  
  ### empty kin matrices for all of Focal's kin
  X_children <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_grand_children <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_great_grand_children <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_parents <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_grand_parents <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_sibs <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_sibs <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_niece_nephew <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_niece_nephew <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_aunt_uncle <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_aunts_uncles <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_older_cousins <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  X_younger_cousins <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  
  ### Initial distributions for kin with non-zero deterministic initial conditions:
  ## Focal, parents, children, grand+great children, younger siblings, and younger nieces/nehpews
  X_Focal[,1] <- IC_Focal
  X_parents[, 1] <- mothers_age_stage
  ### projection all above kin with deterministic initial conditions
  for(i in 1 : (na-1)){
    X_Focal[,i+1] <- G_tilde %*% previous_kin_Focal[,i]
    X_parents[, i+1] <- Uproj %*% prev_kin_parents[, i]
    X_younger_sibs[,i+1] <- Uproj %*% prev_kin_younger_sibs[,i] + Fprojstar %*% prev_kin_parents[,i]
    X_younger_niece_nephew[,i+1] <- Uproj %*% prev_kin_younger_niece_nephew[,i] + Fproj %*% prev_kin_younger_sibs[,i]
    X_children[,i+1] <- Uproj %*% prev_kin_children[,i] + Fproj %*% previous_kin_Focal[,i]
    X_grand_children[,i+1] <- Uproj %*% prev_kin_grandchildren[,i] + Fproj %*% prev_kin_children[,i]
    X_great_grand_children[,i+1] <- Uproj %*% prev_kin_greatgrandchildren[,i] + Fproj %*% prev_kin_grandchildren[,i]
  }
  
  ### IC for kin which are derived from above kin (Focal, parents, children, grand+great children, younger siblings, and younger nieces/nehpews):
  # grand parents, older sibs, younger aunts/uncles, older nieces/nephews
  IC_f_grand_pars <- mothers_age_dist
  IC_m_grand_pars <- fathers_age_dist
  IC_younger_aunts_uncles_f <- mothers_age_dist
  IC_younger_aunts_uncles_m <- fathers_age_dist
  IC_older_sibs_f <- mothers_age_dist
  IC_older_niece_nephew_f <- mothers_age_dist
  for(ic in 1 : (na)){
    X_grand_parents[, 1] <- X_grand_parents[, 1] + (IC_f_grand_pars[ic] + IC_m_grand_pars[ic])*prev_kin_parents[,ic] ## IC the sum of parents of Focal's parents,
    X_older_sibs[,1] <- X_older_sibs[,1] + IC_older_sibs_f[ic]*prev_kin_children[,ic]
    X_older_niece_nephew[,1] <- X_older_niece_nephew[,1] + IC_older_niece_nephew_f[ic]*prev_kin_grandchildren[,ic]
    X_younger_aunts_uncles[,1] <- X_younger_aunts_uncles[,1] + (IC_younger_aunts_uncles_f[ic] + IC_younger_aunts_uncles_m[ic])*prev_kin_younger_sibs[,ic]
  }
  
  ### Projections of older sibs, younger aunts/uncles, older nieces/nephews
  for(i in 1: (na-1)){
    X_grand_parents[, i+1] <- Uproj %*% prev_kin_grand_parents[, i]
    X_older_sibs[,i+1] <- Uproj %*% prev_kin_older_sibs[,i]
    X_older_niece_nephew[,i+1] <- Uproj %*% prev_kin_older_niece_nephew[,i] + Fproj %*% prev_kin_older_sibs[,i]
    X_younger_aunts_uncles[,i+1] <- Uproj %*% prev_kin_younger_aunts_uncles[,i] + Fprojstar %*% prev_kin_grand_parents[,i]
  }
  
  ### IC for kin which are derived from above kin (older sibs, younger aunts/uncles, older nieces/nephews):
  ## older unts/uncles, older cousins, younger cousins
  IC_older_aunt_uncle_f <- mothers_age_dist
  IC_older_aunt_uncle_m <- fathers_age_dist
  IC_older_cousins_f <- mothers_age_dist
  IC_older_cousins_m <- fathers_age_dist
  IC_younger_cousins_f <- mothers_age_dist
  IC_younger_cousins_m <- fathers_age_dist
  for(ic in 1 : (na-1)){
    X_older_aunt_uncle[,1] <- X_older_aunt_uncle[,1] + (IC_older_aunt_uncle_f[ic] + IC_older_aunt_uncle_m[ic])*prev_kin_older_sibs[,ic]
    X_older_cousins[,1] <- X_older_cousins[,1] + (IC_older_cousins_f[ic] + IC_older_cousins_m[ic])*prev_kin_older_niece_nephew[,ic]
    X_younger_cousins[,1] <- X_younger_cousins[,1] + (IC_younger_cousins_f[ic] + IC_younger_cousins_m[ic])*prev_kin_younger_niece_nephew[,ic]
  }
  
  ## Projections of older unts/uncles, older cousins, younger cousins
  for(i in 1: (na-1)){
    X_older_aunt_uncle[,i+1] <- Uproj %*% prev_kin_older_aunts_uncles[,i]
    X_older_cousins[,i+1] <- Uproj %*% prev_kin_older_cousins[,i] + Fproj %*% prev_kin_older_aunts_uncles[,i]
    X_younger_cousins[,i+1] <- Uproj %*% prev_kin_younger_cousins[,i] + Fproj %*% prev_kin_younger_aunts_uncles[,i]
  }
  
  return(list("Focal" = X_Focal,
              "d" = X_children,
              "gd" = X_grand_children,
              "ggd" = X_great_grand_children,
              "m" = X_parents,
              "gm" = X_grand_parents,
              "os" = X_older_sibs,
              "ys" = X_younger_sibs,
              "nos" = X_older_niece_nephew,
              "nys" = X_younger_niece_nephew,
              "oa" = X_older_aunt_uncle,
              "ya" = X_younger_aunts_uncles,
              "coa" = X_older_cousins,
              "cya" = X_younger_cousins,
              "ps" = population_age_stage_structure_next))
}

################## Create data frame output

## Use of "pipe" (don't understand the name, but hey)
`%>%` <- magrittr::`%>%`



#' Title Accumulated kin by each age of Focal, for each time period, and cohort of birth
#'
#' @param kin_matrix_lists list of lists of kin matrices: list( list(X_focal), list(X_parents), ... ). Outer list is length 14  = number of kin. Inner lists have lenght = timescale
#'                 so list(X_focal) = list(X_focal[year1],X_focal[year2],...,X_focal[yearlast])
#' @param kin_names list of characters. Corresponding to above lists: list("F","m",....)
#' @param years vector. The timescale on which we implement the kinship model.
#' @param start_year . First year of varying vital rates (e.g., if years = 1990:2000 then start_year = 1990)
#' @param na numeric. Number of ages.
#' @param ns numeric. Number of stages.
#' @param n_inc numeric. The size of the age/time increment (if abridged). NULL corresponds to 1 year intervals.
#' @param specific_kin character. names of kin we wish to analyse, e.g., list("os","ys"). If null returns all 14.
#'
#' @return A data frame which gives for each age of Focal at each year in the timescale, Focal's experienced number kin demarcated by stages (summed over all ages)
#' @export

create_cumsum_df <- function(kin_matrix_lists,
                             kin_names,
                             years,
                             start_year,
                             na,
                             ns,
                             n_inc,
                             specific_kin){
  df_year_list <- list()
  for(j in years){
    ii <- as.numeric(j) - start_year + 1
    df_list <- list()
    for(i in 1 : length(kin_names)){
      kin_member <- kin_names[[i]]
      kin_data <- kin_matrix_lists[[i]]
      kin_data <- kin_data[[ii]]
      df <- as.data.frame(as.matrix(kin_data))
      dims <- dim( kin_data)
      nr <- dims[1]
      nc <- dims[2]
      female_kin <- df[1:(nr/2), 1:nc]
      male_kin <- df[ (1+nr/2) : nr, 1:nc]
      female_kin$stage <- rep(seq(1, ns), na)
      male_kin$stage <- rep(seq(1, ns), na)
      female_kin$age <- rep(seq(0, (na-1), n_inc), each = ns)
      male_kin$age <- rep(seq(0, (na-1), n_inc), each = ns)
      female_kin$Sex <- "Female"
      male_kin$Sex <- "Male"
      both_kin <- rbind(female_kin, male_kin)
      both_kin <- both_kin %>% reshape2::melt(id = c("age","stage","Sex")) %>%
        dplyr::group_by(variable, stage, Sex) %>%
        dplyr::summarise(num = sum(value)) %>%
        dplyr::ungroup()
      both_kin <- both_kin %>% dplyr::transmute(age_focal = variable,
                                                stage_kin = as.factor(stage),
                                                count = num,
                                                sex_kin = Sex)
      both_kin$age_focal <- as.numeric(gsub("[^0-9.-]", "", both_kin$age_focal)) - 1
      df <- both_kin
      df$year <- j
      df$group <- kin_member
      df_list[[length(df_list)+1]] <- df
    }
    df_list <- do.call("rbind", df_list)
    df_year_list[[(1+length(df_year_list))]] <- df_list
  }
  df_year_list <- do.call("rbind", df_year_list)
  df_year_list <- df_year_list %>% dplyr::mutate(cohort = as.numeric(year) - as.numeric(age_focal),
                                                 cohort_factor = as.factor(cohort))
  if(specific_kin != FALSE){
    df_year_list <- df_year_list %>% dplyr::filter(group %in% specific_kin)
  }
  return(df_year_list)
}

#' Title joint age*stage distributions of kin by each age of Focal, for each time period, and cohort of birth
#'
#' @param kin_matrix_lists list of lists of kin matrices: list( list(X_focal), list(X_parents), ... ). Outer list is length 14  = number of kin. Inner lists have lenght = timescale
#'                 so list(X_focal) = list(X_focal[year1],X_focal[year2],...,X_focal[yearlast])
#' @param kin_names list of characters. Corresponding to above lists: list("F","m",....)
#' @param years vector. The timescale on which we implement the kinship model.
#' @param start_year . First year of varying vital rates (e.g., if years = 1990:2000 then start_year = 1990)
#' @param na numeric. Number of ages.
#' @param ns numeric. Number of stages.
#' @param n_inc numeric. The size of the age/time increment (if abridged). NULL corresponds to 1 year intervals.
#' @param specific_kin character. names of kin we wish to analyse, e.g., list("os","ys"). If null returns all 14.
#'
#' @return A data frame which gives for each age of Focal at each year in the timescale, the full age*stage dist of kin
#' @export
create_full_dists_df <- function(kin_matrix_lists,
                                 kin_names,
                                 years,
                                 start_year,
                                 na,
                                 ns,
                                 n_inc,
                                 specific_kin){
  df_year_list <- list()
  for(j in years){
    ii <- as.numeric(j) - start_year + 1
    df_list <- list()
    for(i in 1 : length(kin_names)){
      kin_member <- kin_names[[i]]
      kin_data <- kin_matrix_lists[[i]]
      kin_data <- kin_data[[ii]]
      df <- as.data.frame(as.matrix(kin_data))
      dims <- dim( kin_data)
      nr <- dims[1]
      nc <- dims[2]
      female_kin <- df[1:(nr/2), 1:nc]
      male_kin <- df[ (1+nr/2) : nr, 1:nc]
      female_kin$stage <- rep(seq(1, ns), na)
      male_kin$stage <- rep(seq(1, ns), na)
      female_kin$age <- rep(seq(0, (na-1), n_inc), each = ns)
      male_kin$age <- rep(seq(0, (na-1), n_inc), each = ns)
      female_kin$Sex <- "Female"
      male_kin$Sex <- "Male"
      both_kin <- rbind(female_kin, male_kin)
      both_kin <- both_kin %>% reshape2::melt(id = c("age","stage","Sex")) %>%
        dplyr::transmute(age_focal = variable,
                         age_kin = age,
                         stage_kin = as.factor(stage),
                         count = value,
                         sex_kin = Sex)
      both_kin$age_focal <- as.numeric(gsub("[^0-9.-]", "", both_kin$age_focal))-1
      df <- both_kin
      df$year <- j
      df$group <- kin_member
      df_list[[length(df_list)+1]] <- df
    }
    df_list <- do.call("rbind", df_list)
    df_year_list[[(1+length(df_year_list))]] <- df_list
  }
  df_year_list <- do.call("rbind", df_year_list)
  df_year_list <- df_year_list %>% dplyr::mutate(cohort = as.numeric(year) - as.numeric(age_focal),
                                                 cohort_factor = as.factor(cohort))
  if(specific_kin != FALSE){
    df_year_list <- df_year_list %>% dplyr::filter(group %in% specific_kin)
  }
  return(df_year_list)
}






