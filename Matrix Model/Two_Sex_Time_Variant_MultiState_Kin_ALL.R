

#' Estimate kin counts by age, stage, and sex, in a time variant framework

#' @description Implementation of combined formal demographic models: Caswell II,III,IV. 

#' @param U_list_females list with matrix entries: period-specific female survival probabilities. Age in rows and states in columns. 
#' @param U_list_males list with matrix entries: period-specific male survival probabilities. Age in rows and states in columns.
#' @param F_list_females list with matrix with elements: period-specific female fertility (age in rows and states in columns).
#' @param F_list_males list with matrix entries: period-specific male fertility (age in rows and states in columns). 
#' @param T_list_females list with matrix entries: period-specific female probability of transferring stage. 
#' @param T_list_males list with matrix entries: period-specific male probability of transferring stage
#' @param H_list list with matrix entries: redistribution of newborns across each stage to a specific age-class 
#' @param alpha numeric. ratio of males to females in population
#' @param parity logical. parity states imply age distribution of mothers re-scaled to not have parity 0 when Focal born. Default `TRUE`.
#' @param specific_kin vector. A vector of particular kin one wishes to obtain results for, e.g., c("m","d","oa"). Default is all kin types.
#' @param dist_output logical. Results as a data frame of accumulated kin by age of Focal if FALSE, and kin by their age*stage distribution by age of Focal if TRUE.
#' @param sex_Focal character. Female or Male as the user requests
#' @param stage_Focal Numeric in Natural number set {1,2,...,}. The stage which Focal is born into (e.g., 1 for parity 0)
#' @param nc numeric. The age/time-increment used in the discretisation of the continuum. 
#' @param time_series vector. The times at which we wish to count kin: start year = time_series[1], and end year = time_series[length.]
#' 
#' @return A data frame with focal´s age, related ages, stages, sexes, and types of kin for each time-period

kin_multi_stage_TV_2_sex_ <- function(U_list_females = NULL,
                                     U_list_males = NULL,
                                     F_list_females = NULL,
                                     F_list_males = NULL,
                                     T_list_females = NULL,
                                     T_list_males = NULL,
                                     H_list = NULL,
                                     alpha = 0.51, ## Sex ration -- UK value default
                                     parity = FALSE,
                                     specific_kin = FALSE,
                                     dist_output = FALSE,
                                     sex_Focal = "Female",
                                     stage_Focal = NULL,
                                     nc = NULL,
                                     time_series){
  
  no_years <- length(U_list_females)
  no_ages <- nrow(U_list_females[[1]])
  no_stages <- ncol(U_list_females[[1]])
  
  # Ensure inputs are lists of matrices
  if(!is.list(U_list_females) | !is.list(U_list_males)) stop("U's must be a list with time-series length. Each list entry should be an age*stage dimensional matrix")
  if(!is.list(F_list_females) | !is.list(F_list_males)) stop("F's must be a list with time-series length. Each list entry should be an age*stage dimensional matrix")
  if(!is.list(T_list_females) | !is.list(T_list_males)) stop("T's must be a list with time-series length. Each list entry should be an age*stage dimensional matrix")
  
  ### Define empty lists for the accumulated kin of Focals's life-course -- each list entry will reflect a time-period
  changing_pop_struct <- list()
  Focal_array <- list()
  mom_array <- list()
  gran_array <- list()
  daughter_array <- list()
  younger_sis_array <- list()
  grand_daughter_array <-list()
  great_grand_daughter_array <- list()
  older_sister_array <- list()
  younger_aunt_array <- list()
  older_aunt_array <- list()
  younger_niece_array <- list()
  older_niece_array <- list()
  younger_cousin_array <- list()
  older_cousin_array <- list()
  
  ### At each time-period we: 1) -- construct the time-variant projection matrices:
  ###                               U_tilde : transfers across stage and advances age
  ###                               F_tilde : makes newborns from stage/age; puts them to stage/age
  ###                         2) -- project Focal and kin using above projection matrices
  
  pb <- progress::progress_bar$new(
    format = "Timescale [:bar] :percent",
    total = no_years + 1, clear = FALSE, width = 60)
  tictoc::tic()
  for(year in 1:no_years){
    pb$tick()
    T_data_f <- T_list_females[[year]] ## For each year we have na number of Transfer matrices
    T_data_m <- T_list_males[[year]]   ## which give probabilities of age-dep movement from stage to stage
    T_f_list <- list()
    T_m_list <- list()
    F_f_list <- list()
    F_m_list <- list()
    U_f_list <- list()
    U_m_list <- list()
    H_list2 <- list()
    
    for(stage in 1:ns){
      Uf <- Matrix::Matrix(nrow = na, ncol = na, data = 0, sparse = TRUE)
      Matrix::diag(Uf[-1,-ncol(Uf)]) <- U_list_females[[year]][1:(na-1),stage]
      Uf[na,na] <- U_list_females[[year]][na,stage]
      Um <- Matrix::Matrix(nrow = na, ncol = na, data = 0, sparse = TRUE)
      Matrix::diag(Um[-1,-ncol(Um)]) <- U_list_males[[year]][1:(na-1),stage]
      Um[na,na] <- U_list_males[[year]][na,stage]
      U_f_list[[(1+length(U_f_list))]] <- Uf
      U_m_list[[(1+length(U_m_list))]] <- Um
      H_mat <- Matrix::Matrix(nrow = na, ncol = na, data = 0, sparse = TRUE)
      H_mat[1,] <- 1
      H_list2[[(1+length(H_list2))]] <- H_mat
    }
    for(age in 1:na){
      T_f <- T_data_f[[age]]
      T_m <- T_data_m[[age]]
      T_f_list[[(1+length(T_f_list))]] <- T_f
      T_m_list[[(1+length(T_m_list))]] <- T_m
      F_f <- Matrix::Matrix(nrow = ns, ncol = ns, data = 0, sparse = TRUE)
      F_m <- Matrix::Matrix(nrow = ns, ncol = ns, data = 0, sparse = TRUE)
      F_f[1,] <- F_list_females[[year]][age,]
      F_m[1,] <- F_list_males[[year]][age,]
      F_f_list[[(1+length(F_f_list))]] <- F_f
      F_m_list[[(1+length(F_m_list))]] <- F_m
    }
    ## create the appropriate block-diagonal matrices
    U_f_BDD <- block_diag_function(U_f_list) ## direct sum of female survivorship, independent over stage (ns diagonal blocks)
    U_m_BDD <- block_diag_function(U_m_list) ## direct sum of male survivorship, independent over stage (ns diagonal blocks)
    H_BDD <- block_diag_function(H_list2) ## direct sum of which age newborns enter, independent over stage (ns diagonal blocks)
    T_f_BDD <- block_diag_function(T_f_list) ## direct sum of female stage transitions, independent over age (na diagonal blocks)
    T_m_BDD <- block_diag_function(T_m_list) ## direct sum of male stage transitions, independent over age (na diagonal blocks)
    F_f_BDD <- block_diag_function(F_f_list) ## direct sum of female stage->stage reproductions, independent over age (na blocks)
    F_m_BDD <- block_diag_function(F_m_list) ## direct sum of male stage->stage reproductions, independent over age (na blocks)
    
    ## create the appropriate projection matrices
    U_tilde_females <- Matrix::t(K_perm_mat(ns, na)) %*%
      U_f_BDD %*%
      K_perm_mat(ns, na) %*%
      T_f_BDD
    
    ## create sex-specific age*stage projections
    U_tilde_males <- Matrix::t(K_perm_mat(ns, na)) %*%
      U_m_BDD %*%
      K_perm_mat(ns, na) %*%
      T_m_BDD
    
    F_tilde_females <- Matrix::t(K_perm_mat(ns, na)) %*%
      H_BDD %*%
      K_perm_mat(ns, na) %*%
      F_f_BDD
    
    F_tilde_males <- Matrix::t(K_perm_mat(ns, na)) %*%
      H_BDD %*%
      K_perm_mat(ns, na) %*%
      F_m_BDD
    
    ## if year == 1 we are at the boundary condition t=0 apply time-invariant kinship projections
    if(year == 1){ 
      ## Output of the static model 
      kin_out_1 <- all_kin_dy_(U_tilde_females, 
                              U_tilde_males , 
                              F_tilde_females, 
                              F_tilde_males, 
                              alpha, 
                              no_ages, 
                              no_stages, 
                              parity,
                              sex_Focal,
                              stage_Focal)
      ### Relative lists' first entries
      Focal_array[[(1+length(Focal_array))]] <- kin_out_1[[1]]
      daughter_array[[(1+length(daughter_array))]] <- kin_out_1[[2]]
      grand_daughter_array[[(1+length(grand_daughter_array))]] <- kin_out_1[[3]]
      great_grand_daughter_array[[(1+length(great_grand_daughter_array))]] <- kin_out_1[[4]]
      mom_array[[(1+length(mom_array))]] <- kin_out_1[[5]]
      gran_array[[(1+length(gran_array))]] <- kin_out_1[[6]]
      younger_sis_array[[( 1+length(younger_sis_array))]] <- kin_out_1[[8]]
      older_sister_array[[(1+length(older_sister_array))]] <- kin_out_1[[7]]
      younger_aunt_array[[(1+length(younger_aunt_array))]] <- kin_out_1[[12]]
      older_aunt_array[[(1+length(older_aunt_array))]] <- kin_out_1[[11]]
      younger_niece_array[[(1+length(younger_niece_array))]] <- kin_out_1[[10]]
      older_niece_array[[(1+length(older_niece_array))]] <- kin_out_1[[9]]
      younger_cousin_array[[(1+length(younger_cousin_array))]] <- kin_out_1[[14]]
      older_cousin_array[[(1+length(older_cousin_array))]] <- kin_out_1[[13]]
      changing_pop_struct[[(1+length(changing_pop_struct))]] <- kin_out_1[[15]]
      
    }
    updating_Focal <- Focal_array[[year]]
    updating_daughter <- daughter_array[[year]]
    updating_grand_daughter <- grand_daughter_array[[year]]
    updating_great_grand_daughter <- great_grand_daughter_array[[year]]
    updating_mom <- mom_array[[year]]
    updating_gran <- gran_array[[year]]
    updating_younger_sis <- younger_sis_array[[year]]
    updating_older_sis <- older_sister_array[[year]]
    updating_youner_aunt <- younger_aunt_array[[year]]
    updating_older_aunt <- older_aunt_array[[year]]
    updating_younger_niece <- younger_niece_array[[year]]
    updating_older_niece <- older_niece_array[[year]]
    updating_younger_cousin <- younger_cousin_array[[year]]
    updating_older_cousin <- older_cousin_array[[year]]
    updating_pop_struct <- changing_pop_struct[[year]]
    
    ## Output of the time-variant model 
    kin_out <- all_kin_dy_TV_(U_tilde_females, 
                             U_tilde_males, 
                             F_tilde_females, 
                             F_tilde_males, 
                             alpha, 
                             no_ages, 
                             no_stages,  
                             parity,
                             sex_Focal,
                             stage_Focal,
                             updating_Focal,
                             updating_daughter, 
                             updating_grand_daughter,
                             updating_great_grand_daughter,
                             updating_mom,
                             updating_gran,
                             updating_older_sis,
                             updating_younger_sis,
                             updating_older_niece,
                             updating_younger_niece,
                             updating_older_aunt,
                             updating_youner_aunt,
                             updating_older_cousin,
                             updating_younger_cousin, 
                             updating_pop_struct)
    ## Relative lists entries correspond to timescale periods (each entry an kin age*stage*2 by Focal age matrix)
    Focal_array[[(1+length(Focal_array))]] <- kin_out[[1]]
    daughter_array[[(1+length(daughter_array))]] <- kin_out[[2]]
    grand_daughter_array[[(1+length(grand_daughter_array))]] <- kin_out[[3]]
    great_grand_daughter_array[[(1+length(great_grand_daughter_array))]] <- kin_out[[4]]
    mom_array[[(1+length(mom_array))]] <- kin_out[[5]]
    gran_array[[(1+length(gran_array))]] <- kin_out[[6]]
    younger_sis_array[[(1+length(younger_sis_array))]] <- kin_out[[8]]
    older_sister_array[[(1+length(older_sister_array))]] <- kin_out[[7]]
    younger_aunt_array[[(1+length(younger_aunt_array))]] <- kin_out[[12]]
    older_aunt_array[[(1+length(older_aunt_array))]] <- kin_out[[11]]
    younger_niece_array[[(1+length(younger_niece_array))]] <- kin_out[[10]]
    older_niece_array[[(1+length(older_niece_array))]] <- kin_out[[9]]
    younger_cousin_array[[(1+length(younger_cousin_array))]] <- kin_out[[14]]
    older_cousin_array[[(1+length(older_cousin_array))]] <- kin_out[[13]]
    changing_pop_struct[[(1+length(changing_pop_struct))]] <- kin_out[[15]]
  }
  tictoc::toc()
  ## create a list of output kin -- each element a time-period specific list of matrices
  relative_data <- list(Focal_array,
                        daughter_array,
                        grand_daughter_array,
                        great_grand_daughter_array,
                        mom_array,
                        gran_array,
                        younger_sis_array,
                        older_sister_array,
                        younger_aunt_array,
                        older_aunt_array,
                        younger_niece_array,
                        older_niece_array,
                        younger_cousin_array,
                        older_cousin_array)
  ## label the kin names to match DemoKin: 
  ## "coa", "cya", "d", "gd", "ggd", "ggm", "gm", "m", "nos", "nys", "oa", "ya", "os", "ys")
  relative_names <- list("Focal",
                         "d",
                         "gd", 
                         "ggd",
                         "m", 
                         "gm" , 
                         "ys", 
                         "os", 
                         "ya" , 
                         "oa",
                         "nys",
                         "nos",
                         "cya",
                         "coa")
  ## create a nice data frame output
  if(!dist_output){
    kin_out <- create_cumsum_df_(relative_data,
                                relative_names, 
                                time_series[1]:time_series[length(time_series)], 
                                time_series[1], 
                                no_ages, 
                                no_stages, 
                                nc,
                                specific_kin)}
  else{
    kin_out <- create_full_dists_df_(relative_data,
                                    relative_names, 
                                    time_series[1]:time_series[length(time_series)], 
                                    time_series[1], 
                                    no_ages, 
                                    no_stages, 
                                    nc,
                                    specific_kin)}
  
  return(kin_out)
}





source(here::here("Matrix Model", "Functions_required.R" ))

##################### Not time varying (boundary time = 0 kin output)

all_kin_dy_ <- function(Uf, 
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
    stage_Focal <- 1
    
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
    entry <- 1 + (stage_Focal-1)*na
    IC_Focal[entry] <- 1}
  else{
    entry <- n + 1 + (stage_Focal-1)*na
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

all_kin_dy_TV_ <- function(Uf, 
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
  X_Focal  <- Matrix::Matrix(nrow = (2*n), ncol = na, data = 0, sparse = TRUE)
  IC_Focal <- rep(0, 2*n)
  if(sex_Focal == "Female"){
    entry <- 1 + (stage_Focal-1)*na
    IC_Focal[entry] <- 1}
  else{
    entry <- n + 1 + (stage_Focal-1)*na
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



################## Create data frame output

################### Turn the above output into nice data frames.

## Use of "pipe" (don't understand the name, but hey)

`%>%` <- magrittr::`%>%`

##### When running the model over a time-series, i use list comprehension: 
##### each additional list entry will give Focal's kinship network (over all ages of Focal's life) for the next time-period
#### I transform these into period-specific data frames of 1) accumulated kin by age of Focal and 2) age*stage dists of kin by age of Focal

## Arguments in the below functions: 
##        dat_list = a list of lists of matrices: list( list(X_foc), list(X_child), ... ) -- each list entry with year/age/sex specific kin distributions
##        list_dist = a character list of which kin we want to analyse
##        years = sequence of years in the time series
##        start_year = when we begin the time-series
##        na = number of ages
##        ns = number of stages
##        nc = age increment (i.e., 1-year age-classes, or 5-year age-classes)

################################################## 1) Accumulated kin by age of Focal 
create_cumsum_df_ <- function(dat_list, 
                             list_dists,
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
    for(i in 1 : length(list_dists)){
      kin_member <- list_dists[[i]]
      kin_data <- dat_list[[i]]
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
      both_kin <- both_kin %>% dplyr::transmute(Age_Focal = variable, 
                                                Stage_Kin = as.factor(stage), 
                                                pred_no_kin = num,
                                                Sex = Sex)
      both_kin$Age_Focal <- as.numeric(gsub("[^0-9.-]", "", both_kin$Age_Focal)) - 1
      df <- both_kin
      df$year <- j
      df$group <- kin_member
      df_list[[length(df_list)+1]] <- df
    }
    df_list <- do.call("rbind", df_list)
    df_year_list[[(1+length(df_year_list))]] <- df_list
  }
  df_year_list <- do.call("rbind", df_year_list)
  df_year_list <- df_year_list %>% dplyr::mutate(cohort = as.numeric(year) - as.numeric(Age_Focal),
                                                 cohort_factor = as.factor(cohort))
  if(specific_kin != FALSE){
    df_year_list <- df_year_list %>% dplyr::filter(group %in% specific_kin)
  }
  return(df_year_list)
}
################################################## 2) Full age*stage distributions by age of Focal
create_full_dists_df_ <- function(dat_list, 
                                 list_dists,
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
    for(i in 1 : length(list_dists)){
      kin_member <- list_dists[[i]]
      kin_data <- dat_list[[i]]
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
        dplyr::transmute(Age_Focal = variable, 
                         Age_Kin = age,
                         Stage_Kin = as.factor(stage), 
                         pred_no_kin = value,
                         Sex = Sex)
      both_kin$Age_Focal <- as.numeric(gsub("[^0-9.-]", "", both_kin$Age_Focal))-1
      df <- both_kin
      df$year <- j
      df$group <- kin_member
      df_list[[length(df_list)+1]] <- df
    }
    df_list <- do.call("rbind", df_list)
    df_year_list[[(1+length(df_year_list))]] <- df_list
  }
  df_year_list <- do.call("rbind", df_year_list)
  df_year_list <- df_year_list %>% dplyr::mutate(cohort = as.numeric(year) - as.numeric(Age_Focal),
                                                 cohort_factor = as.factor(cohort))
  if(specific_kin != FALSE){
    df_year_list <- df_year_list %>% dplyr::filter(group %in% specific_kin)
  }
  return(df_year_list)
}





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
get_G <- function(U, no_ages, no_stages){
  sig <- Matrix::t(rep(1,no_ages*no_stages)) %*% U
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
marg_stage_dist <- function(no_ages, no_stages, full_dist){
  return(kronecker( Matrix::t(rep(1, no_ages)) , Matrix::diag(no_stages) ) %*% full_dist)}

# The marginal age dist (i.e., summing over all stages)
marg_age_dist <- function(no_ages, no_stages, full_dist){
  return(kronecker( Matrix::diag(no_ages) , Matrix::t(rep(1, no_stages)) ) %*% full_dist)
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







