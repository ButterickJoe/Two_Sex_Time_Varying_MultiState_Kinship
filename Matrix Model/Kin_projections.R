

source(here::here("Matrix Model", "Functions_required.R" ))

##################### Not time varying (boundary time = 0 kin output)

all_kin_dy <- function(Uf, 
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
  ###############################################################################################################################
  
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
  X_Focal[,1] <- IC_Focal
  ## Project
  for(i in 1 : (na -1) ){
    X_Focal[,i+1] <- G_tilde%*%X_Focal[,i] 
  }
  
  ######### Category 1) Project Descendants #########################
  
  ################################################### Children
  X_children <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_children <- rep(0, 2*n)
  X_children[,1] <- IC_children
  foreach(i = 1 : (na-1))%do%{
    X_children[,i+1] <- Uproj%*%X_children[,i] + Fproj%*%X_Focal[,i] # e_vector(n+i,2*n) if focal is male
  }
  
  ################################################## Grand Children
  X_grand_children <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  recruit_grand_children <- X_children ## recruitment from children's reproduction
  IC_grand_children <- rep(0, 2*n) ## No grand-children at birth IC = (0,0,....,0) (2n length)
  X_grand_children[,1] <- IC_grand_children
  ## Project
  foreach(i = 1 : (na-1))%do%{
    X_grand_children[,i+1] <- Uproj%*%X_grand_children[,i] + Fproj%*%recruit_grand_children[,i]
  }
  
  ################################################ Great grand children 
  X_great_grand_children <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  recruit_great_grand_children <- X_grand_children ## recruitment from children's reproduction
  IC_great_grand_children <- rep(0, 2*n) ## No grand-children at birth IC = (0,0,....,0) (2n length)
  X_great_grand_children[,1] <- IC_great_grand_children
  ## Project
  foreach(i = 1 : (na-1))%do%{
    X_great_grand_children[,i+1] <- Uproj%*%X_great_grand_children[,i] + Fproj%*%recruit_great_grand_children[,i]
  }
  
  ######### Category 2) Project Ancestors #########################
  
  #################################################### Parents
  X_parents <- as(matrix(0, nrow = 2*n, ncol = na, byrow = T),"sparseMatrix")
  IC_parents <- mothers_age_stage
  X_parents[, 1] <- IC_parents
  ## Project
  foreach(i = 1 : (na-1))%do%{
    X_parents[, i+1] <- Uproj%*%X_parents[, i] 
  }
  ################################################### Grand parents
  X_grand_parents <- as( matrix(0, nrow = 2*n, ncol = na, byrow = T) ,"sparseMatrix")
  ## Initial conditions (female and male grand parents)
  IC_gps <- rep(0, 2*n)
  IC_f_grand_pars <- mothers_age_dist
  IC_m_grand_pars <- fathers_age_dist
  foreach(ic = 1 : (na))%do%{
    IC_gps <- IC_gps + (IC_f_grand_pars[ic] + IC_m_grand_pars[ic])*X_parents[,ic] ## IC the sum of parents of Focal's parents, given their age
  }
  X_grand_parents[, 1] <- IC_gps
  ## Project
  foreach(i = 1 : (na-1))%do%{
    X_grand_parents[, i+1] <- Uproj%*%X_grand_parents[, i] 
  }
  
  #### Category 3). Siblings ########################
  
  ################################################# Older sibs
  X_older_sibs  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  ## Initial conditions (female and male older siblings)
  IC_os <- rep(0, 2*n)
  IC_older_sibs_f <- mothers_age_dist
  foreach(ic = 1 : (na))%do%{
    IC_os <- IC_os + IC_older_sibs_f[ic]*X_children[,ic]
  }
  X_older_sibs[, 1] <- IC_os
  ## Project
  foreach(i = 1 : (na-1))%do%{
    X_older_sibs[,i+1] <- Uproj%*%X_older_sibs[,i] 
  }
  ################################################ Younger siblings
  X_younger_sibs  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  ## Initial conditions (female and male younger siblings)
  IC_younger_sibs <- rep(0, n*2) ## No younger sisters at birth 
  X_younger_sibs[,1] <- IC_younger_sibs
  ## Recruitment
  recruit_younger_sibs <- X_parents
  foreach(i = 1 : (na-1))%do%{
    X_younger_sibs[,i+1] <- Uproj%*%X_younger_sibs[,i] + Fprojstar%*%recruit_younger_sibs[,i]
  }
  
  #### Category 4). Nieces and nephews ##########
  
  ################################################ Older Nieces and Nephews
  X_older_niece_nephew  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  ## Initial conditions (female and male older niece/nephew) expected number of grandchildren of mother * age specific fertility
  IC_onn <- rep(0, n*2) 
  IC_older_niece_nephew_f <- mothers_age_dist
  foreach(ic = 1 : (na-1))%do%{
    IC_onn <- IC_onn + IC_older_niece_nephew_f[ic]*X_grand_children[,ic]
  }
  X_older_niece_nephew[,1] <- IC_onn
  ## Recruitment
  recruit_older_niece_nephew <- X_older_sibs
  foreach(i = 1 : (na-1))%do%{
    X_older_niece_nephew[,i+1] <- Uproj%*%X_older_niece_nephew[,i] + Fproj%*%recruit_older_niece_nephew[,i]
  }
  ############################################### Younger Nieces and Nephews
  X_younger_niece_nephew  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  ## recruitment from younger siblings
  recruit_younger_niece_nephew <- X_younger_sibs
  IC_younger_niece_nephew <- rep(0, n*2) ## At birth, Nieces/Nephews produced by nonexistent younger sibs impossible
  X_younger_niece_nephew[,1] <- IC_younger_niece_nephew
  foreach(i = 1 : (na-1))%do%{
    X_younger_niece_nephew[,i+1] <- Uproj%*%X_younger_niece_nephew[,i] + Fproj%*%recruit_younger_niece_nephew[,i]
  }
  
  #### Category 5). Aunts and uncles ##############
  
  ################################################ Older Aunts and Uncles
  X_older_aunt_uncle  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  ## IC
  IC_au <- rep(0, n*2) ## At birth, older sibs of mother and father over all their possible ages:
  pi_older_aunt_uncle_f <- mothers_age_dist
  pi_older_aunt_uncle_m <- fathers_age_dist
  foreach(ic = 1 : (na-1))%do%{
    IC_au <- IC_au + (pi_older_aunt_uncle_f[ic] + pi_older_aunt_uncle_m[ic])*X_older_sibs[,ic]
  }
  X_older_aunt_uncle[,1] <- IC_au
  foreach(i = 1 : (na-1))%do%{
    X_older_aunt_uncle[,i+1] <- Uproj%*%X_older_aunt_uncle[,i] 
  }
  ################################################ Younger aunts and uncles
  X_younger_aunts_uncles <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  ## IC
  IC_yau <- rep(0, n*2) ## At birth, younger sibs of mother and father over all ages
  pi_younger_aunts_uncles_f <- mothers_age_dist
  pi_younger_aunts_uncles_m <- fathers_age_dist
  foreach(ic = 1 : (na-1))%do%{
    IC_yau <- IC_yau + (pi_younger_aunts_uncles_f[ic] + pi_younger_aunts_uncles_m[ic])*X_younger_sibs[,ic]
  }
  ## recruitment comes from X_grand_parents here
  X_younger_aunts_uncles[,1] <- IC_yau
  foreach(i = 1 : (na-1))%do%{
    X_younger_aunts_uncles[,i+1] <- Uproj%*%X_younger_aunts_uncles[,i] + Fprojstar%*%X_grand_parents[,i]
  }
  
  #### Last Category I consider (but please add more) is Cousins ###########
  
  ################################################# Older cousins
  X_older_cousins  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  ## Initial conditions: 
  IC_oc <- rep(0, n*2) ## At birth, older cousins come from mother and father's older niece/nephews
  IC_older_cousins_f <- mothers_age_dist
  IC_older_cousins_m <- fathers_age_dist
  foreach(ic = 1 : (na-1))%do%{
    IC_oc <- IC_oc + (IC_older_cousins_f[ic] + IC_older_cousins_m[ic])*X_older_niece_nephew[,ic]
  }
  X_older_cousins[,1] <- IC_oc
  ## Recruitment from older_aunts_uncles
  foreach(i = 1 : (na-1))%do%{
    X_older_cousins[,i+1] <- Uproj%*%X_older_cousins[,i] + Fproj%*%X_older_aunt_uncle[,i]
  }
  
  ################################################## Younger cousins
  X_younger_cousins  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  ## Initial conditions: 
  IC_yc <- rep(0, n*2) ## At birth, younger cousins come from mother/father's younger niece/nephew
  IC_younger_cousins_f <- mothers_age_dist
  IC_younger_cousins_m <- fathers_age_dist
  foreach(ic = 1 : (na-1))%do%{
    IC_yc <- IC_yc + (IC_younger_cousins_f[ic] + IC_younger_cousins_m[ic])*X_younger_niece_nephew[,ic]
  }
  X_younger_cousins[,1] <- IC_yc
  ## recruitment from younger_aunts_uncle
  foreach(i = 1 : (na-1))%do%{
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
##########################################################################################################
########################################## Now a time varying kin output -- all annotations the same as above

all_kin_dy_TV <- function(Uf, 
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
  ########################################### Matrix projections for kinship
  ##############################################################################################################
  
  ### Focal 
  f_t <- get_G(Uf, na, ns)
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
  X_Focal[,1] <- IC_Focal
  for(i in 1 : (na -1) ){
    X_Focal[,i+1] <- G_tilde%*%previous_kin_Focal[,i] 
  }
  # children
  X_children <-as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_children <- rep(0, n*2)
  X_children[,1] <- IC_children
  foreach(i = 1 : (na-1))%do%{
    X_children[,i+1] <- Uproj%*%prev_kin_children[,i] + Fproj%*%previous_kin_Focal[,i] 
  }
  ## Grand Children
  X_grand_children <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_grand_children <- rep(0, n*2)
  X_grand_children[,1] <- IC_grand_children
  foreach(i = 1 : (na-1))%do%{
    X_grand_children[,i+1] <- Uproj%*%prev_kin_grandchildren[,i] + Fproj%*%prev_kin_children[,i]
  }
  ## Great grand children 
  X_great_grand_children <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_great_grand_children <- rep(0, n*2) 
  X_great_grand_children[,1] <- IC_great_grand_children
  foreach(i = 1 : (na-1))%do%{
    X_great_grand_children[,i+1] <- Uproj%*%prev_kin_greatgrandchildren[,i] + Fproj%*%prev_kin_grandchildren[,i]
  }
  ## Ancestors
  ############################################ Parents
  X_parents <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_parents <- mothers_age_stage
  X_parents[, 1] <- IC_parents
  foreach(i = 1 : (na-1))%do%{
    X_parents[, i+1] <- Uproj%*%prev_kin_parents[, i] 
  }
  ########################################### Grand parents
  X_grand_parents <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_f_grand_parents <- mothers_age_dist
  IC_m_grand_parents <- fathers_age_dist
  IC_gps <- rep(0, 2*n)
  foreach(i = 1: (na-1))%do%{
    IC_gps <- IC_gps + (IC_f_grand_parents[i] + IC_m_grand_parents[i])*prev_kin_parents[,i]
  }
  X_grand_parents[, 1] <- IC_gps
  foreach(i = 1 : (na-1))%do%{
    X_grand_parents[, i+1] <- Uproj%*%prev_kin_grand_parents[, i] 
  }
  ############################################## Older sibs
  X_older_sibs  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_older_sibs_f <- mothers_age_dist
  IC_os <- rep(0, 2*n)
  foreach(ic = 1 : (na-1))%do%{
    IC_os <- IC_os + IC_older_sibs_f[ic]*prev_kin_children[,ic]
  }
  X_older_sibs[, 1] <- IC_os
  foreach(i = 1 : (na-1))%do%{
    X_older_sibs[,i+1] <- Uproj%*%prev_kin_older_sibs[,i] 
  }
  ############################################### Younger sibs
  X_younger_sibs  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_younger_sibs <- rep(0, n*2) 
  X_younger_sibs[,1] <- IC_younger_sibs
  foreach(i = 1 : (na-1))%do%{
    X_younger_sibs[,i+1] <- Uproj%*%prev_kin_younger_sibs[,i] + Fprojstar%*%prev_kin_parents[,i]
  }
  ############################################## Older Nieces and Nephews
  X_older_niece_nephew  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_onn <- rep(0, n*2) 
  IC_older_niece_nephew_f <- mothers_age_dist
  foreach(ic = 1 : (na-1))%do%{
    IC_onn <- IC_onn + IC_older_niece_nephew_f[ic]*prev_kin_grandchildren[,ic]
  }
  X_older_niece_nephew[,1] <- IC_onn
  foreach(i = 1 : (na-1))%do%{
    X_older_niece_nephew[,i+1] <- Uproj%*%prev_kin_older_niece_nephew[,i] + Fproj%*%prev_kin_older_sibs[,i]
  }
  ########################################### Younger Nieces and Nephews
  X_younger_niece_nephew  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_younger_niece_nephew <- rep(0, n*2) 
  X_younger_niece_nephew[,1] <- IC_younger_niece_nephew
  foreach(i = 1 : (na-1))%do%{
    X_younger_niece_nephew[,i+1] <- Uproj%*%prev_kin_younger_niece_nephew[,i] + Fproj%*%prev_kin_younger_sibs[,i]
  }
  ########################################## Older Aunts and Uncles
  X_older_aunt_uncle  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_au <- rep(0, n*2) 
  pi_older_aunt_uncle_f <- mothers_age_dist
  pi_older_aunt_uncle_m <- fathers_age_dist
  foreach(ic = 1 : (na-1))%do%{
    IC_au <- IC_au + (pi_older_aunt_uncle_f[ic] + pi_older_aunt_uncle_m[ic])*prev_kin_older_sibs[,ic]
  }
  X_older_aunt_uncle[,1] <- IC_au
  foreach(i = 1 : (na-1))%do%{
    X_older_aunt_uncle[,i+1] <- Uproj%*%prev_kin_older_aunts_uncles[,i] 
  }
  ######################################### Younger aunts and uncles
  X_younger_aunts_uncles <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_yau <- rep(0, n*2) 
  pi_younger_aunts_uncles_f <- mothers_age_dist
  pi_younger_aunts_uncles_m <- fathers_age_dist
  foreach(ic = 1 : (na-1))%do%{
    IC_yau <- IC_yau + (pi_younger_aunts_uncles_f[ic] + pi_younger_aunts_uncles_m[ic])*prev_kin_younger_sibs[,ic]
  }
  X_younger_aunts_uncles[,1] <- IC_yau
  foreach(i = 1 : (na - 1) )%do%{
    X_younger_aunts_uncles[,i+1] <- Uproj%*%prev_kin_younger_aunts_uncles[,i] + Fprojstar%*%prev_kin_grand_parents[,i]
  }
  ######################################## Older cousins
  X_older_cousins  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_oc <- rep(0, n*2) 
  IC_older_cousins_f <- mothers_age_dist
  IC_older_cousins_m <- fathers_age_dist
  foreach(ic = 1 : (na))%do%{
    IC_oc <- IC_oc + (IC_older_cousins_f[ic] + IC_older_cousins_m[ic])*prev_kin_older_niece_nephew[,ic]
  }
  X_older_cousins[,1] <- IC_oc
  foreach(i = 1 : (na - 1) )%do%{
    X_older_cousins[,i+1] <- Uproj%*%prev_kin_older_cousins[,i] + Fproj%*%prev_kin_older_aunts_uncles[,i]
  }
  ########################################### Younger cousins
  X_younger_cousins  <- as(matrix(0, nrow = 2*n, ncol = na, byrow = TRUE),"sparseMatrix")
  IC_yc <- rep(0, n*2) 
  IC_younger_cousins_f <- mothers_age_dist
  IC_younger_cousins_m <- fathers_age_dist
  foreach(ic = 1 : (na-1))%do%{
    IC_yc <- IC_yc + (IC_younger_cousins_f[ic] + IC_younger_cousins_m[ic])*prev_kin_younger_niece_nephew[,ic]
  }
  X_younger_cousins[,1] <- IC_yc
  foreach(i = 1 : (na - 1) )%do%{
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


################### Turn the above output into nice data frames.

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

################################################## 1) Accumulated kin by age of Focal 
create_cumsum_TV <- function(dat_list, 
                             list_dists,
                             years, 
                             start_year, 
                             na, 
                             ns, 
                             n_inc){
  df_year_list <- list()
  foreach(j = years)%do%{
    ii <- as.numeric(j) - start_year + 1
    df_list <- list()
    foreach( i = 1 : length(list_dists))%do%{
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
      both_kin <- both_kin%>%melt(id = c("age","stage","Sex"))%>%
        group_by(variable,stage,Sex)%>%
        summarise(num = sum(value))%>%
        ungroup()
      both_kin <- both_kin%>%transmute(Age_Focal = variable, 
                                       Stage = as.factor(stage), 
                                       pred_no_kin = num,
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
  df_year_list <- df_year_list%>% mutate(cohort = as.numeric(year) - as.numeric(Age_Focal),
                                         cohort_factor = as.factor(cohort))
  return(df_year_list)
}
################################################## 2) Full age*stage distributions by age of Focal
create_full_dists_TV <- function(dat_list, 
                                 list_dists,
                                 years, 
                                 start_year, 
                                 na, 
                                 ns, 
                                 n_inc){
  df_year_list <- list()
  foreach(j = years)%do%{
    ii <- as.numeric(j) - start_year + 1
    df_list <- list()
    foreach( i = 1 : length(list_dists))%do%{
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
      both_kin <- both_kin%>%melt(id = c("age","stage","Sex"))%>%
        transmute(Age_Focal = variable, 
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
  df_year_list <- df_year_list%>% mutate(cohort = as.numeric(year) - as.numeric(Age_Focal),
                                         cohort_factor = as.factor(cohort))
  return(df_year_list)
}






