


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
create_cumsum_df <- function(dat_list, 
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
create_full_dists_df <- function(dat_list, 
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






