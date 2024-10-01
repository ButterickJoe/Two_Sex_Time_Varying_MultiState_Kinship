

#' Estimate kin counts by age and stage over two sexes, in a time variant framework

#' @description Implementation of combined formal demographic models of Caswell II,III,IV. 

#' @param U_list_females -- list of matrices with elements: female probability of advancing age-class. Age in rows and states in columns. Each list entry ~ time-period.
#' @param U_list_males -- list of matrices with elements: male probability of advancing age-class. Age in rows and states in columns. Each list entry ~ time-period.
#' @param F_list_females list of matrices with elements: female state-specific fertility (age in rows and states in columns). Each list entry ~ time-period.
#' @param F_list_males list of matrices with elements: male state-specific fertility (age in rows and states in columns). Each list entry ~ time-period.
#' @param T_list_females list of matrices with elements: female probability of moving from one stage to another. Each list entry ~ time-period
#' @param T_list_males list of matrices with elements: male probability of moving from one stage to another. Each list entry ~ time-period
#' @param H_list redistributes newborns across each stage to a specific age-class 
#' @param output_kin character. kin to return. For example "m" for mother, "d" for daughter. See the `vignette` for all kin types.
#' @param alpha numeric. Female portion at birth.
#' @param parity logical. parity states imply age distribution of mothers re-scaled to not have parity 0 when Focal born. Default `TRUE`.
#' @param list_output logical. Results as a list. Default `FALSE`.
#' @param sex_Focal character. Female or Male as the user requests
#' @param stage_Focal Numeric in Natural number set {1,2,...,}. The stage which Focal is born into (e.g., 1 for parity 0)
#' @param time_series vector. The times at which we wish to count kin: start year = time_series[1], and end year = time_series[length.]
#' 
#' @return A data frame with focal´s age, related ages, stages, sexes, and types of kin for each time-period
 
## Import required functions and matrix operations
source(here::here("Matrix Model", "Functions_required.R" ))
source(here::here("Matrix Model", "Kin_projections.R" ))
source(here::here("Matrix Model", "Creating_data_frames.R"))

kin_multi_stage_TV_2_sex <- function(U_list_females = NULL,
                                     U_list_males = NULL,
                                     F_list_females = NULL,
                                     F_list_males = NULL,
                                     T_list_females = NULL,
                                     T_list_males = NULL,
                                     H_list = NULL,
                                     alpha = 0.51, ## Sex ration -- UK value default
                                     output_kin = NULL,
                                     parity = TRUE,
                                     dist_output = FALSE,
                                     sex_Focal = "Female",
                                     stage_Focal = NULL,
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
    format = "Running over input years [:bar] :percent",
    total = no_years + 1, clear = FALSE, width = 60)
  
  foreach(year = 1:no_years)%do%{
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
    
    foreach(stage = 1:ns)%do%{
      Uf <- as(matrix(0, nrow = na, ncol = na),"sparseMatrix")
      diag(Uf[-1,-ncol(Uf)]) <- U_list_females[[year]][1:(na-1),stage]
      Uf[na,na] <- U_list_females[[year]][na,stage]
      Um <- as(matrix(0, nrow = na, ncol = na),"sparseMatrix")
      diag(Um[-1,-ncol(Um)]) <- U_list_males[[year]][1:(na-1),stage]
      Um[na,na] <- U_list_males[[year]][na,stage]
      U_f_list[[(1+length(U_f_list))]] <- Uf
      U_m_list[[(1+length(U_m_list))]] <- Um
      H_mat <- as(matrix(0, nrow = na, ncol = na), "sparseMatrix")
      H_mat[1,] <- 1
      H_list2[[(1+length(H_list2))]] <- H_mat
    }
    foreach(age = 1:na)%do%{
      T_f <- T_data_f[[age]]
      T_m <- T_data_m[[age]]
      T_f_list[[(1+length(T_f_list))]] <- T_f
      T_m_list[[(1+length(T_m_list))]] <- T_m
      F_f <- as(matrix(0, nrow = ns, ncol = ns), "sparseMatrix")
      F_m <- as(matrix(0, nrow = ns, ncol = ns), "sparseMatrix")
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
    U_tilde_females <- t(K_perm_mat(ns, na))%*%U_f_BDD%*%K_perm_mat(ns, na)%*%T_f_BDD ## create sex-specific age*stage projections
    U_tilde_males <- t(K_perm_mat(ns, na))%*%U_m_BDD%*%K_perm_mat(ns, na)%*%T_m_BDD
    F_tilde_females <- t(K_perm_mat(ns, na))%*%H_BDD%*%K_perm_mat(ns, na)%*%F_f_BDD
    F_tilde_males <- t(K_perm_mat(ns, na))%*%H_BDD%*%K_perm_mat(ns, na)%*%F_m_BDD
    
    ## if year == 1 we are at the boundary condition t=0 --> time-invariant kinship projections
    if(year == 1){ 
      ## Output of the static model 
      kin_out_1 <- all_kin_dy(U_tilde_females, 
                              U_tilde_males , 
                              F_tilde_females, 
                              F_tilde_males, 
                              alpha, 
                              no_ages, 
                              no_stages, 
                              parity = TRUE,
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
    kin_out <- all_kin_dy_TV(U_tilde_females, 
                             U_tilde_males, 
                             F_tilde_females, 
                             F_tilde_males, 
                             alpha, 
                             no_ages, 
                             no_stages,  
                             parity = TRUE,
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
                             updating_younger_cousin, updating_pop_struct)
    ## Relative lists entries for t>0
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
  ## label the kin names
  relative_names <- list("Focal",
                         "offspring",
                         "grand offspring", 
                         "great grand offspring",
                         "parents", 
                         "grand parents" , 
                         "younger siblings", 
                         "older siblings", 
                         "younger aunt/unlces" , 
                         "older aunt/unlce",
                         "younger niece/nephews",
                         "older niece/nephews",
                         "younger cousin",
                         "older cousin")
  ## create a nice data frame output
  if(!dist_output){
    kin_out <- create_cumsum_TV(relative_data,
                                relative_names, 
                                time_series[1]:time_series[length(time_series)], 
                                time_series[1], 
                                no_ages, 
                                no_stages, 
                                1)}
  else{
  kin_out <- create_full_dists_TV(relative_data,
                              relative_names, 
                              time_series[1]:time_series[length(time_series)], 
                              time_series[1], 
                              no_ages, 
                              no_stages, 
                              1)}
  
  return(kin_out)
}
  
  