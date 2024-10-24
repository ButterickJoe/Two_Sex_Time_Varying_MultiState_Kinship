


source(here::here("Matrix Model", "matrix_operations.R" ))
source(here::here("Matrix Model", "kin_projections_and_df_construction.R" ))

#' Estimate kin counts by age, stage, and sex, in a time variant framework

#' @description Implementation of combined formal demographic models: Caswell II,III,IV.

#' @param U_list_females list with matrix entries: period-specific female survival probabilities. Age in rows and states in columns.
#' @param U_list_males list with matrix entries: period-specific male survival probabilities. Age in rows and states in columns.
#' @param F_list_females list with matrix with elements: period-specific female fertility (age in rows and states in columns).
#' @param F_list_males list with matrix entries: period-specific male fertility (age in rows and states in columns).
#' @param T_list_females list of lists with matrix entries: each outer list entry is period-specific, and composed of
#'                     a list of stochastic matrices which describe age-specific female probabilities of transferring stage
#' @param T_list_males list of lists with matrix entries: each outer list entry is period-specific, and composed of
#'                     a list of stochastic matrices which describe age-specific male probabilities of transferring stage
#' @param H_list list with matrix entries: redistribution of newborns across each stage to a specific age-class
#' @param birth_female numeric. ratio of males to females in population
#' @param parity logical. parity states imply age distribution of mothers re-scaled to not have parity 0 when Focal born. Default `TRUE`.
#' @param output_kin vector. A vector of particular kin one wishes to obtain results for, e.g., c("m","d","oa"). Default is all kin types.
#' @param summary_kin logical. Results as a data frame of accumulated kin by age of Focal if FALSE, and kin by their age*stage distribution by age of Focal if TRUE.
#' @param sex_Focal character. Female or Male as the user requests
#' @param initial_stage_Focal Numeric in Natural number set {1,2,...,}. The stage which Focal is born into (e.g., 1 for parity 0)
#' @param n_inc numeric. The age/time-increment used in the discretisation of the continuum.
#' @param output_years vector. The times at which we wish to count kin: start year = output_years[1], and end year = output_years[length.]
#'
#' @return A data frame with focal´s age, related ages, stages, sexes, and types of kin for each time-period

#' @export
#'
kin_multi_stage_time_variant_2sex <- function(U_list_females = NULL,
                                              U_list_males = NULL,
                                              F_list_females = NULL,
                                              F_list_males = NULL,
                                              T_list_females = NULL,
                                              T_list_males = NULL,
                                              H_list = NULL,
                                              birth_female = 0.49, ## Sex ratio -- note is 1 - alpha
                                              parity = FALSE,
                                              output_kin = FALSE,
                                              summary_kin = TRUE, # Set to FALSE if we want a full age*stage distribution of kin
                                              sex_Focal = "Female",
                                              initial_stage_Focal = NULL,
                                              n_inc = NULL, ## n_inc is the age-class, time-class increment (e.g., 1year,5year,10year)
                                              output_years){
  
  no_years <- length(U_list_females)
  na <- nrow(U_list_females[[1]])
  ns <- ncol(U_list_females[[1]])
  
  # Ensure inputs are lists of matrices and that the timescale same length
  if(length(U_list_females)!=length(output_years)){stop("Timescale inconsistancy")} ## this is due to my struggles with counting! ( e.g., seq(10, 20, 1) != list(1 : 10) )
  if(!is.list(U_list_females) | !is.list(U_list_males)){stop("U's must be a list with time-series length. Each list entry should be an age*stage dimensional matrix")}
  if(!is.list(F_list_females) | !is.list(F_list_males)){stop("F's must be a list with time-series length. Each list entry should be an age*stage dimensional matrix")}
  if(!is.list(T_list_females) | !is.list(T_list_males)){stop("T's must be a list with time-series length. Each list entry should be an age*stage dimensional matrix")}
  
  ### Define empty lists for the accumulated kin of Focals's life-course -- each list entry will reflect a time-period
  changing_pop_struct <- list()
  Focal_array <- list()
  mom_array <- list()
  gran_array <- list()
  great_gran_array <- list()
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
      kin_out_1 <- all_kin_dy(U_tilde_females,
                              U_tilde_males ,
                              F_tilde_females,
                              F_tilde_males,
                              1-birth_female,
                              na,
                              ns,
                              parity,
                              sex_Focal,
                              initial_stage_Focal)
      ### Relative lists' first entries
      Focal_array[[(1+length(Focal_array))]] <- kin_out_1[["Focal"]]
      daughter_array[[(1+length(daughter_array))]] <- kin_out_1[["d"]]
      grand_daughter_array[[(1+length(grand_daughter_array))]] <- kin_out_1[["gd"]]
      great_grand_daughter_array[[(1+length(great_grand_daughter_array))]] <- kin_out_1[["ggd"]]
      mom_array[[(1+length(mom_array))]] <- kin_out_1[["m"]]
      gran_array[[(1+length(gran_array))]] <- kin_out_1[["gm"]]
      great_gran_array[[(1+length(great_gran_array))]] <- kin_out_1[["ggm"]]
      younger_sis_array[[( 1+length(younger_sis_array))]] <- kin_out_1[["ys"]]
      older_sister_array[[(1+length(older_sister_array))]] <- kin_out_1[["os"]]
      younger_aunt_array[[(1+length(younger_aunt_array))]] <- kin_out_1[["ya"]]
      older_aunt_array[[(1+length(older_aunt_array))]] <- kin_out_1[["oa"]]
      younger_niece_array[[(1+length(younger_niece_array))]] <- kin_out_1[["nys"]]
      older_niece_array[[(1+length(older_niece_array))]] <- kin_out_1[["nos"]]
      younger_cousin_array[[(1+length(younger_cousin_array))]] <- kin_out_1[["cya"]]
      older_cousin_array[[(1+length(older_cousin_array))]] <- kin_out_1[["coa"]]
      changing_pop_struct[[(1+length(changing_pop_struct))]] <- kin_out_1[["ps"]]
      
    }
    updating_Focal <- Focal_array[[year]]
    updating_daughter <- daughter_array[[year]]
    updating_grand_daughter <- grand_daughter_array[[year]]
    updating_great_grand_daughter <- great_grand_daughter_array[[year]]
    updating_mom <- mom_array[[year]]
    updating_gran <- gran_array[[year]]
    updating_great_gran <- great_gran_array[[year]]
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
                             1-birth_female,
                             na,
                             ns,
                             parity,
                             sex_Focal,
                             initial_stage_Focal,
                             updating_Focal,
                             updating_daughter,
                             updating_grand_daughter,
                             updating_great_grand_daughter,
                             updating_mom,
                             updating_gran,
                             updating_great_gran,
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
    Focal_array[[(1+length(Focal_array))]] <- kin_out[["Focal"]]
    daughter_array[[(1+length(daughter_array))]] <- kin_out[["d"]]
    grand_daughter_array[[(1+length(grand_daughter_array))]] <- kin_out[["gd"]]
    great_grand_daughter_array[[(1+length(great_grand_daughter_array))]] <- kin_out[["ggd"]]
    mom_array[[(1+length(mom_array))]] <- kin_out[["m"]]
    gran_array[[(1+length(gran_array))]] <- kin_out[["gm"]]
    great_gran_array[[(1+length(great_gran_array))]] <- kin_out[["ggm"]]
    younger_sis_array[[(1+length(younger_sis_array))]] <- kin_out[["ys"]]
    older_sister_array[[(1+length(older_sister_array))]] <- kin_out[["os"]]
    younger_aunt_array[[(1+length(younger_aunt_array))]] <- kin_out[["ya"]]
    older_aunt_array[[(1+length(older_aunt_array))]] <- kin_out[["oa"]]
    younger_niece_array[[(1+length(younger_niece_array))]] <- kin_out[["nys"]]
    older_niece_array[[(1+length(older_niece_array))]] <- kin_out[["nos"]]
    younger_cousin_array[[(1+length(younger_cousin_array))]] <- kin_out[["cya"]]
    older_cousin_array[[(1+length(older_cousin_array))]] <- kin_out[["coa"]]
    changing_pop_struct[[(1+length(changing_pop_struct))]] <- kin_out[["ps"]]
  }
  tictoc::toc()
  ## create a list of output kin -- each element a time-period specific list of matrices
  ## label the kin names to match DemoKin:
  relative_data <- list("Focal" = Focal_array,
                        "d" = daughter_array,
                        "gd" = grand_daughter_array,
                        "ggd" = great_grand_daughter_array,
                        "m" = mom_array,
                        "gm" = gran_array,
                        "ggm" = great_gran_array,
                        "ys" = younger_sis_array,
                        "os" = older_sister_array,
                        "ya" = younger_aunt_array,
                        "oa" = older_aunt_array,
                        "nys" = younger_niece_array,
                        "nos" = older_niece_array,
                        "cya" = younger_cousin_array,
                        "coa" = older_cousin_array)
  
  relative_names <- names(relative_data)
  ## create a nice data frame output
  if(summary_kin){
    kin_out <- create_cumsum_df(relative_data,
                                relative_names,
                                output_years[1]:output_years[length(output_years)],
                                output_years[1],
                                na,
                                ns,
                                n_inc,
                                output_kin)}
  else{
    kin_out <- create_full_dists_df(relative_data,
                                    relative_names,
                                    output_years[1]:output_years[length(output_years)],
                                    output_years[1],
                                    na,
                                    ns,
                                    n_inc,
                                    output_kin)}
  
  return(kin_out)
}

