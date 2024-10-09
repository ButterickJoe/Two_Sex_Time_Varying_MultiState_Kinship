
### This file prepares data for the UK parity example


`%>%` <- magrittr::`%>%`

################# Historic sex-specific mortality rates
## Each is a list of matrices, of dimension: age*age, and having sub-diagonal survival probs

U_list_F <- readr::read_rds(here::here("Data",  "mort_mats_fem_1938_2070.Rds"))
U_list_M <- readr::read_rds(here::here("Data",  "mort_mats_male_1938_2070.Rds"))

plot(diag(U_list_F[[1]][-1,-ncol(U_list_F[[1]])]))
lines(diag(U_list_M[[1]][-1,-ncol(U_list_M[[1]])]))


################## Parity specific fertility rates 
dat <- readxl::read_excel(here::here("Data", "Parity_births.xlsx"), sheet = "Table", skip=5)
dat <- dat[,c(2,3,16:20)]
dat <- as.data.frame(dat)
dat%>%head()
ncol(dat)
colnames(dat) <- c("Year", "Age", "m1", "m2", "m3", "m4", "m5")

dat%>%head()
dat %>% reshape2::melt(id=c("Year","Age")) %>%
  ggplot2::ggplot(ggplot2::aes(x = Age, y = value, color = Year)) + 
  ggplot2::geom_line(ggplot2::aes(group = factor(Year))) + 
  ggplot2::facet_wrap(~variable) + 
  ggplot2::theme_bw() +
  ggplot2::scale_fill_viridis_c() + 
  ggplot2::scale_colour_viridis_c()

dat <- dat %>% dplyr::transmute(Year = Year,
                                Age = Age,
                                m1 = as.numeric(m1)/1000,
                                m2 = as.numeric(m2)/1000,
                                m3 = as.numeric(m3)/1000,
                                m4 = as.numeric(m4)/1000,
                                m5 = as.numeric(m5)/1000)
dat <- dat %>% dplyr::filter(Year > 1964)
dat%>%head()


start_time_matrix_list <- 1965 - 1938
duration_time <- 2022 - 1965

U_list_M_truncated <- U_list_M[start_time_matrix_list: (start_time_matrix_list + duration_time)] ## start 1970
U_list_F_truncated <- U_list_F[start_time_matrix_list: (start_time_matrix_list + duration_time)]
length(U_list_F_truncated)

## number of ages and stages
na <- 101
ns <- 6
## fertility interval
age_min <- 14
age_max <- 45

U_mat_fem <- list()
U_mat_male <- list()
F_mat_fem <- list()
F_mat_male <- list()
T_mat_fem <- list()
T_mat_male <- list()
H_mat <- list()

for(year in 1965:2022){
  indx_year <- year - 1964
  
  fert_rates_fem <- Matrix::Matrix(nrow = na, ncol = ns, data = 0, sparse = TRUE)
  fert_rates_male <- Matrix::Matrix(nrow = na, ncol = ns, data = 0, sparse = TRUE)
  mort_rates_fem <- Matrix::Matrix(nrow = na, ncol = ns, data = 0, sparse = TRUE)
  mort_rates_male <- Matrix::Matrix(nrow = na, ncol = ns, data = 0, sparse = TRUE)
  
  H_redistribute <- Matrix::Matrix(nrow = na, ncol = ns, data = 0, sparse = TRUE)
  H_redistribute[1,] <- 1
  
  T_trans_list_f <- list()
  T_trans_list_m <- list()
  
  T_no_ptrans <- Matrix::Matrix(nrow = ns, ncol = ns, data = 0, sparse = TRUE)
  diag(T_no_ptrans) <- rep(1, ns) ## diagonal for ages < 14 and > 45 means no leaving P0 or P1
  
  for(stage in 1:ns){
    mort_rates_fem[,stage] <- c(Matrix::diag(U_list_F[[indx_year]][-1,-ncol(U_list_F[[indx_year]])]) , 
                                U_list_F[[indx_year]][na,na])
    mort_rates_male[,stage] <- c(Matrix::diag(U_list_M[[indx_year]][-1,-ncol(U_list_M[[indx_year]])]) , 
                                 U_list_M[[indx_year]][na,na])
  }
  
  for(age in 1:na){
    
    T_female <- Matrix::Matrix(nrow = ns, ncol = ns, data = 0, sparse = TRUE)
    T_male <- Matrix::Matrix(nrow = ns, ncol = ns, data = 0, sparse = TRUE)
    
    if(age >= age_min & age <= age_max){
      dat_temp <- dat %>% dplyr::filter(Year == year, Age == age)
      
      p01 <- (dat_temp$m1/(1+(1-0.5)*dat_temp$m1))
      p12 <- (dat_temp$m2/(1+(1-0.5)*dat_temp$m2))
      p23 <- (dat_temp$m3/(1+(1-0.5)*dat_temp$m3))
      p34 <- (dat_temp$m4/(1+(1-0.5)*dat_temp$m4))
      p45 <- (dat_temp$m5/(1+(1-0.5)*dat_temp$m5))
      
      fert_rates_fem[age,] <- c(p01,p12,p23,p34,p45,p45)
      fert_rates_male[age,] <- c(p01,p12,p23,p34,p45,p45)
      
      
      T_female[1,1] <- 1 - p01
      T_female[2,2] <- 1 - p12
      T_female[3,3] <- 1 - p23
      T_female[4,4] <- 1 - p34
      T_female[5,5] <- 1 - p45
      T_female[6,6] <- 1
      T_female[2,1] <- p01
      T_female[3,2] <- p12
      T_female[4,3] <- p23
      T_female[5,4] <- p34
      T_female[6,5] <- p45
      T_male <- T_female
      
      T_trans_list_f[[(1+length(T_trans_list_f))]] <- T_female
      T_trans_list_m[[(1+length(T_trans_list_m))]] <- T_male
      
    }
    
  }
  
  padd_T_female_juvenile <- list()
  padd_T_male_juvenile <- list()
  padd_T_female_juvenile <- lapply(1:(-1+age_min), function(x){padd_T_female_juvenile[[x]] <- T_no_ptrans})
  padd_T_male_juvenile <- lapply(1:(-1+age_min), function(x){padd_T_male_juvenile[[x]] <- T_no_ptrans})
  padd_T_female_menapausal <- list()
  padd_T_male_menapausal <- list()
  padd_T_female_menapausal <- lapply((1+age_max):na, function(x){padd_T_female_menapausal[[x]] <- T_no_ptrans})
  padd_T_male_menapausal <- lapply((1+age_max):na, function(x){padd_T_male_menapausal[[x]] <- T_no_ptrans})
  
  F_mat_fem[[(1+length(F_mat_fem))]] <- fert_rates_fem
  F_mat_male[[(1+length(F_mat_male))]] <- fert_rates_male
  
  T_mat_fem[[(1+length(T_mat_fem))]] <- c(padd_T_female_juvenile, T_trans_list_f, padd_T_female_menapausal)
  
  U_mat_fem[[(1+length(U_mat_fem))]] <- mort_rates_fem
  U_mat_male[[(1+length(U_mat_male))]] <- mort_rates_male
  
  H_mat[[(1+length(H_mat))]] <- H_redistribute
  
}

dat_out <- here::here("Data")
fs::dir_create(dat_out)

saveRDS(H_mat, paste0(dat_out , "/" , "H_list.Rds"))
saveRDS(F_mat_fem, paste0(dat_out , "/" , "FF_list.Rds"))
saveRDS(F_mat_male, paste0(dat_out , "/" , "FM_list.Rds"))
saveRDS(U_mat_fem, paste0(dat_out , "/" , "UF_list.Rds"))
saveRDS(U_mat_male, paste0(dat_out , "/" , "UM_list.Rds"))
saveRDS(T_mat_fem, paste0(dat_out , "/" , "T_list.Rds"))

