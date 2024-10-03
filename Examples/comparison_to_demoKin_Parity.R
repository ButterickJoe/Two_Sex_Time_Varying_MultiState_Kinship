

source(here::here("Matrix Model" , "Two_Sex_Time_Variant_MultiState_Kinship_tandem.R"))

################# Historic sex-specific mortality rates
U_list_F <- read_rds(here::here("Data",  "mort_mats_fem_1938_2070.Rds"))
U_list_M <- read_rds(here::here("Data",  "mort_mats_male_1938_2070.Rds"))

################## Parity specific fertility rates 
dat <- read_excel(here::here("Data", "Parity_births.xlsx"), sheet = "Table", skip=5)
dat <- dat[,c(2,3,16:20)]
dat <- as.data.frame(dat)
dat%>%head()
ncol(dat)
colnames(dat) <- c("Year", "Age", "m1", "m2", "m3", "m4", "m5")

dat%>%head()
dat%>%melt(id=c("Year","Age"))%>%ggplot(aes(x = Age, y = value, color = Year)) + 
  geom_line(aes(group = factor(Year))) + facet_wrap(~variable) + theme_bw() +
  scale_fill_viridis_c() + scale_colour_viridis_c()

dat <- dat%>%transmute(Year = Year,
                       Age = Age,
                       m1 = as.numeric(m1)/1000,
                       m2 = as.numeric(m2)/1000,
                       m3 = as.numeric(m3)/1000,
                       m4 = as.numeric(m4)/1000,
                       m5 = as.numeric(m5)/1000)
dat <- dat%>%filter(Year > 1964)
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

foreach(year = 1965:1966)%do%{
  indx_year <- year - 1964

  fert_rates_fem <- as(matrix(0, nrow = na, ncol = ns), "sparseMatrix")
  fert_rates_male <- as(matrix(0, nrow = na, ncol = ns), "sparseMatrix")
  mort_rates_fem <- as(matrix(0, nrow = na, ncol = ns), "sparseMatrix")
  mort_rates_male <- as(matrix(0, nrow = na, ncol = ns), "sparseMatrix")
  
  H_redistribute <- as(matrix(0, nrow = na, ncol = ns), "sparseMatrix")
  H_redistribute[1,] <- 1
  
  T_trans_list_f <- list()
  T_trans_list_m <- list()
  
  T_no_ptrans <- as(matrix(0, nrow = ns, ncol = ns), "sparseMatrix")
  diag(T_no_ptrans) <- rep(1, ns) ## diagonal for ages < 14 and > 45 means no leaving P0 or P1
  
  foreach(stage = 1:ns)%do%{
    mort_rates_fem[,stage] <- c(diag(U_list_F[[indx_year]][-1,-ncol(U_list_F[[indx_year]])]), U_list_F[[indx_year]][na,na])
    mort_rates_male[,stage] <- c(diag(U_list_M[[indx_year]][-1,-ncol(U_list_M[[indx_year]])]), U_list_M[[indx_year]][na,na])
  }
  
  foreach(age = 1:na)%do%{
    
    T_female <- as(matrix(0, nrow = ns, ncol = ns), "sparseMatrix")
    T_male <- as(matrix(0, nrow = ns, ncol = ns), "sparseMatrix")
    
    if(age >= age_min & age <= age_max){
      dat_temp <- dat%>%filter(Year == year, Age == age)
      
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
rm(dat)
rm(U_list_F)
rm(U_list_M)
rm(U_list_F_truncated)
rm(U_list_M_truncated)
gc()

kin_out_1965_new <- 
  kin_multi_stage_TV_2_sex_tandem(U_mat_fem[1:1], ## 40 years time-series of demographic rates (1965-2005) see last input...
                         U_mat_male[1:1],
                         F_mat_fem[1:1],
                         F_mat_male[1:1],
                         T_mat_fem[1:1],
                         T_mat_fem[1:1],
                         H_mat[1:1],
                         alpha = 0.51, ## Sex ratio -- UK here
                         parity = TRUE, ## this example uses parity
                         dist_output = FALSE, ## here we only want to know accumulated kin by age of Focal
                         sex_Focal = "Female", ##  define Focal's sex at birth 
                         stage_Focal = 1, ## Define Focals stage at birth -- if parity automatically set to 1
                         seq(1965,(1965+1)))



kin_out_1965_new$group%>%unique()
kin_out_1965_new%>%filter(group == "parents", Sex == "Female", year == 1965)%>%
  dplyr::select(Age_Focal,Stage,pred_no_kin)%>%
  ggplot(aes(x = Age_Focal, y = pred_no_kin, color = Stage, fill = Stage)) +
  geom_bar(position = "stack", stat = "identity")




df_export <- kin_out_1965_new%>%filter(Sex == "Female", year == 1965)%>%
  dplyr::select(Age_Focal,Stage,pred_no_kin,group)%>%transmute(stage_kin = Stage,
                                                               age_focal = Age_Focal,
                                                               living = pred_no_kin, 
                                                               kin = group)
df_export$kin%>%unique()
df_export$kin <- ifelse(df_export$kin == "Focal", "focal",
                        ifelse(df_export$kin == "parents", "m",
                               ifelse(df_export$kin == "offspring", "d",
                                      ifelse(df_export$kin == "younger siblings", "ys",
                                             ifelse(df_export$kin == "older siblings", "os",
                                                    ifelse(df_export$kin == "younger aunt/uncle", "ya",
                                                           ifelse(df_export$kin == "older aunt/unlce", "oa", 
                                                                  ifelse(df_export$kin == "older cousin", "coa", NA))))))))


saveRDS(df_export,"C:/Users/jb4u23/OneDrive - University of Southampton/Butterick year 1/R code/First files/Implementing Caswell/DemoKin/DemoKin-main/Butterick_parity_comparison/Joes_mother_parity_new.Rds")



