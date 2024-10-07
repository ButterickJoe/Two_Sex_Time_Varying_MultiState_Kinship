

source(here::here("Matrix Model" , "Two_Sex_Time_Variant_MultiState_Kinship_tandem_NL.R"))

`%>%` <- magrittr::`%>%`

################# Historic sex-specific mortality rates
U_list_F <- readr::read_rds(here::here("Data",  "mort_mats_fem_1938_2070.Rds"))
U_list_M <- readr::read_rds(here::here("Data",  "mort_mats_male_1938_2070.Rds"))

################## Parity specific fertility rates 
dat <- readxl::read_excel(here::here("Data", "Parity_births.xlsx"), sheet = "Table", skip=5)
dat <- dat[,c(2,3,16:20)]
dat <- as.data.frame(dat)
dat %>% head()
ncol(dat)
colnames(dat) <- c("Year", "Age", "m1", "m2", "m3", "m4", "m5")

dat %>% head()
dat %>% reshape2::melt(id=c("Year","Age")) %>%
  ggplot2::ggplot(ggplot2::aes(x = Age, y = value, color = Year)) + 
  ggplot2::geom_line(ggplot2::aes(group = factor(Year))) + 
  ggplot2::facet_wrap(~variable) + ggplot2::theme_bw() +
  ggplot2::scale_fill_viridis_c() + ggplot2::scale_colour_viridis_c()

dat <- dat %>% dplyr::transmute(Year = Year,
                       Age = Age,
                       m1 = as.numeric(m1)/1000,
                       m2 = as.numeric(m2)/1000,
                       m3 = as.numeric(m3)/1000,
                       m4 = as.numeric(m4)/1000,
                       m5 = as.numeric(m5)/1000)
dat <- dat %>% dplyr::filter(Year > 1964)
dat %>% head()


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
rm(dat)
rm(U_list_F)
rm(U_list_M)
rm(U_list_F_truncated)
rm(U_list_M_truncated)
gc()

kin_out_1965_TEST_NL <- 
  kin_multi_stage_TV_2_sex_tandem_NL(U_mat_fem[1:40], ## 40 years time-series of demographic rates (1965-2005) see last input...
                         U_mat_male[1:40],
                         F_mat_fem[1:40],
                         F_mat_male[1:40],
                         T_mat_fem[1:40],
                         T_mat_fem[1:40],
                         H_mat[1:40],
                         alpha = 0.51, ## Sex ratio -- UK here
                         parity = TRUE, ## this example uses parity
                         dist_output = FALSE, ## here we only want to know accumulated kin by age of Focal
                         sex_Focal = "Female", ##  define Focal's sex at birth 
                         stage_Focal = 1, ## Define Focals stage at birth -- if parity automatically set to 1
                         seq(1965,(1965+40)))



fig_out <- here::here("Figures","Cumulated_kin")
fs::dir_create(fig_out)

parents_fig <- kin_out_1965_TEST_NL %>% 
  dplyr::filter(group == "parents", year %in% c(1965, 1966, 1967, 1968)) %>%
  ggplot2::ggplot(ggplot2::aes(x = Age_Focal, y = pred_no_kin, color = Stage, fill = Stage)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") + 
  ggplot2::facet_grid(Sex~year) +
  ggplot2::scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) + 
  ggplot2::ylab("Parents") +
  ggplot2::theme(legend.position = "top")
parents_fig
ggsave(paste0(fig_out, "/parents_example.png"), parents_fig, width = 8, height = 6)


g_parents_fig <- kin_out_1965_TEST_NL %>% 
  dplyr::filter(group == "grand parents", year %in% c(1965, 1966, 1967, 1968)) %>%
  ggplot2::ggplot(ggplot2::aes(x = Age_Focal, y = pred_no_kin, color = Stage, fill = Stage)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") + 
  ggplot2::facet_grid(Sex~year) +
  ggplot2::scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) + 
  ggplot2::ylab("Grand parents") +
  ggplot2::theme(legend.position = "top")
g_parents_fig

ggsave(paste0(fig_out, "/grand_parents_example.png"), g_parents_fig, width = 8, height = 6)


offspring_fig <- kin_out_1965_TEST_NL %>% 
  dplyr::filter(group == "offspring", year %in% c(1965, 1966, 1967, 1968))%>%
  ggplot2::ggplot(ggplot2::aes(x = Age_Focal, y = pred_no_kin, color = Stage, fill = Stage)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") + 
  ggplot2::facet_grid(Sex~year) +
  ggplot2::scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) + 
  ggplot2::ylab("Offspring") +
  ggplot2::theme(legend.position = "top")
offspring_fig
ggsave(paste0(fig_out, "/offspring_example.png"), offspring_fig, width = 8, height = 6)

g_offspring_fig <- kin_out_1965_TEST_NL %>% 
  dplyr::filter(group == "grand offspring", year %in% c(1965, 1966, 1967, 1968))%>%
  ggplot2::ggplot(ggplot2::aes(x = Age_Focal, y = pred_no_kin, color = Stage, fill = Stage)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") + 
  ggplot2::facet_grid(Sex~year) +
  ggplot2::scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) + 
  ggplot2::ylab("Grand offspring")+
  ggplot2::theme(legend.position = "top")
g_offspring_fig
ggsave(paste0(fig_out, "/grand_offspring_example.png"), g_offspring_fig, width = 8, height = 6)



o_au_fig <- kin_out_1965_TEST_NL %>% 
  dplyr::filter(group == "older aunt/unlce", year %in% c(1965, 1966, 1967, 1968))%>%
  ggplot2::ggplot(ggplot2::aes(x = Age_Focal, y = pred_no_kin, color = Stage, fill = Stage)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") + 
  ggplot2::facet_grid(Sex~year) +
  ggplot2::scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) + 
  ggplot2::ylab("Older aunts and uncles") +
  ggplot2::theme(legend.position = "top")
o_au_fig
ggsave(paste0(fig_out, "/older_aunt_uncle_example.png"), o_au_fig, width = 8, height = 6)

y_au_fig <- kin_out_1965_TEST_NL %>% 
  dplyr::filter(group == "younger aunt/unlces", year %in% c(1965, 1966, 1967, 1968))%>%
  ggplot2::ggplot(ggplot2::aes(x = Age_Focal, y = pred_no_kin, color = Stage, fill = Stage)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") + 
  ggplot2::facet_grid(Sex~year) +
  ggplot2::scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) + 
  ggplot2::ylab("Younger aunts and uncles") +
  ggplot2::theme(legend.position = "top")
y_au_fig
ggsave(paste0(fig_out, "/younger_aunt_uncle_example.png"), y_au_fig, width = 8, height = 6)

o_nn_fig <- kin_out_1965_TEST_NL %>% 
  dplyr::filter(group == "older niece/nephews", year %in% c(1965, 1966, 1967, 1968))%>%
  ggplot2::ggplot(ggplot2::aes(x = Age_Focal, y = pred_no_kin, color = Stage, fill = Stage)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") + 
  ggplot2::facet_grid(Sex~year) +
  ggplot2::scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) + 
  ggplot2::ylab("Older nieces and nephews") +
  ggplot2::theme(legend.position = "top")
o_nn_fig
ggsave(paste0(fig_out, "/older_niece_nephew_example.png"), o_nn_fig, width = 8, height = 6)

y_nn_fig <- kin_out_1965_TEST_NL %>% 
  dplyr::filter(group == "younger niece/nephews", year %in% c(1965, 1966, 1967, 1968))%>%
  ggplot2::ggplot(ggplot2::aes(x = Age_Focal, y = pred_no_kin, color = Stage, fill = Stage)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") + 
  ggplot2::facet_grid(Sex~year) +
  ggplot2::scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) + 
  ggplot2::ylab("Younger nieces and nephews")
y_nn_fig
ggsave(paste0(fig_out, "/younger_niece_nephew_example.png"), y_nn_fig)


oc_fig <- kin_out_1965_TEST_NL %>% 
  dplyr::filter(group == "older cousin", year %in% c(1965, 1966, 1967, 1968))%>%
  ggplot2::ggplot(ggplot2::aes(x = Age_Focal, y = pred_no_kin, color = Stage, fill = Stage)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") + 
  ggplot2::facet_grid(Sex~year) +
  ggplot2::scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) + 
  ggplot2::ylab("Older cousins") +
  ggplot2::theme(legend.position = "top")

oc_fig
ggsave(paste0(fig_out, "/older_cousin_example.png"), oc_fig, width = 8, height = 6)

yc_fig <- kin_out_1965_TEST_NL %>% 
  dplyr::filter(group == "younger cousin", year %in% c(1965, 1966, 1967, 1968))%>%
  ggplot2::ggplot(ggplot2::aes(x = Age_Focal, y = pred_no_kin, color = Stage, fill = Stage)) +
  ggplot2::geom_bar(position = "stack", stat = "identity") + 
  ggplot2::facet_grid(Sex~year) +
  ggplot2::scale_x_continuous(breaks = c(0,10,20,30,40,50,60,70,80,90,100)) + 
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5)) + 
  ggplot2::ylab("Younger cousins") +
  ggplot2::theme(legend.position = "top")
yc_fig
ggsave(paste0(fig_out, "/younger_cousin_example.png"), yc_fig, width = 8, height = 6)
