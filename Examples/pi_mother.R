source(here::here("Matrix Model" , "Functions_required.R"))

################## Parity specific stuff 

dat <- read_excel(here::here("Data", "Parity_births.xlsx") ,sheet = "Table", skip=5)
head(dat)
dim(dat)
dat <- dat[,c(2,3,16:20)]
dat <- as.data.frame(dat)
dat%>%head()
ncol(dat)
colnames(dat) <- c("Year","Age","m1","m2","m3","m4","m5")

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
dat%>%filter(Year == 2000)
dat <- dat%>%filter(Year > 1964)

dat$Year%>%unique()

U_list_F <- read_rds(here::here("Data",  "mort_mats_fem_1938_2070.Rds"))
U_list_M <- read_rds(here::here("Data",  "mort_mats_male_1938_2070.Rds"))
length(U_list_F)

start_time_matrix_list <- 1965 - 1938
duration_time <- 2022 - 1965

U_list_M_truncated <- U_list_M[start_time_matrix_list: (start_time_matrix_list + duration_time)] ## start 1970
U_list_F_truncated <- U_list_F[start_time_matrix_list: (start_time_matrix_list + duration_time)]
length(U_list_F_truncated)

year_Ftilde_F_list <- list()
year_Ftilde_M_list <- list()
year_Utilde_M_list <-list()
year_Utilde_F_list <-list()
foreach(year = 1965:1965)%do%{
  indx_year <- year - 1964
  block_diag_fert_mat <- list()
  block_diag_age_mat <- list()
  foreach(age = 1:101)%do%{
    fert_mat <- as(matrix(0,6,6), "sparseMatrix")
    fert_mat[1,] <- 0
    age_mat <- as(matrix(0,6,6), "sparseMatrix")
    diag(age_mat) <- c(1,1,1,1,1,1) ## diagonal for ages < 14 and > 45 means no leaving P0 or P1
    if(age >= 14 & age <= 45){
      dat_temp <- dat%>%filter(Year == year, Age == age)
      
      p01 <- (dat_temp$m1/(1+(1-0.5)*dat_temp$m1))
      p12 <- (dat_temp$m2/(1+(1-0.5)*dat_temp$m2))
      p23 <- (dat_temp$m3/(1+(1-0.5)*dat_temp$m3))
      p34 <- (dat_temp$m4/(1+(1-0.5)*dat_temp$m4))
      p45 <- (dat_temp$m5/(1+(1-0.5)*dat_temp$m5))
      
      age_mat <- as(matrix(0,6,6),"sparseMatrix")
      
      age_mat[1,1] <- 1 - p01
      age_mat[2,2] <- 1 - p12
      age_mat[3,3] <- 1 - p23
      age_mat[4,4] <- 1 - p34
      age_mat[5,5] <- 1 - p45
      age_mat[6,6] <- 1
      age_mat[2,1] <- p01
      age_mat[3,2] <- p12
      age_mat[4,3] <- p23
      age_mat[5,4] <- p34
      age_mat[6,5] <- p45
      
      
      fert_mat <- as(matrix(0,6,6),"sparseMatrix")
      fert_mat[1,] <- c(p01,p12,p23,p34,p45,p45)
    }
    
    block_diag_age_mat[[age]] <- age_mat
    block_diag_fert_mat[[age]] <- fert_mat
  }
  
  U_tran <- block_diag_function(block_diag_age_mat)
  F_tran <- block_diag_function(block_diag_fert_mat)
  
  Um <- U_list_M_truncated[[indx_year]]
  Uf <- U_list_F_truncated[[indx_year]]
  
  male_parity_dep_suv <- block_diag_function(list(Um,
                                                  Um,
                                                  Um,
                                                  Um,
                                                  Um,
                                                  Um))
  
  female_parity_dep_suv <- block_diag_function(list(Uf,
                                                    Uf,
                                                    Uf,
                                                    Uf,
                                                    Uf,
                                                    Uf))
  
  U_tilde_M <- t(K_perm_mat(6,101))%*%male_parity_dep_suv%*%K_perm_mat(6,101)%*%U_tran
  U_tilde_F <- t(K_perm_mat(6,101))%*%female_parity_dep_suv%*%K_perm_mat(6,101)%*%U_tran
  
  year_Utilde_M_list[[(1+length(year_Utilde_M_list))]] <- as(U_tilde_M,"sparseMatrix")
  year_Utilde_F_list[[(1+length(year_Utilde_F_list))]] <- as(U_tilde_F,"sparseMatrix")
  
  H_mat <- matrix(0,101,101)
  H_mat[1,] <- 1
  H_mat <- block_diag_function(list(H_mat,H_mat,H_mat,H_mat,H_mat,H_mat))
  H_mat <- as(H_mat,"sparseMatrix")
  
  F_tilde_M <- t(K_perm_mat(6,101))%*%H_mat%*%K_perm_mat(6,101)%*%F_tran
  #F_tilde_F <- t(K_perm_mat(6,101))%*%H_mat%*%K_perm_mat(6,101)%*%F_tran
  
  year_Ftilde_M_list[[(1+length(year_Ftilde_M_list))]] <- as(F_tilde_M,"sparseMatrix")
  year_Ftilde_F_list[[(1+length(year_Ftilde_F_list))]] <- as(F_tilde_M,"sparseMatrix")
  
}

year_Ftilde_F_list%>%length()
year_Ftilde_M_list%>%length()
year_Utilde_F_list%>%length()
year_Utilde_M_list%>%length()

year_Ftilde_F_list[[1]]%>%dim()
year_Ftilde_M_list[[1]]%>%dim()
year_Utilde_F_list[[1]]%>%dim()
year_Utilde_M_list[[1]]%>%dim()

pi_dists <- pi_mix_parity(year_Utilde_F_list[[1]],year_Utilde_M_list[[1]],year_Ftilde_F_list[[1]],year_Ftilde_M_list[[1]],0.51,101,6)

plot(pi_dists[[2]])

P_gen_fer <- function(U, F, no_ages, no_stages){
  
  A <- U + F
  w <- SD(A)
  lambda <- lambda(A)
  
  P_fer <- as(matrix(0, nrow = no_ages*no_stages, ncol = no_ages*no_stages),"sparseMatrix")
  
  foreach(i = 1 : (no_ages*no_stages))%do%{
    foreach(j = 1 : (no_ages*no_stages))%do%{
      
      P_fer[i,j] <- w[j]*F[i,j]/(w[i]*lambda)
      if(is.na(P_fer[i,j]) == TRUE){P_fer[i,j] = 0}
    }
  }
  return(P_fer)
}

W <- P_gen_fer(year_Utilde_F_list[[1]], (1-0.51)*year_Ftilde_F_list[[1]],101,6)

focal_0 <- rep(0,606)
focal_0[1] <- 1

mom_before <- t(t(focal_0)%*%W)

age_mat_star <- as(matrix(0,6,6),"sparseMatrix")
age_mat_star[1,1] <- 0
age_mat_star[2,2] <- 0
age_mat_star[3,3] <- 0
age_mat_star[4,4] <- 0
age_mat_star[5,5] <- 0
age_mat_star[6,6] <- 1
age_mat_star[2,1] <- 1
age_mat_star[3,2] <- 1
age_mat_star[4,3] <- 1
age_mat_star[5,4] <- 1
age_mat_star[6,5] <- 1
age_mat_star_list <- list()
age_mat_star_list <- lapply(1:101, function(x){age_mat_star_list[[x]] <- age_mat_star})

U_tran_star <- block_diag_function(age_mat_star_list)
Ustar <- t(K_perm_mat(6,101))%*%male_parity_dep_suv%*%K_perm_mat(6,101)%*%U_tran_star
Ustar

mom_gen <- Ustar%*%mom_before
mom_demoKin <- pi_dists[[2]]

mom_gen_df <- data.frame(dist = as.matrix(mom_gen),
                         age = rep(seq(0,100),each=6),
                         stage = rep(seq(1,6),101),
                         method = "Gen: U_star %*% W_1")

sum(mom_gen)

mom_demokin_df <- data.frame(dist = as.matrix(mom_demoKin),
                             age = rep(seq(0,100),each=6),
                             stage = rep(seq(1,6),101),
                             method = "demoKin: pi_f[-parity 0]")

com_df <- rbind(mom_gen_df,mom_demokin_df)
com_df%>%
  ggplot(aes(x = age, y = dist, color = factor(stage), fill = factor(stage))) +
  geom_bar(position = "stack", stat = "identity") + facet_grid(~method)

pi_unequal <- com_df%>%
  ggplot(aes(x = age, y = dist, color = factor(stage))) +
  geom_point(aes(shape = method), size = 1) + facet_wrap(~stage) +  theme_bw() +
  theme(text = element_text(size = 12)) + theme(legend.position = "top")


pi_unequal
ggsave(paste0(fig_out, "/unequal_pi.png"), pi_unequal, width  =8, height = 6)

pi_distsNP <- pi_mix(year_Utilde_F_list[[1]],year_Utilde_M_list[[1]],year_Ftilde_F_list[[1]],year_Ftilde_M_list[[1]],0.51,101,6)
mom_np <- pi_distsNP[[2]]

mom_before_df <- data.frame(dist = as.matrix(mom_before),
                            age = rep(seq(0,100),each=6),
                            stage = rep(seq(1,6),101),
                            method = "Gen: W_1, at t=-1")
mom_pi_np_df <- data.frame(dist = as.matrix(mom_np),
                           age = rep(seq(0,100),each=6),
                           stage = rep(seq(1,6),101),
                           method = "demoKin: pi_f, at t=0")
cm_df <- rbind(mom_before_df, mom_pi_np_df) 

pi_equal <- cm_df%>%
  ggplot(aes(x = age, y = dist, color = factor(stage))) +
  geom_point(aes(shape = method), size = 1) + facet_wrap(~stage)  + theme_bw() +
  theme(text = element_text(size = 12)) + theme(legend.position = "top")
fig_out <- here::here("Figures","Musing on pi")
ggsave(paste0(fig_out, "/equal_pi.png"), pi_equal, width = 8, height = 6)
