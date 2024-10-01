#####################

#################### Create a data frame from the time-varying output

## Args : dat_list = a list of the matrices X_foc, X_child, ... , each with year/age/sex specific kin distributions

##        list_dist = a character list of which kin we want to analyse, e.g., list = "parents", "cousins", "grandpars"

##        years = sequence of years in the time series, e.g., 1950:2050 or seq(1950,2050,1)

##        Start year of inputs

## na = number of ages and ns = number of stages

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

