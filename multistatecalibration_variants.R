# Variant testing 

# Variant data  ###########################################################

variant_data <- read.csv("data/variants/SARS-CoV-2_Variant_Proportions_20241105.csv")

variant_data_clean <- variant_data %>%
  mutate(week = as.Date(str_sub(week_ending, 1, 10), format = "%Y-%m-%d"),
         published_week = as.Date(str_sub(published_date, 1, 10), format = "%Y-%m-%d"),
         variant_group = case_when(variant %in% c("BA.1.1", "B.1.1.529") ~ "BA.1.1",
                                   variant %in% c("BQ.1", "BQ.1.1") ~ "BQ",
                                   variant %in% c("BA.2", "BA.2.12.1") ~ "BA.2",
                                   variant %in% c("BA.4", "BA.5") ~ "BA.4/5",
                                   variant %in% c("XBB.1.16", "XBB.1.9.1", "XBB.2.3", "XBB.1.5") ~ "XBB",
                                   variant %in% c("EG.5", "HV.1") ~ "EG/HV",
                                   variant %in% c("JN.1", "JN.1.7") ~ "JN.1",
                                   variant %in% c("KP.2", "KP.2.3", "KP.3", "KP.3.1.1", "LB.1") ~ "KP/LB",
                                   TRUE ~ "Low-Circulating Variants"),
         variant_group = if_else(share < 0.1  & !(variant_group %in% c("B.1.1.529", "BA.1.1", "BA.2", "BA.4/5", "BQ", "EG.5", "HV.1", "JN.1", "XBB")),  "Low-Circulating Variants\n< 0.1 variant share", variant_group),
         weeks_since = as.numeric(interval(as.Date('2022-01-01'), week) %/% weeks(1)),
         modeltype = if_else(weeks_since %in% c(144:147), "weighted", modeltype)) %>%
  filter(usa_or_hhsregion %in% c(1, 2, 3, 4, 5, 6, 8, 9, 10),
         share > 0.1,
         variant_group != "Low-Circulating Variants",
         week >= as.Date("2021-12-01") & week <= as.Date("2024-10-31"),
         time_interval == "biweekly",
         modeltype == "weighted",
         share > 0) %>%
  group_by(usa_or_hhsregion, variant, week, weeks_since) %>%
  mutate(num_vars = if_else(variant_group == "Low-Circulating Variants\n< 0.1 variant share", n(), 1)) %>% ungroup() %>%
  group_by(usa_or_hhsregion, variant, week, weeks_since) %>% filter(published_week == max(published_week)) %>% dplyr::summarise(share = sum(share)) %>%
  ungroup() 
head(variant_data_clean)

# Disaggregate biweekly -> weekly
variant_data_filled <- variant_data_clean %>%
  group_by(usa_or_hhsregion, variant) %>%
  mutate(lag_share = lag(share),
         weeks_since_lag = weeks_since - 1,
         weeks_lag = week %m-% weeks(1),
         average_share = (lag_share + share) / 2) %>%
  drop_na() %>%
  dplyr::select(usa_or_hhsregion, variant, weeks_lag, weeks_since_lag, average_share) %>%
  dplyr::rename(week = weeks_lag, weeks_since = weeks_since_lag, share = average_share)

combined_weekly <- as.data.frame(rbind(variant_data_clean, variant_data_filled)) %>%
  group_by(usa_or_hhsregion, week, weeks_since) %>%
  mutate(total_shares_by_week = sum(share),
         rescaled_share = share/total_shares_by_week) 
head(combined_weekly)

ggplot(data = combined_weekly, aes(x = week, y = share, col = variant)) +
  geom_line(lwd = 1.4) +
 # geom_line(aes(x = week, y = rescaled_share)) +
  facet_wrap(vars(usa_or_hhsregion)) +
  theme_bw()

# Translate frequencies into single mulipliter reflecting variance dominance. 
# at every single week point 
frequency= combined_weekly %>%
  group_by(usa_or_hhsregion, week) %>% 
  mutate(max_frequency = max(share)) %>%
  distinct(week, max_frequency) %>%
  mutate(frequency_complement = 1 - max_frequency) %>%
  mutate(odds_frequency = max_frequency/frequency_complement) %>%
  ungroup() %>%
  mutate(rescaled_lof = odds_frequency) %>%
  mutate(week= as.Date(week)) %>%
  group_by(usa_or_hhsregion) %>%
#  mutate(immunity_scalar = ifelse(max_frequency > .60, 1.02, 1)) %>%
  mutate(
    above60 = max_frequency > 0.70,
    consecutive_weeks_above60 = {
      rl <- rle(above60)
      # for TRUE runs produce 1..n, for FALSE runs produce 0 repeated
      unlist(mapply(function(val, len) if (val) seq_len(len) else rep(0, len),
                    rl$values, rl$lengths, SIMPLIFY = FALSE))
    }) %>%
  ungroup() %>%
  group_by(usa_or_hhsregion) %>%
  mutate(slope = max_frequency - lag(max_frequency)) %>%
  ungroup() %>%
  mutate(infusion = ifelse(slope > 0 & consecutive_weeks_above60 > 1, 1.01, 1))

head(frequency)

ggplot(data = frequency, aes(x = week, y = infusion, col = usa_or_hhsregion )) +
  geom_line(lwd = 1.5) +
  theme_bw() +
  facet_wrap(vars(usa_or_hhsregion)) +
  ylab("Consecutive weeks above 60%") + xlab("Date") +
  scale_color_viridis_d() +
  ggtitle('CDC variant date')

# Visualize 
ggplot(data = frequency, aes(x = week, y = max_frequency, col = usa_or_hhsregion )) +
  geom_line(lwd = 1.5) +
  theme_bw() +
  facet_wrap(vars(usa_or_hhsregion)) +
  ylab("Max share") + xlab("Date") +
  scale_color_viridis_d() +
  ggtitle('CDC variant date')


Location = c("CA", "CO", "CT", "GA", "MD",
             "MI", "MN", "NM", "NY", "OH", 
             "OR", "TN", "UT")
usa_or_hhsregion = c(9, 8, 1, 4, 3, 5, 5, 
                     6, 2, 5, 10, 4, 8)
location_hhs = data.frame(Location, usa_or_hhsregion) %>%
  mutate(usa_or_hhsregion = as.character(usa_or_hhsregion))


variant_dataset = left_join(location_hhs, frequency, by = "usa_or_hhsregion", 
                            relationship = "many-to-many") %>%
  dplyr::select(Location, infusion, week, max_frequency) %>% 
  mutate(date = week) %>% dplyr::select(-week) %>%
  mutate(week = week(date), year = year(date) )

# Visualize 
ggplot(data = variant_dataset, aes(x = week/52 + year, y = infusion, col = Location )) +
  geom_line(lwd = 1.5) +
  theme_bw() +
  ylab("Max share") + xlab("Date") +
  scale_color_viridis_d() +
  ggtitle('Max frequency of variant observed in a given week') +
  facet_wrap(vars(Location))


 # let's try this derivative thing 
#head(variant_dataset)
#variant_dataset_slope = variant_dataset %>%
#  group_by(Location) %>%
#  mutate(freq_lag1 = lag(max_frequency)) %>%
#  mutate(rise =  max_frequency - freq_lag1 ) %>%
# # mutate(rise_smooth = rollmean(rise, k = 2, na.pad = TRUE)) %>%
#  mutate(boost = ifelse(rise > 0.0 & lag(rise) >  0.0 & lag(rise, 2) > 0
#                        & lag(rise, 3) > 0, 1.2, 1 )) %>%
#  ungroup() 


#variant_dataset_slope <- variant_dataset %>%
#  group_by(Location) %>%
#  mutate(
#    freq_lag1 = lag(max_frequency),
#    rise = max_frequency - freq_lag1,
#    rise_smooth = rollmean(rise, k = 2, na.pad = TRUE)
#  ) %>%
#  mutate(
#    pos_run = with(rle(rise > 0), rep(seq_along(lengths), lengths)),
#    run_length = with(rle(rise > 0), rep(lengths, lengths)),
#    boost = ifelse(rise > 0 & run_length >= 4, 1.2, 1)
#  ) %>%
#  ungroup()

#ggplot(data = variant_dataset_slope, aes(x = week/52 + year, y = rise, col = Location )) +
#  geom_line(lwd = 1.5) +
#  theme_bw() +
#  ylab("Slope") + xlab("Date") +
#  scale_color_viridis_d() +
#  ggtitle('Slope') +
#  facet_wrap(vars(Location)) +
#  geom_hline(yintercept = 0, lty = "dashed")

#ggplot(data = variant_dataset_slope, aes(x = week/52 + year, y = boost, col = Location )) +
#  geom_line(lwd = 1.5) +
#  theme_bw() +
#  ylab("Slope") + xlab("Date") +
#  scale_color_viridis_d() +
#  ggtitle('Max frequency of variant observed in a given week') +
#  facet_wrap(vars(Location))# +
 # geom_hline(yintercept = 0, lty = "dashed")


####################################################
# lets only look at specific variants 
variant_data_filled <- variant_data_clean %>%
  group_by(usa_or_hhsregion, variant) %>%
  mutate(lag_share = lag(share),
         weeks_since_lag = weeks_since - 1,
         weeks_lag = week %m-% weeks(1),
         average_share = (lag_share + share) / 2) %>%
  drop_na() %>%
  dplyr::select(usa_or_hhsregion, variant, weeks_lag, weeks_since_lag, average_share) %>%
  dplyr::rename(week = weeks_lag, weeks_since = weeks_since_lag, share = average_share)

combined_weekly <- as.data.frame(rbind(variant_data_clean, variant_data_filled)) %>%
  group_by(usa_or_hhsregion, week, weeks_since) %>%
  mutate(total_shares_by_week = sum(share),
         rescaled_share = share/total_shares_by_week) %>%
  filter(variant %in% c("BA.5", "XBB.1.5", "JN.1")) %>%
  group_by(variant, usa_or_hhsregion ) %>%
  mutate(slope = share - lag(share)) %>%
  mutate(slope_pos = ifelse(slope > 0, 1, 0)) %>%
  mutate(across(slope_pos, ~ replace(., is.na(.), 0))) %>%
  mutate(slope_pos = replace(slope_pos, week > "2022-09-03" & variant == "BA.5", 0)) %>%
  mutate(slope_pos = replace(slope_pos, week > "2023-06-03" & variant == "XBB.1.5", 0)) %>%
  mutate(slope_pos = replace(slope_pos, week > "2024-02-01" & variant == "JN.1", 0)) 



# Plot variant 
ggplot(data = combined_weekly, aes(x = week, y = share, col = variant)) +
  geom_line(lwd = 1.4) +
  # geom_line(aes(x = week, y = rescaled_share)) +
  facet_wrap(vars(usa_or_hhsregion)) +
  theme_bw() 

ggplot(data = combined_weekly , aes(x = week, y = slope_pos , col = variant)) +
  geom_line(lwd = 1.4) +
  # geom_line(aes(x = week, y = rescaled_share)) +
  facet_wrap(vars(usa_or_hhsregion)) +
  theme_bw() 


Location = c("CA", "CO", "CT", "GA", "MD",
             "MI", "MN", "NM", "NY", "OH", 
             "OR", "TN", "UT")
usa_or_hhsregion = c(9, 8, 1, 4, 3, 5, 5, 
                     6, 2, 5, 10, 4, 8)
location_hhs = data.frame(Location, usa_or_hhsregion) %>%
  mutate(usa_or_hhsregion = as.character(usa_or_hhsregion))


variant_dataset = left_join(location_hhs, combined_weekly, by = "usa_or_hhsregion", 
                            relationship = "many-to-many") %>%
  dplyr::select(Location, week, slope_pos) %>% 
  mutate(date = week) %>% dplyr::select(-week) %>%
  mutate(week = week(date), year = year(date) )



#####################################################################
# try to run it in the model 
setwd("/Users/sambents/Desktop/Lo/covid_seasonality/data/processed")

# DATA ###################################################################################
Location = c( "CT", "GA", "MI",  "MN", "NY",  "OH", "TN", "CO")
pop_size_2023 = c( 3643023,  11064432, 10083356, 5753048, 19737367, 11824034,  7148304 ,  5901339)
#Location = c( "CA", "CO", "NM", "OR", "UT")
#pop_size_2023 = c( 39198693,   5901339,  2121164, 4253653, 3443222)
census_dat = data.frame(Location, pop_size_2023) 
#Location = c( "CA", "OR", "UT")
#pop_size_2023 = c( 39198693, 4253653, 3443222)
#census_dat = data.frame(Location, pop_size_2023) 

# Load in data 
full_calibration_dat = read.csv("processed_calibration_dat.csv")  %>%
  filter(Location %in% c( "CT", "MN", "OH", "MI", "NY", "TN", "GA", "CO" )) %>%
 # filter(Location %in% c("CA", "CO", "NM", "OR", "UT" )) %>%
 # filter(Location %in% c("CA",  "OR", "UT" )) %>%
  left_join(census_dat, by = "Location")  %>%
  left_join(variant_dataset, by = c("Location", "week", "year")) %>%
  mutate(across(slope_pos, ~ replace(., is.na(.), 0)))  %>%
  dplyr::select(-date)
  
head(full_calibration_dat )

ggplot(data = full_calibration_dat %>% group_by(Location) %>% filter(Location != "CO"), aes(x = t, y = slope_pos)) +
  facet_wrap(vars(Location)) + 
  # geom_line(aes(x = t, y = hosp/max(hosp)), col = "red") +
  geom_line() +
  theme_bw() #+
 # geom_hline(yintercept  = 40)


# MCMC ###################################################################################

run_hierarchical_mcmc_calibration <- function(full_calibration_dat, iterations = 18000, burn_in = 8000, thin = 2, n_chains = 3) {
  
  # check states to be fit in the calibration
  states <- unique(full_calibration_dat$Location)
  n_states <- length(states)
  
  # state populations as of the 2023 census 
  state_populations <- list(
       "CT" = 3643023, 
       "GA" = 11064432, 
       "MN" =  5753048, 
       "OH" = 11824034,
       "MI" = 10083356, 
       "NY" = 19737367, 
       "TN" = 7148304,
  #  "CA" =  39198693,
     "CO" = 5901339 
  #  "NM" = 2121164, 
 #   "OR" = 4253653, 
  #   "UT" = 3443222
    )
  
  # SEIRS model code 
  seirs_model_hybrid <- function(t, state, parms) {
    
    # name state variables 
    S <- state[1]
    E <- state[2] 
    I <- state[3]
    R <- state[4]
    V <- state[5]
    S2 <- state[6]
    E2 <- state[7]
    I2 <- state[8]
    R2 <- state[9]
    V2 <- state[10]
    
    # model parameters 
   # var <- parms[["var"]]
    var_scalar = parms[["var_scalar"]]
    hum <- parms[["hum"]] 
    temp <- parms[["temp"]]
    lp <- parms[["lp"]]
    gamma <- parms[["gamma"]]
    imm <- parms[["imm"]]
    imm_vax <- parms[["imm_vax"]]
    VE <- parms[["VE"]]
    hr1 <- parms[["hr1"]]
    hr2 <- parms[["hr2"]]
    beta <- parms[["beta"]]
    N <- parms[["N"]]
    
    # state calibration data 
    state_data <- parms[["state_data"]]
    
    # time, specify it must be less than the number of rows in calibration data 
    t_idx <- round(t)
    if(t_idx < 1) t_idx <- 1
    if(t_idx > nrow(state_data)) t_idx <- nrow(state_data)
    
    # pull in time series influencing beta 
    vax_rate_at_t <- state_data$incident[t_idx]
    humidity_at_t <- state_data$humid_smooth[t_idx]
    temp_at_t <- state_data$inverse_temp[t_idx]
 #   variant_at_t <- state_data$odds_frequency[t_idx]
    
    # seasonal function
    seasonal <- beta * ( hum * humidity_at_t + temp * temp_at_t  )
    
    # make seasonal value something small and positive if negative- this should not occur to start with 
    if(!is.finite(seasonal) || seasonal < 0) seasonal <- 0.01
    
    # make time varying immunity parameter 
    
    state_data_imm = state_data %>%
      mutate(immunity_scalar_fitted = ifelse(infusion > 1.005, infusion*var_scalar, infusion))
    
    
    imm_vec = imm* state_data_imm$immunity_scalar_fitted[t_idx] 
    
    #  diff eq's 
    dS <- -seasonal * S * (I + I2) / N -  VE * vax_rate_at_t * S + imm_vax * V
    dE <- seasonal * S * (I + I2) / N - lp * E
    dI <- lp * E - gamma * I
    dR <- I * gamma - R * imm_vec
    dV <- VE * vax_rate_at_t * S - imm_vax * V
    
    dS2 <- -seasonal * S2 * (I + I2) / N + imm_vec * R + imm_vec * R2 + imm_vax * V2 - VE * vax_rate_at_t * S2
    dE2 <- seasonal * S2 * (I + I2) / N - lp * E2
    dI2 <- lp * E2 - gamma * I2
    dR2 <- I2 * gamma - R2 * imm_vec
    dV2 <- VE * vax_rate_at_t * S2 - imm_vax * V2
    
    return(list(c(dS, dE, dI, dR, dV, dS2, dE2, dI2, dR2, dV2)))
  }
  
  # likelihood function 
  likelihood_func <- function(params) {
    
    # initialize likelihood
    total_log_lik <- 0
    
    # loop to fit each state 
    for(i in 1:n_states) {
      
      state_name <- states[i]
      state_data <- full_calibration_dat[full_calibration_dat$Location == state_name, ]
      state_data <- state_data[order(state_data$t), ]
      
      # state-specific parameters
      beta_param <- paste0("beta_", state_name)
      
      beta_value <- params[[beta_param]]
      
      S0_value <- params[["S0"]] 
      
      N_val <- state_populations[[state_name]]
      
      # state-specific initial conditions 
      state_init <- c(
        S = S0_value * N_val,
        E = 0.01 * N_val,
        I = 0.001 * N_val,
        R = (1 - S0_value - 0.01 - 0.001 - 0.049) * N_val,
        V = 0.049 * N_val,
        S2 = 0, E2 = 0, I2 = 0, R2 = 0, V2 = 0
      )
      
      # parmraters used in state-specific model 
      all_params <- list(
      #  var = params["var"], 
        var_scalar = params["var_scalar"],
        hum = params["hum"], temp = params["temp"],
        imm = params["imm"], imm_vax = params["imm_vax"],
        hr1 = params["hr1"], hr2 = params["hr2"],
        beta = beta_value, N = N_val,
        lp = 1/(4/7), gamma = 1/(4/7), VE = 0.30,
        state_data = state_data
      )
      
      # time vector 
      time_vec <- 1:nrow(state_data)
      
      # solve ODE, if there is an error, return very large error 
      out <- tryCatch({
        ode(y = state_init, times = time_vec, func = seirs_model_hybrid, parms = all_params)
      }, error = function(e) {
        return(NULL)
      })
      
      if(is.null(out)) {
        return(-1e6)
      }
      
      # if model is solved, save output and calculate hospitalizations 
      out_df <- as.data.frame(out)
      predicted_hosp <- out_df$I * all_params$hr1 + out_df$I2 * all_params$hr2
      
      # observed data 
      observed_hosp <- state_data$hosp
      
      # discard first few samples 
      
      #   pred_subset = predicted_hosp[11:length(time_vec)] 
      #   obs_subset =  observed_hosp[11:length(time_vec)] 
      #  
      # discard first few samples in case of initial dynamics 
      # discard first 15 time points when calculating likelihood
      start_idx <- 11  # start from the 16th element (discarding first 15)
      end_idx <- length(observed_hosp)
      
      if(start_idx >= end_idx || length(predicted_hosp) < end_idx || length(observed_hosp) < end_idx) {
        return(-1e6)
      }
      
      pred_subset <- pmax(predicted_hosp[start_idx:end_idx], 1e-6)
      obs_subset <- observed_hosp[start_idx:end_idx]
      
      if(any(is.na(pred_subset)) || any(is.na(obs_subset))) {
        return(-1e6)
      }
      
      # calculate likelihood using poisson distribution 
      state_log_lik <- sum(dpois(obs_subset, lambda = pred_subset, log = TRUE))
      
      # adjust 
      # state_log_lik = state_log_lik * unique(state_data$adjust)
      
      
      # iteratively add likelihoods from across states 
      total_log_lik <- total_log_lik + state_log_lik
    }
    
    print(total_log_lik)
    return(total_log_lik)
  }
  
  
  
  
  
  #   start_idx <- min(10, length(observed_hosp) - 50)
  #    end_idx <- length(observed_hosp)
  
  #    if(start_idx >= end_idx || length(predicted_hosp) < end_idx || length(observed_hosp) < end_idx) {
  #      return(-1e6)
  #    }
  
  #    pred_subset <- pmax(predicted_hosp[start_idx:end_idx], 1e-6)
  #    obs_subset <- observed_hosp[start_idx:end_idx]
  #    
  #    if(any(is.na(pred_subset)) || any(is.na(obs_subset))) {
  #      return(-1e6)
  #    }
  
  # calculate likelihood using poission distribution 
  #    state_log_lik <- sum(dpois(obs_subset, lambda = pred_subset, log = TRUE))
  
  # iteratively add likeliloods from across states 
  #  total_log_lik <- total_log_lik + state_log_lik
  #   }
  #    
  #    print(total_log_lik)
  #    return(total_log_lik)
  #  }
  
  # prior function 
  log_prior <- function(params) {
    
    # global parameters 
    global_bounds <- list(
    #  var = c(0.05, .09),
      var_scalar = c(.991, 3), 
      hum = c(0.00001, .30), 
      temp = c(0.0, 1.5),
      imm = c(1/50, 1/12),
      imm_vax = c(1/20, 1/6),
      hr1 = c(0.0001, 0.05),
      hr2 = c(0.0001, 0.05), 
      S0 = c(0.30, 0.80)
    )
    
    # state-specific parameters 
    state_bounds <- list(
      beta = c(0.1, 0.70)
    )
    
    # Check global parameters
    for(param_name in names(global_bounds)) {
      if(param_name %in% names(params)) {
        bounds <- global_bounds[[param_name]]
        if(params[param_name] < bounds[1] || params[param_name] > bounds[2]) {
          return(-Inf)
        }
      }
    }
    
    # Check state-specific parameters
    for(state_name in states) {
      # Check beta
      beta_param <- paste0("beta_", state_name)
      if(beta_param %in% names(params)) {
        bounds <- state_bounds[["beta"]]
        if(params[beta_param] < bounds[1] || params[beta_param] > bounds[2]) {
          return(-Inf)
        }
      }
      
    }
    
    return(0)  # Uniform priors within bounds
  }
  
  
  log_posterior <- function(params) {
    
    # prior 
    # likelihood of initial guess of data 
    lp <- log_prior(params)
    
    if(is.infinite(lp)) {
      return(-Inf)
    }
    
    # posterior - likelihood of model given data 
    ll <- likelihood_func(params)
    
    # sum them to get the posterior 
    return(lp + ll)
  }
  
  # initial parameter guesses 
  global_initial <- list(
    var_scalar = 1.7,
    #var = .07,
    hum = .11,   #0.10,
    temp =.05,   # 0.06,
    imm = 1/26,
    imm_vax = 1/8,
    hr1 = 0.008,
    hr2 = 0.003, 
    S0 = 0.34
  )
  
  # global_initial <- list(
  #   var = 0.03,
  #   hum = 0.14,
  #   temp = 0.22,
  #   imm = 1/24,
  #   imm_vax = 1/8,
  #   hr1 = 0.008,
  #   hr2 = 0.003, 
  #   S0 = 0.33
  # )
  
  #  global_initial <- list(
  #    var = 0.02,
  #    hum = 0.19,
  #    temp = 0.0001,
  #    imm = 1/24,
  #    imm_vax = 1/16,
  #    hr1 = 0.008,
  #    hr2 = 0.003, 
  #    S0 = 0.50
  #  )
  
  # state-specific initial values (only fitted parameters)
  state_initial <- list()
  for(state_name in states) {
    state_initial[[paste0("beta_", state_name)]] <- runif(1, 0.25, .35)
  }
  
  # combine initial parameters 
  initial_params <- c(global_initial, state_initial)
  param_names <- names(initial_params)
  
  cat("=== PARAMETER SUMMARY ===\n")
  cat("Total FITTED parameters:", length(param_names), "\n")
  cat("  Global fitted parameters:", length(global_initial), "\n")
  cat("  State-specific fitted parameters per state:", (length(initial_params) - length(global_initial)) / n_states, "\n")
  cat("    - beta (transmission rate): FITTED\n") 
  cat("    - S0 (initial susceptible fraction): FITTED\n")
  cat("Total FIXED parameters:", length(states), "\n")
  cat("  Population sizes (N): FIXED (not fitted)\n")
  cat("  Other fixed: lp, gamma, VE\n")
  cat("==========================\n\n")
  
  # Test the likelihood function with initial parameters
  cat("\n=== TESTING LIKELIHOOD FUNCTION ===\n")
  test_params <- unlist(initial_params)
  names(test_params) <- param_names
  test_ll <- likelihood_func(test_params)
  cat("Test likelihood with initial parameters:", test_ll, "\n")
  cat("=== END TEST ===\n\n")
  
  # Function to run a single MCMC chain with JOINT PROPOSALS AT FIXED 16%
  run_chain <- function(chain_id, iterations, burn_in, thin) {
    
    current_params <- unlist(initial_params)
    names(current_params) <- param_names
    
    # Perturb initial values for different chains
    if(chain_id > 1) {
      current_params <- current_params * runif(length(current_params), 0.90, 1.10)
    }
    
    samples <- matrix(NA, nrow = ceiling((iterations - burn_in) / thin), 
                      ncol = length(param_names))
    colnames(samples) <- param_names
    
    
    # posterior of current parametres 
    current_log_post <- log_posterior(current_params)
    cat(sprintf("Chain %d: Initial log posterior = %.2f\n", chain_id, current_log_post))
    
    if(!is.finite(current_log_post)) {
      cat(sprintf("Warning: Chain %d has non-finite initial log posterior. Adjusting initial values.\n", chain_id))
      attempts <- 0
      while(!is.finite(current_log_post) && attempts < 10) {
        current_params <- current_params * runif(length(current_params), 0.8, 1.2)
        current_log_post <- log_posterior(current_params)
        attempts <- attempts + 1
      }
      
      if(!is.finite(current_log_post)) {
        stop(sprintf("Chain %d: Cannot find valid initial parameters after 10 attempts", chain_id))
      }
      cat(sprintf("Chain %d: Found valid initial log posterior = %.2f after %d attempts\n", 
                  chain_id, current_log_post, attempts))
    }
    
    # FIXED 16% PROPOSAL: Initialize fixed covariance matrix
    n_params <- length(param_names)
    
    # Fixed diagonal covariance matrix using 16% of parameter values
    # this means: 
    
    #  Takes each parameter's current value
    #  Multiplies by 0.16 to get the standard deviation for proposals
    # Example: if a parameter = 0.5, then proposal std = 0.5 × 0.16 = 0.08
    
    fixed_scales <- abs(current_params) * 0.06  # FIXED 16% proposal
    fixed_scales[fixed_scales < 1e-5] <- 1e-4
    proposal_cov <- diag(fixed_scales^2)
    rownames(proposal_cov) <- param_names
    colnames(proposal_cov) <- param_names
    
    # Fixed Cholesky decomposition - no adaptation
    chol_cov <- chol(proposal_cov)
    
    # Initialize counters
    accept_count <- 0
    total_proposals <- 0
    
    sample_idx <- 1
    
    cat(sprintf("Chain %d: Using FIXED 16%% joint proposals (no adaptation)\n", chain_id))
    
    for(i in 1:iterations) {
      
      # JOINT PROPOSAL: Propose all parameters at once using FIXED 16% rate
      tryCatch({
        # Generate multivariate normal proposal with FIXED covariance
        z <- rnorm(n_params)
        proposal_increment <- as.vector(t(chol_cov) %*% z)
        proposed_params <- current_params + proposal_increment
        
        # Ensure proposed parameters have correct names
        names(proposed_params) <- param_names
        
        # Evaluate posterior at proposed point
        proposed_log_post <- log_posterior(proposed_params)
        
        total_proposals <- total_proposals + 1
        
        if (i %% 100 == 0) {
          cat(sprintf("Iter %d: Current log post = %.2f, Proposed = %.2f\n", 
                      i, current_log_post, proposed_log_post))
        }
        
        # Accept/reject decision
        if(!is.na(proposed_log_post) && is.finite(proposed_log_post)) {
          
          # posterior of new/ posterior of current 
          # log posterior of new - log posterior current 
          log_accept_ratio <- proposed_log_post - current_log_post
          
          # If proposal is better (log_accept_ratio > 0), always accept
          # If proposal is worse (log_accept_ratio < 0), accept with probability exp(log_accept_ratio)
          if(is.finite(log_accept_ratio) && log(runif(1)) < log_accept_ratio) {
            current_params <- proposed_params
            current_log_post <- proposed_log_post
            accept_count <- accept_count + 1
          }
        }
        
        # NO ADAPTATION - covariance matrix remains fixed at 16%
        
      }, error = function(e) {
        cat(sprintf("Error in iteration %d: %s\n", i, e$message))
        # Continue with current parameters
      })
      
      # Store samples after burn-in and thinning
      if(i > burn_in && (i - burn_in) %% thin == 0) {
        samples[sample_idx, ] <- current_params
        sample_idx <- sample_idx + 1
      }
      
      # print progress
      if(i %% 500 == 0) {
        overall_accept <- if(total_proposals > 0) accept_count / total_proposals else 0
        cat(sprintf("Chain %d, Iteration %d: Acceptance rate = %.3f (FIXED 16%% proposals)\n", 
                    chain_id, i, overall_accept))
      }
    }
    
    final_accept_rate <- if(total_proposals > 0) accept_count / total_proposals else 0
    
    return(list(
      samples = samples,
      acceptance_rate = final_accept_rate,
      final_proposal_cov = proposal_cov
    ))
  }
  
  # run multiple chains
  cat("Starting hierarchical MCMC with FIXED 16% JOINT PROPOSALS,", n_chains, "chains,", iterations, "iterations\n")
  
  if(n_chains > 1 && requireNamespace("parallel", quietly = TRUE)) {
    chain_results <- mclapply(
      1:n_chains, 
      function(id) run_chain(id, iterations, burn_in, thin),
      mc.cores = min(n_chains, detectCores() - 1)
    )
  } else {
    chain_results <- lapply(
      1:n_chains, 
      function(id) run_chain(id, iterations, burn_in, thin)
    )
  }
  
  # combine results
  all_samples <- list()
  acceptance_rates <- numeric(n_chains)
  
  for(i in 1:n_chains) {
    all_samples[[i]] <- chain_results[[i]]$samples
    acceptance_rates[i] <- chain_results[[i]]$acceptance_rate
  }
  
  # create MCMC objects
  mcmc_objects <- lapply(all_samples, mcmc)
  mcmc_list <- mcmc.list(mcmc_objects)
  
  # diagnostics - with error handling
  gelman_diag <- NULL
  if(n_chains > 1) {
    tryCatch({
      valid_samples <- TRUE
      for(i in 1:n_chains) {
        chain_samples <- all_samples[[i]]
        if(any(!is.finite(chain_samples)) || any(apply(chain_samples, 2, var, na.rm = TRUE) == 0)) {
          valid_samples <- FALSE
          break
        }
      }
      
      if(valid_samples) {
        gelman_diag <- gelman.diag(mcmc_list, autoburnin = FALSE)
        cat("Gelman-Rubin diagnostics calculated successfully\n")
      } else {
        cat("Warning: Cannot calculate Gelman-Rubin diagnostics due to numerical issues\n")
      }
    }, error = function(e) {
      cat("Warning: Could not calculate Gelman-Rubin diagnostics:", e$message, "\n")
    })
  }
  
  summary_stats <- summary(mcmc_list)
  
  return(list(
    samples = all_samples,
    mcmc_list = mcmc_list,
    acceptance_rates = acceptance_rates,
    gelman_diag = gelman_diag,
    summary = summary_stats,
    states = states,
    param_names = param_names,
    final_proposal_covs = lapply(chain_results, function(x) x$final_proposal_cov)
  ))
}




#############################################################
### Run the hierarchical MCMC calibration 
results <- run_hierarchical_mcmc_calibration(
  full_calibration_dat = full_calibration_dat,  # your data
  iterations = 700,     # start small for testing
  burn_in = 300,         # burn-in period
  thin = 1,              # thinning interval
  n_chains = 1 )          # number of MCMC chains


#############################################################
### Post-processing MCMC results 

sink("mcmc_summary_sep19.txt")


# Look at convergence diagnostics
if(!is.null(results$gelman_diag)) {
  print("Gelman-Rubin diagnostics (values close to 1.0 indicate convergence):")
  print(results$gelman_diag)
}

# Look at acceptance rates
print("Acceptance rates by chain:")
print(results$acceptance_rates)

# Summary statistics
print("Parameter estimates:")
print(results$summary$statistics)

# Close the sink
sink()
dev.off()

# Plot trace plots and posterior distributions
pdf("trace_posterior_plots_sep19.pdf", width = 14, height = 20)  

mcmc_df <- as.data.frame(as.matrix(results$mcmc_list))

# Add iteration index
mcmc_df$iteration <- 1:nrow(mcmc_df)

# Make ggplot traces
temp_plot <- ggplot(mcmc_df, aes(x = iteration, y = temp)) +
  geom_line() +
  ggtitle("Trace plot for 'temp' parameter") +
  theme_minimal()
temp_plot

hum_plot <- ggplot(mcmc_df, aes(x = iteration, y = hum)) +
  geom_line() +
  ggtitle("Trace plot for 'hum' parameter") +
  theme_minimal()

hum_plot

imm_plot <- ggplot(mcmc_df, aes(x = iteration, y = imm)) +
  geom_line() +
  ggtitle("Trace plot for 'imm' parameter") +
  theme_minimal()

imm_plot

imm_vax_plot <- ggplot(mcmc_df, aes(x = iteration, y = imm_vax)) +
  geom_line() +
  ggtitle("Trace plot for 'imm' parameter") +
  theme_minimal()
imm_vax_plot

s0_plot <- ggplot(mcmc_df, aes(x = iteration, y = S0)) +
  geom_line() +
  ggtitle("Trace plot for 'SO' parameter") +
  theme_minimal()
s0_plot


global_plot = plot_grid(temp_plot, hum_plot,  imm_plot, imm_vax_plot, s0_plot)
global_plot 


ggsave("global_trace_plots.png", plot = global_plot, width = 14, height = 10, dpi = 300)


# Plot trace plots for key parameters (optional)
# You can plot traces to check convergence visually
if(requireNamespace("coda", quietly = TRUE)) {
  
  # Plot trace for a state-specific parameter
  first_state <- results$states[1]
  beta_param_name <- paste0("beta_", first_state)
  if(beta_param_name %in% colnames(results$mcmc_list[[1]])) {
    plot(results$mcmc_list[, beta_param_name], 
         main = paste("Trace plot for", beta_param_name))
  }
  
  ##  second_state <- results$states[2]
  #  beta_param_name <- paste0("beta_", second_state)
  #  if(beta_param_name %in% colnames(results$mcmc_list[[1]])) {
  #    plot(results$mcmc_list[, beta_param_name], 
  #         main = paste("Trace plot for", beta_param_name))
  #  }
  
  
}

dev.off() 

##############################
# loop to plot all outputs 
model_outputs = left_join(full_calibration_dat, census_dat, by = c("Location", "pop_size_2023"))
head(model_outputs)

tab = data.frame(results$summary$statistics)
head(tab)

# try to add more data 
#snippet = model_outputs %>%
#  filter(t < 31)

#model_outputs = left_join(full_calibration_dat, census_dat, by = c("Location", "pop_size_2023")) %>%
#  group_by(Location) %>%
#  mutate(t = t + 30) %>%
#  ungroup()

#model_outputs = rbind(snippet, model_outputs)



# initalize plots

plot_list <- list()

state = print(unique(model_outputs$Location))

for(i in state){
  
  dat = model_outputs %>%
    filter(Location == i)
  
  tab_state <- results$summary$statistics[grepl(paste0(i, "$"), rownames(results$summary$statistics)), ] %>%
    as.data.frame()
  
  
  # states 
  N <- print(unique(dat$pop_size_2023))
  
  # S0 <- .13*N
  # S0 <- (tab_state[2][i]) **N
  # S0 <- as.numeric(tab_state[grepl("S0_", rownames(tab_state)), "Mean"]) * N

   #  S0 =  results$summary$statistics[8] * N
   S0 = .36*N
  # print(S0/N)
#  S0 =  .50* N
  E0 <- 0.01 *N
  I0 <- 0.001 * N
  R0 <- .94*N - S0 
  V0 <- 0.049 * N
  
  S20 = 0.00 * N
  E20 <- 0.00*N 
  I20 <- 0.00* N
  R20 <- 0.00 * N
  V20 <- 0.00 * N
  
  state <- c(S = S0, E = E0, I = I0, R = R0, V = V0, S2 = S20, E2 = E20, I2 = I20, R2 = R20, V2 = V20)
  
  # Set time for testing simulation
  time = seq(1, nrow( dat), 1)
  
  print(results$param_names)
  # fit 
  param <- c(
    var_scalar = 1.34,    #results$summary$statistics[1],   
    hum = .13,  # results$summary$statistics[2], 
    temp = .056, # results$summary$statistics[3],  
    imm = .049,  #results$summary$statistics[4], 
    imm_vax = .076, # = results$summary$statistics[5],
    hr = .0097, # results$summary$statistics[6],
    hr2 = .0036, #results$summary$statistics[7],
     beta = .43,
   # beta = tab_state["Mean", 1],                       
    #  beta = as.numeric(tab_state[grepl("beta_", rownames(tab_state)), "Mean"]),
    lp = 1/(4/7),                    # latent period 
    gamma = 1/(4/7),              
    VE = .30, 
    N =  print(unique(dat$pop_size_2023))
    
  )
  
  print(param)
  
  seirs_model_hybrid_state = function(t, state, param) {
    with(as.list(state), {
      
      # parameters for seasonality multipliers 
      var = param["var"]
      hum = param["hum"]
      temp = param["temp"]
      
      # parameters for underlying dynamics 
      beta = param["beta"]
      lp = param["lp"]           # latent period
      gamma = param["gamma"]     # infectious period 
      imm = param["imm"]         # immunity from infection 1
      imm_vax = param["imm_vax"] # immunity from vaccination 
      N = param["N"]             # pop
      VE = param["VE"]           # VE
      hr1 = param["hr"]           # hosp rate 1 
      hr2 = param["hr2"]           # hosp rate 1 
      
      # Empirical time series 
      vax_rate_at_t = dat$incident[t] 
      humidity_at_t = dat$humid_smooth[t] 
      temp_at_t = dat$inverse_temp[t]
    #  variant_at_t = dat$odds_frequency[t] 
      
      dat = dat %>%
        mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, 1.24, 1))
      
      imm_vec = imm  * dat$immunity_scalar_fitted[t]
      
      # Seasonal function 
      seasonal = beta*( hum*humidity_at_t + temp*temp_at_t  ) 
      # .00001* humidity_at_t*temp_at_t)
      
      dS = -seasonal*S*(I + I2)/N - VE*vax_rate_at_t*S + imm_vax*V
      dE =  seasonal*S*(I+ I2)/N - lp*E
      dI =  lp*E - gamma*I
      dR = I*gamma -  R*imm_vec  
      dV = VE*vax_rate_at_t*S  - imm_vax*V   
      
      dS2 = -seasonal*S2*(I + I2)/N  + imm_vec*R  + imm_vec*R2 + imm_vax*V2 - VE*vax_rate_at_t*S2
      dE2 =  seasonal*S2*(I + I2)/N - lp*E2
      dI2 =  lp*E2 - gamma*I2
      dR2 = I2*gamma -  R2*imm_vec
      dV2 = VE*vax_rate_at_t*S2  - imm_vax*V2  
      
      # Return the rates of change
      list(c(dS, dE, dI, dR, dV, dS2, dE2, dI2, dR2, dV2))
    })
  }
  
  # Solve the ODE system
  out_check <- ode(y = state, times = time, func = seirs_model_hybrid_state, parms = param)
  
  # hop
  out_df_check <- as.data.frame(out_check) %>%
    mutate(predicted_hosp = I* .0097 + I2* .0036 )
  
  print(out_df_check$predicted_hosp)
  
  # plot 
  plot =  ggplot() +
    geom_line(data = out_df_check %>% filter(time > 12), aes(x = time, y = predicted_hosp), col = "red",lwd = 1.5) +
    geom_line(data = dat %>% filter(t > 12), aes(x = t, y = rollmean(hosp, k = 3, na.pad = TRUE)), col = "black", lwd = 1.5) +
    ggtitle(paste("State:", i)) +
    ylab("Hospitalizations") +
    theme_bw()
  
  plot_list[[i]] <- plot
  
}



ny = plot_list[["NY"]]
ct = plot_list[["CT"]]

mn = plot_list[["MN"]]
oh = plot_list[["OH"]]
mi = plot_list[["MI"]]


tn = plot_list[["TN"]]
ga = plot_list[["GA"]]

co = plot_list[["CO"]]

plot_grid(ny, mn, oh, mi, tn , ga, ct, co)









# initalize plots

plot_list <- list()

state = print(unique(model_outputs$Location))

for(i in state){
  
  dat = model_outputs %>%
    filter(Location == i)  %>%
    mutate(immunity_scalar_fitted = ifelse(immunity_scalar > 1.005, immunity_scalar*results$summary$statistics[1]*2, immunity_scalar))
  
  #imm_vec = imm* state_data1$immunity_scalar_fitted[t_idx] #%>%
  # filter(t > 2)
  
  tab_state <- results$summary$statistics[grepl(paste0(i, "$"), rownames(results$summary$statistics)), ] %>%
    as.data.frame()
  
  # print(head(tab_state))
  
  # states 
  N <- print(unique(dat$pop_size_2023))
  
  # S0 <- .13*N
  # S0 <- (tab_state[2][i]) **N
  # S0 <- as.numeric(tab_state[grepl("S0_", rownames(tab_state)), "Mean"]) * N
   S0 =  results$summary$statistics[8] * N
  #  S0 = .37*N
  # print(S0/N)
 # S0 =  .50* N
  E0 <- 0.01 *N
  I0 <- 0.001 * N
  R0 <- .94*N - S0 
  V0 <- 0.049 * N
  
  S20 = 0.00 * N
  E20 <- 0.00*N 
  I20 <- 0.00* N
  R20 <- 0.00 * N
  V20 <- 0.00 * N
  
  state <- c(S = S0, E = E0, I = I0, R = R0, V = V0, S2 = S20, E2 = E20, I2 = I20, R2 = R20, V2 = V20)
  
  # Set time for testing simulation
  time = seq(1, nrow( dat), 1)
  
  print(results$param_names)
  # fit 
  param <- c(
     var_scalar = results$summary$statistics[1],   
    hum =   results$summary$statistics[2], 
    temp =  results$summary$statistics[3],  
    imm = results$summary$statistics[4], 
    imm_vax  = results$summary$statistics[5],
    hr =  results$summary$statistics[6],
    hr2 = results$summary$statistics[7],
  #  beta = .31,
     beta = tab_state["Mean", 1],                       
    #  beta = as.numeric(tab_state[grepl("beta_", rownames(tab_state)), "Mean"]),
    lp = 1/(4/7),                    # latent period 
    gamma = 1/(4/7),              
    VE = .30, 
    N =  print(unique(dat$pop_size_2023))
    
  )
  
  #param <- c(
  #    var = 0.07,   
  #    hum = .099, 
  #    temp = .01,  
  #    imm =.0537,     #results$summary$statistics[4], 
  #    imm_vax  = .11789 , # results$summary$statistics[5],
  #    hr =  results$summary$statistics[6], 
  #    hr2 = results$summary$statistics[7],
  #    beta =   .80,                       
  #    #  beta = as.numeric(tab_state[grepl("beta_", rownames(tab_state)), "Mean"]),
  #    lp = 1/(4/7),                    # latent period 
  #    gamma = 1/(4/7),              
  #    VE = .30, # 
  #    N =  print(unique(dat$pop_size_2023))
  #    )
  
  print(param)
  
  seirs_model_hybrid_state = function(t, state, param) {
    with(as.list(state), {
      
      # parameters for seasonality multipliers 
      var = param["var"]
      hum = param["hum"]
      temp = param["temp"]
      
      # parameters for underlying dynamics 
      beta = param["beta"]
      lp = param["lp"]           # latent period
      gamma = param["gamma"]     # infectious period 
      imm = param["imm"]         # immunity from infection 1
      imm_vax = param["imm_vax"] # immunity from vaccination 
      N = param["N"]             # pop
      VE = param["VE"]           # VE
      hr1 = param["hr"]           # hosp rate 1 
      hr2 = param["hr2"]           # hosp rate 1 
      
      # Empirical time series 
      vax_rate_at_t = dat$incident[t] 
      humidity_at_t = dat$humid_smooth[t] 
      temp_at_t = dat$inverse_temp[t]
      #  variant_at_t = dat$odds_frequency[t] 
      
    #  imm_vec = imm  #* dat$immunity_scalar[t]
      
      imm_vec = imm #* dat$immunity_scalar_fitted[t]
      
      # Seasonal function 
      seasonal = beta*( hum*humidity_at_t + temp*temp_at_t  ) 
      # .00001* humidity_at_t*temp_at_t)
      
      dS = -seasonal*S*(I + I2)/N - VE*vax_rate_at_t*S + imm_vax*V
      dE =  seasonal*S*(I+ I2)/N - lp*E
      dI =  lp*E - gamma*I
      dR = I*gamma -  R*imm_vec  
      dV = VE*vax_rate_at_t*S  - imm_vax*V   
      
      dS2 = -seasonal*S2*(I + I2)/N  + imm_vec*R  + imm_vec*R2 + imm_vax*V2 - VE*vax_rate_at_t*S2
      dE2 =  seasonal*S2*(I + I2)/N - lp*E2
      dI2 =  lp*E2 - gamma*I2
      dR2 = I2*gamma -  R2*imm_vec
      dV2 = VE*vax_rate_at_t*S2  - imm_vax*V2  
      
      # Return the rates of change
      list(c(dS, dE, dI, dR, dV, dS2, dE2, dI2, dR2, dV2))
    })
  }
  
  # Solve the ODE system
  out_check <- ode(y = state, times = time, func = seirs_model_hybrid_state, parms = param)
  
  # hop
  out_df_check <- as.data.frame(out_check) %>%
    mutate(predicted_hosp = I* results$summary$statistics[6]+ I2* results$summary$statistics[7] )
  
  print(out_df_check$predicted_hosp)
  
  # plot 
  plot =  ggplot() +
    geom_line(data = out_df_check %>% filter(time > 12), aes(x = time, y = predicted_hosp), col = "red",lwd = 1.4) +
    geom_line(data = dat %>% filter(t > 12), aes(x = t, y = rollmean(hosp, k = 2, na.pad = TRUE)), col = "black", lwd = 1.4) +
    ggtitle(paste("State:", i)) +
    ylab("Hospitalizations") +
    theme_bw()
  
  plot_list[[i]] <- plot
  
}

#     state_data1 = state_data %>%
#mutate(immunity_scalar_fitted = ifelse(immunity_scalar > 1.005, immunity_scalar*var_scalar, immunity_scalar))

#imm_vec = imm* state_data1$immunity_scalar_fitted[t_idx] 


ny = plot_list[["NY"]]
ct = plot_list[["CT"]]

mn = plot_list[["MN"]]
oh = plot_list[["OH"]]
mi = plot_list[["MI"]]


tn = plot_list[["TN"]]
ga = plot_list[["GA"]]

plot_grid(ny, mn, oh, mi, tn , ga, ct)

ca = plot_list[["CA"]]
or = plot_list[["OR"]]
ut = plot_list[["UT"]]
nm = plot_list[["NM"]]
co = plot_list[["CO"]]
plot_grid(or, ut, nm, co)

plot_grid(ca, or, ut, nm, co)

plot_grid(ca, or, ut)

dev.off()

mn = plot_list[["MN"]]
oh = plot_list[["OH"]]
mi = plot_list[["MI"]]

plot_grid(mn, oh, mi, nrow = 1)



#############################################################
# based on fit, find best S0 and beta for a given state 

library(lhs)

n <- 1000   # number of samples
k <- 4     # number of parameters
U <- randomLHS(n, k)

sus_frac <- qunif(U[,1], min = 0.26, max = .46)
beta <- qunif(U[,2], min = 0.29, max = 0.45)
h1 <- qunif(U[,3], min = 0.0001, max = 0.012)
h2 <- qunif(U[,4], min = 0.0008, max = 0.004)

print(U)

# states 
dat = full_calibration_dat %>%
  filter(Location == "GA") %>%
  filter(t > 16) %>%
  mutate(t = seq(1, nrow(dat), 1))
N <- print(unique(dat$pop_size_2023))

store_RMSE <- numeric(n)

for(i in seq_len(n)){

S0 = sus_frac[i]*N
E0 <- 0.01 *N
I0 <- 0.001 * N
R0 <- .94*N - S0 
V0 <- 0.049 * N

S20 = 0.00 * N
E20 <- 0.00*N 
I20 <- 0.00* N
R20 <- 0.00 * N
V20 <- 0.00 * N

state <- c(S = S0, E = E0, I = I0, R = R0, V = V0, S2 = S20, E2 = E20, I2 = I20, R2 = R20, V2 = V20)

# Set time for testing simulation
time = seq(1, nrow( dat), 1)

# fit 
param <- c(
  # var = results$summary$statistics[1],   
  hum =   results$summary$statistics[1], 
  temp =  results$summary$statistics[2],  
  imm = results$summary$statistics[3], 
  imm_vax  = results$summary$statistics[4],
  hr = h1[i],
  hr2 =  h2[i],
  beta = beta[i],
  lp = 1/(4/7),                    # latent period 
  gamma = 1/(4/7),              
  VE = .30, 
  N =  unique(dat$pop_size_2023)
  
)

seirs_model_hybrid_state = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    var = param["var"]
    hum = param["hum"]
    temp = param["temp"]
    
    # parameters for underlying dynamics 
    beta = param["beta"]
    lp = param["lp"]           # latent period
    gamma = param["gamma"]     # infectious period 
    imm = param["imm"]         # immunity from infection 1
    imm_vax = param["imm_vax"] # immunity from vaccination 
    N = param["N"]             # pop
    VE = param["VE"]           # VE
    hr1 = param["hr"]           # hosp rate 1 
    hr2 = param["hr2"]           # hosp rate 1 
    
    # Empirical time series 
    vax_rate_at_t = dat$incident[t] 
    humidity_at_t = dat$humid_smooth[t] 
    temp_at_t = dat$inverse_temp[t]
    #  variant_at_t = dat$odds_frequency[t] 
    
    imm_vec = imm  #* dat$immunity_scalar[t]
    
    # Seasonal function 
    seasonal = beta*( hum*humidity_at_t + temp*temp_at_t  ) 
    # .00001* humidity_at_t*temp_at_t)
    
    dS = -seasonal*S*(I + I2)/N - VE*vax_rate_at_t*S + imm_vax*V
    dE =  seasonal*S*(I+ I2)/N - lp*E
    dI =  lp*E - gamma*I
    dR = I*gamma -  R*imm_vec  
    dV = VE*vax_rate_at_t*S  - imm_vax*V   
    
    dS2 = -seasonal*S2*(I + I2)/N  + imm_vec*R  + imm_vec*R2 + imm_vax*V2 - VE*vax_rate_at_t*S2
    dE2 =  seasonal*S2*(I + I2)/N - lp*E2
    dI2 =  lp*E2 - gamma*I2
    dR2 = I2*gamma -  R2*imm_vec
    dV2 = VE*vax_rate_at_t*S2  - imm_vax*V2  
    
    # Return the rates of change
    list(c(dS, dE, dI, dR, dV, dS2, dE2, dI2, dR2, dV2))
  })
}

# Solve the ODE system
out_check <- ode(y = state, times = time, func = seirs_model_hybrid_state, parms = param)

# hop
out_df_check <- as.data.frame(out_check) %>%
  mutate(predicted_hosp = I* h1[i]+ I2* h2[i]) %>%
  mutate(t = time) %>%
 # filter(t > 11) %>%
  left_join(dat, by = "t")

RMSE <- sqrt(sum((out_df_check$predicted_hosp - out_df_check$hosp)^2, na.rm = TRUE))
store_RMSE[i] <- RMSE

} 

head(store_RMSE)
performance <- data.frame(
  sample = seq_len(n),
  RMSE = store_RMSE) %>%
  arrange(RMSE) 


# plot 10 
performance_top10 = performance[1:12,]

samples = performance_top10$sample
print(samples)

params_df <- data.frame(
  sample = numeric(0),
  sus_frac = numeric(0),
  beta = numeric(0),
  h1 = numeric(0),
  h2 = numeric(0), 
  RMSE = numeric(0)
)


plot_top10 = list()

for (k in seq_along(samples)) {
  j <- samples[k]
  
fit1= sus_frac[j] # .33
fit2 = beta[j]
fit3 = h1[j]
fit4 = h2[j]


# plot it

S0 = fit1*N
E0 <- 0.01 *N
I0 <- 0.001 * N
R0 <- .94*N - S0 
V0 <- 0.049 * N

S20 = 0.00 * N
E20 <- 0.00*N 
I20 <- 0.00* N
R20 <- 0.00 * N
V20 <- 0.00 * N

state <- c(S = S0, E = E0, I = I0, R = R0, V = V0, S2 = S20, E2 = E20, I2 = I20, R2 = R20, V2 = V20)

# Set time for testing simulation
time = seq(1, nrow( dat), 1)

# fit 
param <- c(
  # var = results$summary$statistics[1],   
  hum =   results$summary$statistics[1], 
  temp = results$summary$statistics[2],  
  imm = results$summary$statistics[3], 
  imm_vax  = results$summary$statistics[4],
  hr = fit3,
  hr2 =  fit4,
  beta = fit2,
  lp = 1/(4/7),                    # latent period 
  gamma = 1/(4/7),              
  VE = .30, 
  N =  unique(dat$pop_size_2023)
  
)

seirs_model_hybrid_state = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    var = param["var"]
    hum = param["hum"]
    temp = param["temp"]
    
    # parameters for underlying dynamics 
    beta = param["beta"]
    lp = param["lp"]           # latent period
    gamma = param["gamma"]     # infectious period 
    imm = param["imm"]         # immunity from infection 1
    imm_vax = param["imm_vax"] # immunity from vaccination 
    N = param["N"]             # pop
    VE = param["VE"]           # VE
    hr1 = param["hr"]           # hosp rate 1 
    hr2 = param["hr2"]           # hosp rate 1 
    
    # Empirical time series 
    vax_rate_at_t = dat$incident[t] 
    humidity_at_t = dat$humid_smooth[t] 
    temp_at_t = dat$inverse_temp[t]
    #  variant_at_t = dat$odds_frequency[t] 
    
    imm_vec = imm  #* dat$immunity_scalar[t]
    
    # Seasonal function 
    seasonal = beta*( hum*humidity_at_t + temp*temp_at_t  ) 
    # .00001* humidity_at_t*temp_at_t)
    
    dS = -seasonal*S*(I + I2)/N - VE*vax_rate_at_t*S + imm_vax*V
    dE =  seasonal*S*(I+ I2)/N - lp*E
    dI =  lp*E - gamma*I
    dR = I*gamma -  R*imm_vec  
    dV = VE*vax_rate_at_t*S  - imm_vax*V   
    
    dS2 = -seasonal*S2*(I + I2)/N  + imm_vec*R  + imm_vec*R2 + imm_vax*V2 - VE*vax_rate_at_t*S2
    dE2 =  seasonal*S2*(I + I2)/N - lp*E2
    dI2 =  lp*E2 - gamma*I2
    dR2 = I2*gamma -  R2*imm_vec
    dV2 = VE*vax_rate_at_t*S2  - imm_vax*V2  
    
    # Return the rates of change
    list(c(dS, dE, dI, dR, dV, dS2, dE2, dI2, dR2, dV2))
  })
}

# Solve the ODE system
out_check <- ode(y = state, times = time, func = seirs_model_hybrid_state, parms = param)

# hop
out_df_check_state <- as.data.frame(out_check) %>%
  mutate(predicted_hosp = I* fit3+ I2* fit4) %>%
  mutate(t = time) %>%
 # filter(t > 11) %>%
  left_join(dat, by = "t")
head(out_df_check_state)

rmse_val <- sqrt(mean((out_df_check_state$predicted_hosp - out_df_check_state$hosp)^2, na.rm=TRUE))

params_df <- rbind(params_df, data.frame(
  sample = j,
  sus_frac = fit1,
  beta = fit2,
  h1 = fit3,
  h2 = fit4,
  RMSE = rmse_val
))

plot_top = ggplot(data = out_df_check_state) +
 # geom_line(aes(x = date, y = predicted_hosp), col = "aquamarine3", lwd =2.5) +
  geom_line(aes(x = date, y = predicted_hosp), col = "mediumpurple1", lwd =1.75) +
  geom_line(aes(x = date, y = hosp), col = "black", lwd = 1.75) +
  ggtitle(unique(out_df_check_state$Location)) +
  theme_bw() +
  ylab("Hospitalizations") + xlab('Date')
 print(plot_top)

plot_top10[[k]] = plot_top


}
  

plot_all10 <- wrap_plots(plot_top10, ncol=3)
plot_all10

params_df

##################################################


tn = plot_list[["TN"]]
ga = plot_list[["GA"]]
plot_grid(ga, tn )

ny = plot_list[["NY"]]
ct = plot_list[["CT"]]
plot_grid(ny, ct )


print(plot_list[["MN"]])
print(plot_list[["OH"]])
print(plot_list[["MI"]])

library(zoo)


print(plot_list[["CT"]])
print(plot_list[["MN"]])
print(plot_list[["OH"]])


pdf("state_fit_plots_sep19.pdf", width = 24, height =8)  

#wrap_plots(plotlist = plot_list, ncol = 4)

# Convert list to regular list of plots
plots <- plot_list

# Display plots in a grid of 4 columns (adjust number of rows as needed)
do.call(grid.arrange, c(plots, ncol = 4))


##############################################################
# try with variant scalar 

# based on fit, find best S0 and beta for a given state 

library(lhs)

n <- 150   # number of samples
k <- 4     # number of parameters
U <- randomLHS(n, k)

sus_frac <- qunif(U[,1], min = 0.36, max = .46)
beta <- qunif(U[,2], min = 0.38, max = 0.46)
h1 <- qunif(U[,3], min = 0.0001, max = 0.012)
h2 <- qunif(U[,4], min = 0.0008, max = 0.004)

print(U)

# states 
dat = full_calibration_dat %>%
  filter(Location == "CT") %>%
  filter(t > 12) %>%
  mutate(t = seq(1, nrow(dat), 1))
N <- print(unique(dat$pop_size_2023))

store_RMSE <- numeric(n)

for(i in seq_len(n)){
  
  S0 = sus_frac[i]*N
  E0 <- 0.01 *N
  I0 <- 0.001 * N
  R0 <- .94*N - S0 
  V0 <- 0.049 * N
  
  S20 = 0.00 * N
  E20 <- 0.00*N 
  I20 <- 0.00* N
  R20 <- 0.00 * N
  V20 <- 0.00 * N
  
  state <- c(S = S0, E = E0, I = I0, R = R0, V = V0, S2 = S20, E2 = E20, I2 = I20, R2 = R20, V2 = V20)
  
  # Set time for testing simulation
  time = seq(1, nrow( dat), 1)
  
  # fit 
#  param <- c(
#     var_scalar = results$summary$statistics[1],   
#    hum =   results$summary$statistics[2], 
#    temp =  results$summary$statistics[3],  
#    imm = results$summary$statistics[4], 
#    imm_vax  = results$summary$statistics[5],
#    hr = h1[i],
#    hr2 =  h2[i],
#    beta = beta[i],
#    lp = 1/(4/7),                    # latent period 
#    gamma = 1/(4/7),              
#    VE = .30, 
#    N =  unique(dat$pop_size_2023)
    
#  )
  
  # fit 
  param <- c(
    var_scalar = 1.24,    #results$summary$statistics[1],   
    hum = .13,  # results$summary$statistics[2], 
    temp = .056, # results$summary$statistics[3],  
    imm = .049,  #results$summary$statistics[4], 
    imm_vax = .076, # = results$summary$statistics[5],
    hr = .0097, # results$summary$statistics[6],
    hr2 = .0036, #results$summary$statistics[7],
    beta = .43,
    # beta = tab_state["Mean", 1],                       
    #  beta = as.numeric(tab_state[grepl("beta_", rownames(tab_state)), "Mean"]),
    lp = 1/(4/7),                    # latent period 
    gamma = 1/(4/7),              
    VE = .30, 
    N =  print(unique(dat$pop_size_2023)))
    
  
  seirs_model_hybrid_state = function(t, state, param) {
    with(as.list(state), {
      
      # parameters for seasonality multipliers 
      var = param["var"]
      hum = param["hum"]
      temp = param["temp"]
      
      # parameters for underlying dynamics 
      beta = param["beta"]
      lp = param["lp"]           # latent period
      gamma = param["gamma"]     # infectious period 
      imm = param["imm"]         # immunity from infection 1
      imm_vax = param["imm_vax"] # immunity from vaccination 
      N = param["N"]             # pop
      VE = param["VE"]           # VE
      hr1 = param["hr"]           # hosp rate 1 
      hr2 = param["hr2"]           # hosp rate 1 
      
      # Empirical time series 
      vax_rate_at_t = dat$incident[t] 
      humidity_at_t = dat$humid_smooth[t] 
      temp_at_t = dat$inverse_temp[t]
      #  variant_at_t = dat$odds_frequency[t] 
      
   #   dat = dat %>%
   #     mutate(immunity_scalar_fitted = ifelse(infusion > 1.005, infusion*results$summary$statistics[1], infusion))
      
  #    imm_vec = imm  * dat$immunity_scalar_fitted[t]
      
      dat = dat %>%
        mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, 1.24, 1))
      
      imm_vec = imm  * dat$immunity_scalar_fitted[t]
      
      
      # Seasonal function 
      seasonal = beta*( hum*humidity_at_t + temp*temp_at_t  ) 
      # .00001* humidity_at_t*temp_at_t)
      
      dS = -seasonal*S*(I + I2)/N - VE*vax_rate_at_t*S + imm_vax*V
      dE =  seasonal*S*(I+ I2)/N - lp*E
      dI =  lp*E - gamma*I
      dR = I*gamma -  R*imm_vec  
      dV = VE*vax_rate_at_t*S  - imm_vax*V   
      
      dS2 = -seasonal*S2*(I + I2)/N  + imm_vec*R  + imm_vec*R2 + imm_vax*V2 - VE*vax_rate_at_t*S2
      dE2 =  seasonal*S2*(I + I2)/N - lp*E2
      dI2 =  lp*E2 - gamma*I2
      dR2 = I2*gamma -  R2*imm_vec
      dV2 = VE*vax_rate_at_t*S2  - imm_vax*V2  
      
      # Return the rates of change
      list(c(dS, dE, dI, dR, dV, dS2, dE2, dI2, dR2, dV2))
    })
  }
  
  # Solve the ODE system
  out_check <- ode(y = state, times = time, func = seirs_model_hybrid_state, parms = param)
  
  # hop
  out_df_check <- as.data.frame(out_check) %>%
    mutate(predicted_hosp = I* h1[i]+ I2* h2[i]) %>%
    mutate(t = time) %>%
    # filter(t > 11) %>%
    left_join(dat, by = "t")
  
  RMSE <- sqrt(sum((out_df_check$predicted_hosp - out_df_check$hosp)^2, na.rm = TRUE))
  store_RMSE[i] <- RMSE
  
} 

head(store_RMSE)
performance <- data.frame(
  sample = seq_len(n),
  RMSE = store_RMSE) %>%
  arrange(RMSE) 

# plot 10 
performance_top10 = performance[1:12,]

samples = performance_top10$sample
print(samples)

params_df <- data.frame(
  sample = numeric(0),
  sus_frac = numeric(0),
  beta = numeric(0),
  h1 = numeric(0),
  h2 = numeric(0), 
  RMSE = numeric(0)
)




plot_top10 = list()

for (k in seq_along(samples)) {
  j <- samples[k]
  
  fit1= sus_frac[j] # .33
  fit2 = beta[j]
  fit3 = h1[j]
  fit4 = h2[j]
  
  
  # plot it
  
 # S0 = fit1*N
  S0 = .26*N 
  E0 <- 0.01 *N
  I0 <- 0.001 * N
  R0 <- .94*N - S0 
  V0 <- 0.049 * N
  
  S20 = 0.00 * N
  E20 <- 0.00*N 
  I20 <- 0.00* N
  R20 <- 0.00 * N
  V20 <- 0.00 * N
  
  state <- c(S = S0, E = E0, I = I0, R = R0, V = V0, S2 = S20, E2 = E20, I2 = I20, R2 = R20, V2 = V20)
  
  # Set time for testing simulation
  time = seq(1, nrow( dat), 1)
  
  # fit 
 # param <- c(
#     var_scalar = results$summary$statistics[1],   
#    hum =   results$summary$statistics[2], 
#    temp = results$summary$statistics[3],  
#    imm = results$summary$statistics[4], 
#    imm_vax  = results$summary$statistics[5],
#    hr = fit3,
#    hr2 =  fit4,
#    beta = fit2,
#    lp = 1/(4/7),                    # latent period 
#    gamma = 1/(4/7),              
#    VE = .30, 
#    N =  unique(dat$pop_size_2023)
    
#  )
  
  param <- c(
    var_scalar = 1.24,    #results$summary$statistics[1],   
    hum = .13,  # results$summary$statistics[2], 
    temp = .056, # results$summary$statistics[3],  
    imm = .049,  #results$summary$statistics[4], 
    imm_vax = .076, # = results$summary$statistics[5],
    hr = .0097, # results$summary$statistics[6],
    hr2 = .0036, #results$summary$statistics[7],
    beta = .54,   #fit2,
    # beta = tab_state["Mean", 1],                       
    #  beta = as.numeric(tab_state[grepl("beta_", rownames(tab_state)), "Mean"]),
    lp = 1/(4/7),                    # latent period 
    gamma = 1/(4/7),              
    VE = .30, 
    N =  print(unique(dat$pop_size_2023)))
  
  
  seirs_model_hybrid_state = function(t, state, param) {
    with(as.list(state), {
      
      # parameters for seasonality multipliers 
      var = param["var"]
      hum = param["hum"]
      temp = param["temp"]
      
      # parameters for underlying dynamics 
      beta = param["beta"]
      lp = param["lp"]           # latent period
      gamma = param["gamma"]     # infectious period 
      imm = param["imm"]         # immunity from infection 1
      imm_vax = param["imm_vax"] # immunity from vaccination 
      N = param["N"]             # pop
      VE = param["VE"]           # VE
      hr1 = param["hr"]           # hosp rate 1 
      hr2 = param["hr2"]           # hosp rate 1 
      
      # Empirical time series 
      vax_rate_at_t = dat$incident[t] 
      humidity_at_t = dat$humid_smooth[t] 
      temp_at_t = dat$inverse_temp[t]
      #  variant_at_t = dat$odds_frequency[t] 
      
    #  dat = dat %>%
    #    mutate(immunity_scalar_fitted = ifelse(infusion > 1.005, infusion*results$summary$statistics[1], infusion))
      
   #   imm_vec = imm  * dat$immunity_scalar_fitted[t]
      
      
      dat = dat %>%
        mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, 1.24, 1))
      
      imm_vec = imm  * dat$immunity_scalar_fitted[t]

      
      # Seasonal function 
      seasonal = beta*( hum*humidity_at_t + temp*temp_at_t  ) 
      # .00001* humidity_at_t*temp_at_t)
      
      dS = -seasonal*S*(I + I2)/N - VE*vax_rate_at_t*S + imm_vax*V
      dE =  seasonal*S*(I+ I2)/N - lp*E
      dI =  lp*E - gamma*I
      dR = I*gamma -  R*imm_vec  
      dV = VE*vax_rate_at_t*S  - imm_vax*V   
      
      dS2 = -seasonal*S2*(I + I2)/N  + imm_vec*R  + imm_vec*R2 + imm_vax*V2 - VE*vax_rate_at_t*S2
      dE2 =  seasonal*S2*(I + I2)/N - lp*E2
      dI2 =  lp*E2 - gamma*I2
      dR2 = I2*gamma -  R2*imm_vec
      dV2 = VE*vax_rate_at_t*S2  - imm_vax*V2  
      
      # Return the rates of change
      list(c(dS, dE, dI, dR, dV, dS2, dE2, dI2, dR2, dV2))
    })
  }
  
  # Solve the ODE system
  out_check <- ode(y = state, times = time, func = seirs_model_hybrid_state, parms = param)
  
  # hop
  out_df_check_state <- as.data.frame(out_check) %>%
    mutate(predicted_hosp = I* fit3+ I2* fit4) %>%
    mutate(t = time) %>%
    # filter(t > 11) %>%
    left_join(dat, by = "t")
  print(head(out_df_check_state))
  
  rmse_val <- sqrt(mean((out_df_check_state$predicted_hosp - out_df_check_state$hosp)^2, na.rm=TRUE))
  
  params_df <- rbind(params_df, data.frame(
    sample = j,
    sus_frac = fit1,
    beta = fit2,
    h1 = fit3,
    h2 = fit4,
    RMSE = rmse_val
  ))
  
  plot_top = ggplot(data = out_df_check_state) +
    # geom_line(aes(x = date, y = predicted_hosp), col = "aquamarine3", lwd =2.5) +
    geom_line(aes(x = week/52 + year, y = predicted_hosp), col = "aquamarine3", lwd =2) +
    geom_line(aes(x = week/52 + year, y = hosp), col = "black", lwd = 2) +
    ggtitle(unique(out_df_check_state$Location)) +
    theme_bw() +
    ylab("Hospitalizations") + xlab('Date')
  print(plot_top)
  
  plot_top10[[k]] = plot_top
  
  
}


plot_all10 <- wrap_plots(plot_top10, ncol=3)
plot_all10

params_df

##################################################


tn = plot_list[["TN"]]
ga = plot_list[["GA"]]
plot_grid(ga, tn )

ny = plot_list[["NY"]]
ct = plot_list[["CT"]]
plot_grid(ny, ct )


print(plot_list[["MN"]])
print(plot_list[["OH"]])
print(plot_list[["MI"]])

library(zoo)


print(plot_list[["CT"]])
print(plot_list[["MN"]])
print(plot_list[["OH"]])



