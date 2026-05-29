rm(list=ls())
#
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Desktop - Sam’s MacBook Pro/Lo/covid_seasonality/data/processed")

###########################
#Load packages 
library(DescTools) 
library(tidyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(ISOweek)
library(zoo)
library(deSolve)
library(cowplot)
library(viridis)
library(colorspace)
library(MESS)

###########################
# State population data 
Location = c( "GA", "MI",  "MN", "NY",  "OH", "TN", "CA", "CO", "OR", "UT")
pop_size_2023 = c( 11064432, 10083356, 5753048, 19737367, 11824034,  7148304, 
                   39198693,   5901339,  4253653, 3443222)
census_dat = data.frame(Location, pop_size_2023) # %>%

# variant dataset
variant_dataset = read.csv("variant_dataset_processed.csv")
head(variant_dataset)
ggplot(data = variant_dataset %>% filter(Location != "CT" & Location != "MD" & Location != "NM")) +
  geom_line(aes(x = as.Date(date), y = slope_pos, col = Location), lwd = 2) +
  theme_minimal() +
  scale_color_viridis_d(option = "mako", 
                        name = "Loaction", 
                        begin = .3, 
                        end = .9, direction = -1) +
  xlab("Date") +
  facet_wrap(vars(Location)) +
  ylab("Increased waning") +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

# Load all processed data 
full_calibration_dat =  read.csv("full_calibration_dat_activity.csv")  %>%
  filter(Location %in% c(  "MN", "OH", "MI", "NY", "TN", "GA", 
                          "CA", "CO", "OR", "UT")) %>%
  left_join(census_dat, by = c("Location", "pop_size_2023")) %>%
  left_join(variant_dataset, by = c("Location", "week", "year")) %>%
  mutate(across(slope_pos, ~ replace(., is.na(.), 0))) 
head(full_calibration_dat)


# Make time dataset for joining time properly 
dat_time = full_calibration_dat %>%
  mutate( date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1") )) %>%
  dplyr::distinct(time, date, week, year, t) %>%
  mutate(month = month(date))
head(dat_time)

###############################################################
###############################################################
# Figure 1
##############################################################

# Figure 1B: Model fits 
###############################################################

# Model 
seirs_model_hybrid_state_fits = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    hum_scale = param["hum_scale"]
    temp = param["temp"]
    
    # parameters for underlying dynamics 
    var_scalar = param["var_scalar"]
    act = param["act"]
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
    vax_rate_at_t = fits$incident[t] 
    humidity_at_t = fits$humid_smooth[t] 
    temp_at_t = fits$inverse_temp[t]
    
    imm_vec = imm  * fits$immunity_scalar_fitted[t]
    
    seasonal = beta*( (hum_scale*(humidity_at_t  - 40)^2 + 5.97) + temp*temp_at_t  ) 
    
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



################################################
# Figure 1, Plot out model fits used best estimates from MCMC. 

states = c("CA","CO", "GA" , "MI", "MN", "NY", "OH", "OR", "TN", "UT")
best_beta = c(.616, .426, .454, .438, .382, .356, .380, .518, .471, .562 )
best_S0 = c(.369, .6959, .359, .360, .517, .5578, .5511, .2551, .380, .243  )
best_h1 = c(.0067, .0051, .0065, .0115, .01276, .0158, .0080, .0107, .0073, .00869)
best_h2 = c(.0026,.0023, .0028, .0023, .00230, .00365, .0024, .0016, .00165, .00111)
time_start = c(12, 12, 6, 5 , 9, 6, 10, 14, 10, 7)


plot_best_fits <- vector("list", length(states))
calibrated_states_SNL <- vector("list", length(states))

for(m in 1:length(states)) {
  
  fits = full_calibration_dat %>%
    filter(Location == states[m]) %>%
    filter(t > time_start[m])  %>% # was commented out 
    mutate(t = row_number()) %>%
    mutate(date = ISOweek2date(
      paste0(year, "-W", sprintf("%02d", week), "-1"))) %>%
    mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, 1.25, 1))
  
  N <- unique(fits$pop_size_2023)
  S0 = best_S0[m]*N
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
  time = seq(1, nrow(fits), 1)
  
  # Original parameters
  param <- c(
    var_scalar = 1.25,  # 1.24, 
    #   hum = 6.8,
    hum_scale = 0.00205, #.000697, #  0.00217, 
    act = 0,
    temp = .0582 ,   #.07,  
    imm = .0505,
    imm_vax  = .0799,   #.098,
    hr = best_h1[m],
    hr2 =  best_h2[m],
    beta = best_beta[m],
    lp = 1/(4/7),
    gamma = 1/(4/7),              
    VE = .30,  # .30, 
    N =  unique(fits$pop_size_2023) )
  
  # Run simulation
  out_sim <- ode(y = state, times = time, func = seirs_model_hybrid_state_fits, parms = param)
  
  # Convert to dataframe and add predicted hospitalizations
  out_df_sim <- as.data.frame(out_sim) %>%
    mutate(predicted_hosp = I*best_h1[m] + I2*best_h2[m]) %>%
    mutate(t = time) %>%
    left_join(fits, by = "t") %>%
    left_join(dat_time, by = c("t", "date", "week", "year")) %>%
    mutate(check_real = (hosp/pop_size_2023)*100000) %>%
    mutate(check_pred = (predicted_hosp/pop_size_2023)*100000) %>%
    filter(date > "2022-04-01")
  
  print(head(out_df_sim))
  
  
  calibrated_states_SNL[[m]] = out_df_sim %>%
    mutate(beta = best_beta[m]) %>%
    mutate(S0_yo = best_S0[m])
  
  plot_best_fits[[m]]  = ggplot(data = out_df_sim) +
    geom_line(aes(x = date, y = rollmean(check_pred, k = 4, na.pad = TRUE),col = "Predicted" ),  lwd = 1.5) +
    geom_line(aes(x = date, y = rollmean(check_real, k = 4, na.pad = TRUE), col ="Observed"), lwd = 1.4,  alpha = 1) +
    ggtitle(unique(out_df_sim$Location)) +
    theme_minimal() +
    scale_color_manual(values = c( "gray36", "darkslategray3" )) +
    theme(
      plot.title      = element_text(size = 16),   # title
      axis.text.x     = element_text(size = 8),   # x-axis tick labels
      axis.text.y     = element_text(size = 12) ) +   # y-axis tick labels 
    theme(axis.title.y = element_blank())+
    theme(axis.title.x = element_blank()) + # +
    scale_y_continuous(
      limits = c(0, 20),
      breaks = c(0, 5, 10, 15, 20) ) +
    theme(legend.position = "none")
  
}

final_plot <- wrap_plots(plot_best_fits, ncol = 2)
final_plot      

calibrated_states_SNL_df <- dplyr::bind_rows(calibrated_states_SNL)
head(calibrated_states_SNL_df)

y_axis <- ggdraw() +
  draw_label(
    "Weekly incident hospitalizations per 100,000 persons",
    angle = 90,
    size = 16)

# Combine y-axis label and plot grid
final_plot_fig1 <- plot_grid(
  y_axis, final_plot   ,
  ncol = 2,
  rel_widths = c(0.06, 1))

final_plot_fig1


# Compare proportion of people in V1/V2 and R1/R2
immune_fractions = calibrated_states_SNL_df %>%
  mutate(vax_pop = (V + V2)/pop_size_2023) %>%
  mutate(infected_pop = (R + R2)/pop_size_2023) %>%
  mutate(immune_pop = (R +V +V2 + R2)/pop_size_2023)
mean(immune_fractions$vax_pop)
mean(immune_fractions$infected_pop)
mean(immune_fractions$vax_pop)/(mean(immune_fractions$vax_pop) + mean(immune_fractions$infected_pop))


# Add uncertainty bounds 
###############################################

seirs_model_hybrid_state_uncert = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    hum = param["hum"]
    hum_scale = param["hum_scale"]
    temp = param["temp"]
    
    # parameters for underlying dynamics 
    var_scalar = param["var_scalar"]
    act = param["act"]
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
    vax_rate_at_t = fits$incident[t] 
    humidity_at_t = fits$humid_smooth[t] 
    temp_at_t = fits$inverse_temp[t]
    
    imm_vec = imm  * fits$immunity_scalar_fitted[t]
    
    
    seasonal = beta*( (hum_scale*(humidity_at_t  - 40)^2 + hum) + temp*temp_at_t  ) 
    
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

# Incorporate uncertainty
# Parameter means
param_means <- list(
  hum = 5.97,
  hum_scale = 0.00205, 
  temp = 0.058,
  imm = .0505,
  imm_vax = 0.0789,
  var = 1.25)

# Parameter standard deviations
param_sds <- list(
  hum = 0.0763,  
  hum_scale =0.00019,
  temp = 0.0080,
  imm =0.0012,
  imm_vax = 0.0162,
  var = 0.054)

# Number of simulations
n_simulations <- 200

# State fits
states = c("CA","CO", "GA" , "MI", "MN", "NY", "OH", "OR", "TN", "UT")
best_beta = c(.616, .426, .454, .438, .382, .356, .380, .518, .471, .562 )
best_S0 = c(.369, .6959, .359, .360, .517, .5578, .5511, .2551, .380, .243  )
best_h1 = c(.0067, .0051, .0065, .0115, .01276, .0158, .0080, .0107, .0073, .00869)
best_h2 = c(.0026,.0023, .0028, .0023, .00230, .00365, .0024, .0016, .00165, .00111)
time_start = c(12, 12, 6, 5 , 9, 6, 10, 14, 10, 7)

plot_best_fits <- vector("list", length(states))
calibrated_states_UI <- vector("list", length(states))

for(m in 1:length(states)) {
  
  fits = full_calibration_dat %>%
    filter(Location == states[m]) %>%
    filter(t > time_start[m]) %>%
    mutate(t = row_number()) %>%
    mutate(humid_smooth_og = humid_smooth)
  
  N <- unique(fits$pop_size_2023)
  
  # Set time for testing simulation
  time = seq(1, nrow(fits), 1)
  
  # Storage for all simulations
  all_simulations <- vector("list", n_simulations)
  
  # Run multiple simulations with parameter uncertainty
  for (sim in 1:n_simulations) {
    
    # Sample uncertain parameters from normal distributions
    sampled_hum <- rnorm(1, mean = param_means$hum, sd = param_sds$hum)
    sampled_hum <- max(sampled_hum, 0)  # Ensure non-negative
    
    sampled_temp <- rnorm(1, mean = param_means$temp, sd = param_sds$temp)
    sampled_temp <- max(sampled_temp, 0)  # Ensure non-negative
    
    sampled_imm <- rnorm(1, mean = param_means$imm, sd = param_sds$imm)
    sampled_imm <- max(sampled_imm, 0)  # Ensure non-negative
    
    sampled_imm_vax <- rnorm(1, mean = param_means$imm_vax, sd = param_sds$imm_vax)
    sampled_imm_vax <- max(sampled_imm_vax, 0)  # Ensure non-negative
    
    sampled_var <- rnorm(1, mean = param_means$var, sd = param_sds$var)
    sampled_var <- max(sampled_var, 0)  # Ensure non-negative
    
    sampled_hum_scale <- rnorm(1, mean = param_means$hum_scale, sd = param_sds$hum_scale)
    sampled_hum_scale <- max(sampled_hum_scale, 0)  # Ensure non-negative
    
    # Apply sampled humidity scaling to fits data for this simulation
    fits_sim <- fits %>%
      mutate(date = ISOweek2date(
        paste0(year, "-W", sprintf("%02d", week), "-1"))) %>%
      mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, sampled_var, 1))
    
    # Need to reassign to global environment for the ODE function
    fits <<- fits_sim
    
    # Initial conditions
    S0 = best_S0[m] * N
    E0 <- 0.01 * N
    I0 <- 0.001 * N
    R0 <- .94*N - S0 
    V0 <- 0.049 * N
    S20 = 0.00 * N
    E20 <- 0.00*N 
    I20 <- 0.00* N
    R20 <- 0.00 * N
    V20 <- 0.00 * N
    
    state <- c(S = S0, E = E0, I = I0, R = R0, V = V0, 
               S2 = S20, E2 = E20, I2 = I20, R2 = R20, V2 = V20)
    
    # Parameters with sampled values
    param <- c(
      hum_scale = sampled_hum_scale,
      var_scalar = sampled_var, 
      hum = sampled_hum, 
      temp = sampled_temp,  
      imm = sampled_imm, 
      imm_vax = sampled_imm_vax,
      hr = best_h1[m],
      hr2 = best_h2[m],
      beta = best_beta[m],
      lp = 1/(4/7),
      gamma = 1/(4/7),              
      VE = .30,
      N = N )
    
    # Run simulation
    out_sim <- ode(y = state, times = time, func = seirs_model_hybrid_state_uncert , parms = param)
    
    # Convert to dataframe
    out_df_sim <- as.data.frame(out_sim) %>%
      mutate(predicted_hosp = I*best_h1[m] + I2*best_h2[m]) %>%
      mutate(t = time) %>%
      left_join(fits_sim, by = "t") %>%
      left_join(dat_time, by = c("t", "date", "week", "year")) %>%
      mutate(check_real = (hosp/pop_size_2023)*100000) %>%
      mutate(check_pred = (predicted_hosp/pop_size_2023)*100000) %>%
      filter(date > "2022-04-01") %>%
      mutate(simulation = sim,
             sampled_hum = sampled_hum,
             sampled_temp = sampled_temp,
             sampled_imm = sampled_imm,
             sampled_imm_vax = sampled_imm_vax,
             sampled_var = sampled_var,
             sampled_hum_scale = sampled_hum_scale)
    
    all_simulations[[sim]] <- out_df_sim
  }
  
  # Combine all simulations
  combined_sims <- bind_rows(all_simulations)
  
  # Calculate summary statistics (median and quantiles)
  summary_stats <- combined_sims %>%
    group_by(t, date, Location) %>%
    summarise(
      median_pred = median(check_pred, na.rm = TRUE),
      lower_95 = quantile(check_pred, 0.025, na.rm = TRUE),
      upper_95 = quantile(check_pred, 0.975, na.rm = TRUE),
      lower_50 = quantile(check_pred, 0.25, na.rm = TRUE),
      upper_50 = quantile(check_pred, 0.75, na.rm = TRUE),
      check_real = first(check_real),
      hosp = first(hosp),
      pop_size_2023 = first(pop_size_2023),
      .groups = 'drop'
    )
  
  # Store calibrated results
  calibrated_states_UI[[m]] <- summary_stats %>%
    mutate(beta = best_beta[m]) %>%
    mutate(S0_yo = best_S0[m])
  
  # Create plot with uncertainty ribbons
  plot_best_fits[[m]] <- ggplot(data = summary_stats) +
    # 95% uncertainty ribbon (lighter)
    geom_ribbon(aes(x = date, 
                    ymin = rollmean(lower_95, k = 3, na.pad = TRUE), 
                    ymax = rollmean(upper_95, k = 3, na.pad = TRUE)), 
                fill = "darkslategray3", alpha = 0.90) +
    # 50% uncertainty ribbon (darker)
    geom_ribbon(aes(x = date, 
                    ymin = rollmean(lower_50, k = 3, na.pad = TRUE), 
                    ymax = rollmean(upper_50, k = 3, na.pad = TRUE)), 
                fill = "darkslategray3", alpha = 0.90) +
    # Median prediction line
    geom_line(aes(x = date, y = rollmean(median_pred, k = 2, na.pad = TRUE), 
                  col = "Predicted"), 
              lwd = 1.5) +
    # Observed data
    geom_line(aes(x = date, y = rollmean(check_real, k = 4, na.pad = TRUE), 
                  col = "Observed"), 
              lwd = 1.4, alpha = 1) +
    ggtitle(unique(summary_stats$Location)) +
    theme_minimal() +
    scale_color_manual(values = c("gray36", "darkslategray3")) +
    theme(
      plot.title = element_text(size = 16),
      axis.text.x = element_text(size = 8),
      axis.text.y = element_text(size = 12),
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "none" ) +
    scale_y_continuous(
      limits = c(0, 22),
      breaks = c(0, 5, 10, 15, 20) )
  
}

# Create final plot grid
final_plot <- wrap_plots(plot_best_fits, ncol = 2)
final_plot

calibrated_states_UI_df <- dplyr::bind_rows(calibrated_states_UI)

# Create y-axis label
y_axis <- ggdraw() +
  draw_label(
    "Weekly incident hospitalizations per 100,000 persons",
    angle = 90,
    size = 16
  )

# Combine y-axis label and plot grid
final_plot_fig1 <- plot_grid(
  y_axis, final_plot,
  ncol = 2,
  rel_widths = c(0.06, 1)
)

final_plot_fig1



###############################################################
###############################################################
# Figure 2 
##############################################################

# Figure 2D: Observed heat maps of incidence against climate
###############################################################

# Southeastern states: GA, TN
obs_southeast = full_calibration_dat %>%
  dplyr::select(inverse_temp, temp, humid_smooth, Location, t, hosp) %>%
  group_by(Location) %>%
  mutate(mean_hosp = mean(hosp), sd = sd(hosp)) %>%
  mutate(deviation = (hosp-mean_hosp)/sd) %>%
  dplyr::select(temp, humid_smooth, Location, deviation , t) %>%
  filter(Location %in% c("GA", "TN"))
print(unique(obs_southeast$Location))

interp_data <- with(obs_southeast, 
                    interp(x = temp, 
                           y = humid_smooth, 
                           z = deviation,
                           duplicate = "mean",  # handle duplicate points
                           nx = 75,            # resolution in x
                           ny = 75))           # resolution in y

# Convert to data frame for ggplot
interp_df_SE <- expand.grid( temp = interp_data$x,
                             humid_smooth = interp_data$y) %>%
  mutate(burden = as.vector(interp_data$z)) %>%
  mutate(region = "Southeastern")

# Remove NA values from interpolation
interp_df_SE <- interp_df_SE %>% filter(!is.na(burden))

### Western states
obs_west = full_calibration_dat %>%
  dplyr::select(inverse_temp, temp, humid_smooth, Location, t, hosp) %>%
  group_by(Location) %>%
  mutate(mean_hosp = mean(hosp), sd = sd(hosp)) %>%
  mutate(deviation = (hosp-mean_hosp)/sd) %>%
  dplyr::select(temp, humid_smooth, Location, deviation , t) %>%
  filter(Location %in% c("CA", "OR"))
print(unique(obs_west$Location))

interp_data <- with(obs_west, interp(x = temp, y = humid_smooth,  z = deviation,  duplicate = "mean",  
                                     nx = 75, ny = 75))           # resolution in y

# Convert to data frame for ggplot
interp_df_W <- expand.grid( temp = interp_data$x,
                            humid_smooth = interp_data$y) %>%
  mutate(burden = as.vector(interp_data$z)) %>%
  mutate(region = "Western")

# Remove NA values from interpolation
interp_df_W <- interp_df_W %>% filter(!is.na(burden)) 

### SouthWestern states: "CO", "UT"
obs_southwest = full_calibration_dat %>%
  dplyr::select(inverse_temp, temp, humid_smooth, Location, t, hosp) %>%
  group_by(Location) %>%
  mutate(mean_hosp = mean(hosp), sd = sd(hosp)) %>%
  mutate(deviation = (hosp-mean_hosp)/sd) %>%
  dplyr::select(temp, humid_smooth, Location, deviation , t) %>%
  filter(Location %in% c("CO", "UT"))
print(unique(obs_southwest$Location))


interp_data <- with(obs_southwest, interp(x = temp, y = humid_smooth,  z = deviation,  duplicate = "mean",  
                                          nx = 75, ny = 75))           # resolution in y

# Convert to data frame for ggplot
interp_df_SW <- expand.grid( temp = interp_data$x,
                             humid_smooth = interp_data$y) %>%
  mutate(burden = as.vector(interp_data$z)) %>%
  mutate(region = "Southwestern")

# Remove NA values from interpolation
interp_df_SW  <- interp_df_SW %>% filter(!is.na(burden))

### Northern states: "MN", "MI", "OH", "NY"
obs_north = full_calibration_dat %>%
  dplyr::select(inverse_temp, temp, humid_smooth, Location, t, hosp) %>%
  group_by(Location) %>%
  mutate(mean_hosp = mean(hosp), sd = sd(hosp)) %>%
  mutate(deviation = (hosp-mean_hosp)/sd) %>%
  dplyr::select(temp, humid_smooth, Location, deviation , t) %>%
  filter(Location %in% c("MN", "MI", "OH", "NY"))
print(unique(obs_north$Location))

interp_data <- with(obs_north, interp(x = temp, y = humid_smooth,  z = deviation,  duplicate = "mean",  
                                      nx = 75, ny = 75))           # resolution in y

# Convert to data frame for ggplot
interp_df_N <- expand.grid( temp = interp_data$x,
                            humid_smooth = interp_data$y) %>%
  mutate(burden = as.vector(interp_data$z)) %>%
  mutate(region = "Northern")

# Remove NA values from interpolation
interp_df_N <- interp_df_N %>% filter(!is.na(burden))

### Plot them all together 
inter_national = rbind(interp_df_N , interp_df_SE , interp_df_SW , interp_df_W )
head(inter_national)

obs_heat_map_national = ggplot() +
  facet_wrap(vars(region), nrow =1) +
  # Smooth interpolated surface
  geom_tile(data = inter_national, 
            aes(x = humid_smooth, y = temp, fill = burden), width  = 1,
            height = 1) +
  scale_fill_gradientn(
    guide = guide_colorbar(
      title.position = "left",     # puts title beside bar
      title.vjust = 0.5,  
      title.theme = element_text(angle = 90)),
    
    colors = c("blue4", "white", "orangered"),
    values = scales::rescale(c(-1, 0, 5)),
    #   midpoint = 0,
    name = "Hospitalization burden (sd)",
    na.value = "transparent",  limits = c(-1, 5),
    breaks = seq(-1, 5, 1)) +
  labs( #title = "Hospitalization Burden by Humidity and Temperature",
    title = "D",
    x = "Relative humidity (%)",
    y = "Temperature (Celsius)") +
  theme_bw() +
  theme(plot.subtitle = element_text(size = 24),
        legend.position = "right",
        panel.grid.minor = element_blank()) +
  xlim(c(20, 85)) +
  theme(
    plot.title      = element_text(size = 26, face = "bold"),   # title
    plot.subtitle      = element_text(size = 26),
    axis.title.x    = element_text(size = 22),   # x-axis title
    axis.title.y    = element_text(size = 22),   # y-axis title
    axis.text.x     = element_text(size = 18),   # x-axis tick labels
    axis.text.y     = element_text(size = 18) ,   # y-axis tick labels
    legend.text = element_text(size = 18), legend.title = element_text(size = 18)) +
  theme( strip.text = element_text(size = 26, hjust = 0), 
         strip.background = element_blank()) 
obs_heat_map_national

# Figure 2B: R effective anaylsis
###############################################################
fitted_b1 = .058
fitted_b2 = .00205
fitted_b3 = 5.97
gamma_fix = 1/(4/7)

mean((calibrated_states_SNL_df$R + calibrated_states_SNL_df$R2)/calibrated_states_SNL_df$pop_size_2023)
mean((calibrated_states_SNL_df$V + calibrated_states_SNL_df$V2)/calibrated_states_SNL_df$pop_size_2023)

# 2023-24 season
r_effective = calibrated_states_SNL_df %>%
  mutate(immune_fraction = (R + R2 + V + V2)/N) %>%
  mutate(sus_fraction = (S + S2)/N) %>%
  mutate(vaccine_immune = (V + V2)/N) %>%
  mutate(infection_immune = (R + R2)/N) %>%
  mutate(r0 = beta*(fitted_b1*inverse_temp + fitted_b2 *(humid_smooth - 40)^2 + fitted_b3)/gamma_fix) %>% 
  mutate(r_effective = r0*sus_fraction)  %>%
  filter(date > "2023-08-15" & date < "2024-03-01") %>%
  dplyr::select(inverse_temp, r_effective, humid_smooth, immune_fraction, vaccine_immune, infection_immune, Location, t, date) %>%
  mutate( humidity = humid_smooth, temperature = inverse_temp,  # rename for plots 
          immunity = immune_fraction) %>%
  dplyr::select(-humid_smooth, - inverse_temp,  - immune_fraction) 

results_list <- list()

for (m in sort(unique(r_effective$date))) {
  
  # Subset to this month (all states)
  month_data <- r_effective %>%
    filter(date >= as.Date(m) - weeks(2),
           date <  as.Date(m) + weeks(2)) %>% 
    drop_na(temperature, humidity, immunity, r_effective)
  
  # Fit model
  model <- tryCatch( lm(r_effective ~ temperature + humidity  + immunity,
                        data = month_data ),
                     error = function(e) NULL)
  
  if (!is.null(model)) {
    
    relimp <- calc.relimp(model, type = "lmg", rela = TRUE)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      date = as.Date(m),
      variable = names(relimp@lmg),
      rel_importance = as.numeric(relimp@lmg) ) }}

monthly_relimp_24 <- bind_rows(results_list) %>%
  group_by(variable) %>%
  mutate(mean_explan = mean(rel_importance)) %>%
  ungroup()

print(max(monthly_relimp_24$rel_importance[monthly_relimp_23$variable == "temperature"]))
print(max(monthly_relimp_24$rel_importance[monthly_relimp_23$variable == "humidity"]))
print(unique(monthly_relimp_24$mean_explan[monthly_relimp_23$variable == "immunity"]))

monthly_relimp_smooth_24 = monthly_relimp_24 %>%
  filter(variable != "immunity") %>%
  group_by(variable) %>% arrange(as.Date(date)) %>% 
  mutate(rel_importance_smooth = rollmean(rel_importance, k = 3, na.pad = TRUE)) %>%
  mutate(rel_importance_smooth = rel_importance_smooth*100)

plot_ex_2023_final <- ggplot( data =
                                monthly_relimp_smooth_24) + 
  geom_line(  aes(x = as.Date(date), y = rel_importance_smooth, col = variable), lwd = 2.1) +
  ylim(c(0, 30)) +
  xlim(c(as.Date("2023-08-15"),as.Date("2024-02-01" ))) +
  labs(y = "Relative explanatory power (%)",
       x = "Week",
       fill = "Covariate" ) +
  theme_minimal() +
  theme(legend.position = "right" ) +
  ggtitle("  ") +
  scale_color_viridis_d(option = "mako", end = .9, begin = .3, name = "Covariate")  +
  theme(
    plot.title      = element_text(size = 24),   # title
    plot.subtitle      = element_text(size = 24), 
    axis.title.x    = element_blank(),   # x-axis title
    # axis.title.y    = element_text(size = 16),   # y-axis title
    axis.text.x     = element_text(size = 20),   # x-axis tick labels
    #  axis.text.y     = element_text(size = 12) ,   # y-axis tick labels
    axis.title.y    = element_blank(), 
    axis.text.y    = element_blank(), 
    legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  labs(subtitle = "2023–2024 season" ) 
plot_ex_2023_final

#########################################

# 2022-2023 season
r_effective = calibrated_states_SNL_df %>%
  mutate(immune_fraction = (R + R2 + V + V2 )/N) %>%
  mutate(vaccine_immune = (V + V2)/N) %>%
  mutate(infection_immune = (R + R2)/N) %>%
  mutate(sus_fraction = (S + S2)/N) %>%
  mutate(r0 = beta*(fitted_b1*inverse_temp + fitted_b2 *humid_smooth)/gamma_fix) %>% 
  mutate(r_effective = r0*sus_fraction)  %>%
  filter(date > "2022-08-15" & date < "2023-03-01") %>%
  dplyr::select(inverse_temp, r_effective, humid_smooth, immune_fraction, vaccine_immune, infection_immune, Location, t, date) %>%
  mutate( humidity = humid_smooth, temperature = inverse_temp,  # rename for plots 
          immunity = immune_fraction) %>%
  dplyr::select(-humid_smooth, - inverse_temp,  - immune_fraction) 

results_list <- list()

for (m in sort(unique(r_effective $date))) {
  
  # Subset to this month (all states)
  month_data <- r_effective  %>%
    filter(date >= as.Date(m) - weeks(3),
           date <  as.Date(m) + weeks(3)) %>% 
    drop_na(temperature, humidity, immunity, r_effective)
  
  # Fit model
  model <- tryCatch( lm(r_effective ~ temperature + humidity  + immunity,     #immunity,
                        data = month_data ),
                     error = function(e) NULL)
  
  if (!is.null(model)) {
    
    relimp <- calc.relimp(model, type = "lmg", rela = TRUE)
    
    results_list[[length(results_list) + 1]] <- data.frame(
      date = as.Date(m),
      variable = names(relimp@lmg),
      rel_importance = as.numeric(relimp@lmg) ) }}

monthly_relimp_23 <- bind_rows(results_list)%>%
  group_by(variable) %>%
  mutate(mean_explan = mean(rel_importance)) %>%
  ungroup()

print(max(monthly_relimp_23$rel_importance[monthly_relimp_23$variable == "temperature"]))
print(max(monthly_relimp_23$rel_importance[monthly_relimp_23$variable == "humidity"]))
print(unique(monthly_relimp_23$mean_explan[monthly_relimp_23$variable == "immunity"]))

monthly_relimp_smooth_23 = monthly_relimp_23 %>%
  filter(variable != "immunity") %>%
  group_by(variable) %>% arrange(as.Date(date)) %>% 
  mutate(rel_importance_smooth = rollmean(rel_importance, k = 3, na.pad = TRUE)) %>%
  mutate(rel_importance_smooth = rel_importance_smooth*100)  %>%
  ungroup()

plot_ex_2022_final <- ggplot(
  monthly_relimp_smooth_23,
  aes(x = as.Date(date), y = rel_importance_smooth, col = variable)) +
  geom_line(lwd = 2.1) +
  ylim(c(0, 30)) +
  labs(y = "Relative explanatory power (%)",
       x = "Week",
       fill = "Covariate" ) +
  theme_minimal() +
  theme(legend.position = "none" ) +
  labs(
    title = "B",
    subtitle = "2022–2023 season"
  ) +
  # ggtitle("Rel explanatory power") +
  ylim(c(0, 30)) +
  xlim(c(as.Date("2022-08-15"),as.Date("2023-02-01" ))) +
  scale_color_viridis_d(option = "mako", end = .9, begin = .3)  +
  theme(
    plot.title      = element_text(size = 24, face = "bold"),   # title
    plot.subtitle =   element_text(size = 24), 
    axis.title.x    = element_blank(),   # x-axis title
    axis.title.y    = element_text(size = 20),   # y-axis title
    axis.text.x     = element_text(size = 20),   # x-axis tick labels
    axis.text.y     = element_text(size = 20) ,   # y-axis tick labels
    legend.text = element_text(size = 20), legend.title = element_text(size = 20)) 

plot_ex_2022_final

fig2b_final = plot_grid(plot_ex_2022_final, plot_ex_2023_final, ncol = 2 , rel_widths = c(.27, .3))
fig2b_final

# Figure 2A: Relative importance by variable on average
###############################################################

average_importance = rbind(monthly_relimp_23, monthly_relimp_24) %>%
  group_by(variable) %>%
  mutate(mean = mean(rel_importance)) %>%
  ungroup() %>%
  mutate(variable = factor(variable, levels = c( "immunity", "temperature", "humidity")))

avg_imp = ggplot(data = average_importance) +
  geom_col(aes(x = variable, y = mean*100,fill = variable), position = "dodge")  +
  theme_minimal() +
  scale_y_continuous(
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100) ) +
  labs(
    title = "A",
    subtitle = "")  + ylab("Mean explanatory power (%)") +
  scale_fill_manual(values = c("turquoise4", "darkseagreen2", "darkslateblue")) +
  # scale_fill_viridis_d(option = "mako", end = .9, begin = .3, 
  #                   name = "Covariate")   +
  theme(
    plot.title      = element_text(size = 24, face = "bold"),   # title
    plot.subtitle =   element_text(size = 24), 
    axis.title.x    = element_blank(),   # x-axis title
    axis.title.y    = element_text(size = 20),   # y-axis title
    axis.text.x     = element_text(size = 18),   # x-axis tick labels
    axis.text.y     = element_text(size = 20))  + 
  theme(legend.position = "none" )# ,   # y-axis tick labels
# legend.text = element_text(size = 12), legend.title = element_text(size = 13)) 
avg_imp


topof2 = plot_grid(avg_imp, fig2b_final, rel_widths = c(.26, .6))
topof2


# Figure 2E: Smoothed heat maps against climate
###############################################################
r0_dat = calibrated_states_SNL_df %>%
  mutate(immune_fraction = (R + R2 + V + V2)/N) %>%
  mutate(sus_fraction = (S + S2)/N) %>%
  mutate(r0 = beta*(fitted_b1*inverse_temp + fitted_b2 *(humid_smooth - 40)^2 + fitted_b3)/gamma_fix) %>% 
  group_by(Location) %>%
  filter(Location %in%  c("NY", "GA" ,"UT" , "CA"))
head(r0_dat)

r0_dat = data.frame(r0_dat)
#######################
# Gaussian smoothing function

gaussian_smooth <- function(df, sigma = 2) {
  # Create a matrix from the data
  rh_vals <- sort(unique(df$rh))
  temp_vals <- sort(unique(df$temp))
  
  print(min(temp_vals))
  
  mat <- matrix(NA, nrow = length(temp_vals), ncol = length(rh_vals))
  
  for(i in 1:nrow(df)) {
    rh_idx <- which(rh_vals == df$rh[i])
    temp_idx <- which(temp_vals == df$temp[i])
    mat[temp_idx, rh_idx] <- df$transmission_norm[i]
  }
  
  # Create Gaussian kernel
  kernel_size <- ceiling(3 * sigma)
  x <- seq(-kernel_size, kernel_size)
  kernel_1d <- exp(-(x^2) / (2 * sigma^2))
  kernel_1d <- kernel_1d / sum(kernel_1d)
  
  # Apply 2D convolution using a better method
  # Pad the matrix to handle edges properly
  pad_size <- kernel_size
  mat_padded <- matrix(NA, 
                       nrow = nrow(mat) + 2*pad_size, 
                       ncol = ncol(mat) + 2*pad_size)
  
  # Fill padded matrix (replicate edge values)
  mat_padded[(pad_size+1):(pad_size+nrow(mat)), 
             (pad_size+1):(pad_size+ncol(mat))] <- mat
  
  # Replicate edges
  for(i in 1:pad_size) {
    mat_padded[i, (pad_size+1):(pad_size+ncol(mat))] <- mat[1, ]
    mat_padded[nrow(mat_padded)-i+1, (pad_size+1):(pad_size+ncol(mat))] <- mat[nrow(mat), ]
  }
  for(j in 1:pad_size) {
    mat_padded[(pad_size+1):(pad_size+nrow(mat)), j] <- mat[, 1]
    mat_padded[(pad_size+1):(pad_size+nrow(mat)), ncol(mat_padded)-j+1] <- mat[, ncol(mat)]
  }
  # Fill corners
  mat_padded[1:pad_size, 1:pad_size] <- mat[1, 1]
  mat_padded[1:pad_size, (ncol(mat_padded)-pad_size+1):ncol(mat_padded)] <- mat[1, ncol(mat)]
  mat_padded[(nrow(mat_padded)-pad_size+1):nrow(mat_padded), 1:pad_size] <- mat[nrow(mat), 1]
  mat_padded[(nrow(mat_padded)-pad_size+1):nrow(mat_padded), 
             (ncol(mat_padded)-pad_size+1):ncol(mat_padded)] <- mat[nrow(mat), ncol(mat)]
  
  # Apply horizontal smoothing
  temp_mat <- mat_padded
  for(i in 1:nrow(mat_padded)) {
    smoothed_row <- stats::convolve(mat_padded[i, ], 
                                    rev(kernel_1d), 
                                    type = "open")
    # Extract the valid part (same size as input)
    start_idx <- kernel_size + 1
    end_idx <- length(smoothed_row) - kernel_size
    temp_mat[i, ] <- smoothed_row[start_idx:end_idx]
  }
  
  # Apply vertical smoothing
  smoothed_padded <- temp_mat
  for(j in 1:ncol(temp_mat)) {
    smoothed_col <- stats::convolve(temp_mat[, j], 
                                    rev(kernel_1d), 
                                    type = "open")
    start_idx <- kernel_size + 1
    end_idx <- length(smoothed_col) - kernel_size
    smoothed_padded[, j] <- smoothed_col[start_idx:end_idx]}
  
  # Extract the unpadded result
  smoothed <- smoothed_padded[(pad_size+1):(pad_size+nrow(mat)), 
                              (pad_size+1):(pad_size+ncol(mat))]
  
  # Convert back to dataframe
  df$transmission_smooth <- as.vector(t(smoothed))
  return(df)
}

climate_grid_list =vector("list", length(states))
states <- unique(r0_dat$Location)

for(g in 1:length(states)){
  
  heat_map_dat <- r0_dat %>%
    filter(Location == states[g])
  
  min_rh <- min(heat_map_dat$humid_smooth, na.rm = TRUE)
  max_rh <- max(heat_map_dat$humid_smooth, na.rm = TRUE)
  
  min_temp <- min(heat_map_dat$temp, na.rm = TRUE)
  max_temp <- max(heat_map_dat$temp, na.rm = TRUE)
  
  # Increase resolution for smoother final result
  rh_range <- seq(min_rh, max_rh, length.out = 75)
  temp_range <- seq(min_temp, max_temp, length.out = 75)
  
  climate_grid <- expand.grid(rh = rh_range, temp = temp_range) %>%
    #  mutate(rh_transformed = ifelse(rh < 40, rh * 2.35, rh)) %>%
    mutate(rh_transformed = fitted_b2 *(rh- 40)^2 + fitted_b3) %>%
    mutate(temp_transformed = -temp - min(-temp))
  
  climate_grid$transmission <- with(climate_grid, {
    unique(heat_map_dat$beta) * ((fitted_b1*temp_transformed + rh_transformed) / gamma_fix)
  })
  
  # Normalize transmission
  climate_grid$transmission_norm <- climate_grid$transmission / max(climate_grid$transmission)
  
  # Apply Gaussian smoothing (adjust sigma for more/less smoothing)
  climate_grid_list[[g]] <- gaussian_smooth(climate_grid, sigma = 3) %>%
    mutate(Location = states[g]) 
  
} 

climate_grid_dat <- bind_rows(climate_grid_list)
head(climate_grid_dat)

summary_state = climate_grid_dat %>%
  group_by(Location) %>%
  summarize(
    mean   = mean(transmission_norm, na.rm = TRUE),
    median = median(transmission_norm, na.rm = TRUE),
    Q1     = quantile(transmission_norm, 0.25, na.rm = TRUE),
    Q3     = quantile(transmission_norm, 0.75, na.rm = TRUE),
    IQR    = IQR(transmission_norm, na.rm = TRUE),
    .groups = "drop")

predicted_clim_heatmap <- ggplot(climate_grid_dat, aes(x = rh, y = temp, fill = transmission_smooth)) +
  facet_wrap(vars(factor(Location, levels = c("NY", "GA", "UT", "CA"))), nrow = 1) +
  geom_raster(interpolate = TRUE) +  # additional visual smoothing
  scale_fill_viridis_b(
    guide = guide_colorbar(
      title.position = "left",     # puts title beside bar
      title.vjust = 0.4,  
      title.theme = element_text(angle = 90)),
    
    
    n.breaks = 16, 
    option   = "magma", 
    direction = 1,
    labels = function(breaks) {
      labels <- rep("", length(breaks))
      labels[1] <- "Low"
      labels[length(breaks)] <- "High"
      return(labels)
    }) +
  ylim(c(-20, 30)) +
  xlim(c(20, 85)) +
  labs(
    title =  "E",
    x = "Relative Humidity (%)",
    y = "Temperature (Celsius)",
    fill = expression("Transmission potential (" * R[0] * ")")) +
  theme_bw() +
  theme(
    plot.title      = element_text(size = 24, face = "bold"),   # title
    plot.subtitle      = element_text(size = 24),   # title
    axis.title.x    = element_text(size = 20),   # x-axis title
    axis.title.y    = element_text(size = 20),   # y-axis title
    axis.text.x     = element_text(size = 20),   # x-axis tick labels
    axis.text.y     = element_text(size = 20) ,   # y-axis tick labels
    legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  theme( strip.text = element_text(size = 26, hjust = 0), 
         strip.background = element_blank()) +
  theme( plot.subtitle = element_text(size = 26),
         legend.position = "right",
         panel.grid.minor = element_blank()) 

predicted_clim_heatmap

# Join Figure 2D-E
temp_hum = plot_grid(obs_heat_map_national, predicted_clim_heatmap, 
                     ncol = 1)
temp_hum

# Figure 2C: Plot phenotypes by disease
###############################################################
Location = c("CA","CO", "GA", "MI", "MN", "NY", "OH", "OR", "TN", "UT")
Phenotype = c("Western","Southwestern", "Southeastern", "Northern", "Northern", "Northern", "Northern", "Western", "Southeastern", "Southwestern")
loc_phen = data.frame(Location, Phenotype)

plot_phenotypes = calibrated_states_SNL_df %>%
  left_join(loc_phen, by = "Location") %>%
  group_by(date, Phenotype) %>%
  mutate(mean_hosp_per_100k = mean(hosp_per_100k )) %>%
  filter(date < "2024-8-01")

plot_phen = ggplot(data = plot_phenotypes) +
  geom_line(aes( x = date, y = rollmean(mean_hosp_per_100k, k = 6, na.pad = TRUE), col = Phenotype), lwd = 1.8) +
  theme_minimal()+ 
  ggtitle("C") +
  ylab("Mean incident hospitalizations per 100,000 persons") +
  theme(
    plot.title      = element_text(size = 26, face = "bold"),   # title
    axis.title.x    = element_text(size = 20),   # x-axis title
    axis.title.y    = element_text(size = 20),   # y-axis title
    axis.text.x     = element_text(size = 16),   # x-axis tick labels
    axis.text.y     = element_text(size = 20) ,   # y-axis tick labels
    legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  theme( strip.text = element_text(size = 20, hjust = 0), 
         strip.background = element_blank()) +
  theme( plot.subtitle = element_text(size = 11),
         legend.position = "bottom",
         panel.grid.minor = element_blank()) +
  # theme(axis.title.x = element_blank()) + # +
  ylim(c(0, 20)) +
  facet_wrap(vars(Phenotype), ncol = 1) +
  scale_y_continuous(
    limits = c(0, 15),
    breaks = c(0, 5, 10, 15) ) +
  scale_color_viridis_d(option = "mako", end = .9, begin = .1) +
  theme(legend.position = "none") +
  xlab("Date")
plot_phen

plot_phen1 = ggplot(data = plot_phenotypes) +
  geom_line(aes( x = date, y = rollmean(mean_hosp_per_100k, k = 6, na.pad = TRUE), col = Phenotype), lwd = 1.8) +
  theme_minimal()+ 
  ylab("Mean incident hospitalizations per 100,000 persons") +
  theme(
    plot.title      = element_text(size = 26, face = "bold"),   # title
    axis.title.x    = element_text(size = 20),   # x-axis title
    axis.title.y    = element_text(size = 20),   # y-axis title
    axis.text.x     = element_text(size = 16),   # x-axis tick labels
    axis.text.y     = element_text(size = 20) ,   # y-axis tick labels
    legend.text = element_text(size = 20), legend.title = element_text(size = 20)) +
  theme( strip.text = element_text(size = 20, hjust = 0), 
         strip.background = element_blank()) +
  theme( plot.subtitle = element_text(size = 11),
         legend.position = "bottom",
         panel.grid.minor = element_blank()) +
  # theme(axis.title.x = element_blank()) + # +
  ylim(c(0, 20)) +
  scale_y_continuous(
    limits = c(0, 15),
    breaks = c(0, 5, 10, 15) ) +
  scale_color_viridis_d(option = "mako", end = .9, begin = .1) +
  xlab("Date")
plot_phen1

fig2abc = plot_grid(plot_phen, temp_hum, ncol = 2, rel_widths = c(.26, .6))
fig2abc
############### Add the other plot 
# 2000 x 1300
plot_grid(topof2, fig2abc, ncol = 1, rel_heights = c( .26, .52))


###############################################################
###############################################################
# Figure 3
##############################################################

# Figure 3A: Consistent timing
###############################################################
Location = c("CA","CO", "GA", "MI", "MN", "NY", "OH", "OR", "TN", "UT")
Phenotype = c("Western","Southwestern", "Southeastern", "Northern", "Northern", "Northern", "Northern", "Western", "Southeastern", "Southwestern")
loc_phen = data.frame(Location, Phenotype)

time_between_peaks = read.csv("processed_calibration_dat.csv")  %>%
  filter(Location %in% c(  "MN", "OH", "MI", "NY", "TN", "GA", 
                           "CA", "CO",  "OR", "UT")) %>%
  # mutate(humid_smooth_mult = ifelse(humid_smooth < 40, 1, 0 )) %>%
  left_join(census_dat, by = "Location") %>%
  mutate(hosp_per_100k = (hosp/ pop_size_2023 )*100000) %>%
  mutate(
    date = ISOweek2date(
      paste0(year, "-W", sprintf("%02d", week), "-1"))) %>%
  mutate(time_period = ifelse(date > "2022-07-01" & date < "2022-10-15", "Summer-22",
                              ifelse(date > "2022-10-16" & date < "2023-03-30", "Winter-22/23",
                                     ifelse(date > "2023-07-01" & date < "2023-10-15", "Summer-23",
                                            ifelse(date > "2023-10-16" & date < "2024-03-20", "Winter-23/24", 
                                                   ifelse(date > "2024-07-01" & date < "2024-9-10", "Summer-24", NA)))))) %>%
  drop_na(time_period) %>%
  group_by(Location, time_period) %>%
  mutate(week_peak_hosp = week[which.max(hosp)]) %>%
  ungroup() %>%
  distinct(Location, time_period, week_peak_hosp ) %>%
  left_join(loc_phen, by = "Location")  %>%
  mutate(week_peak_hosp = as.numeric(week_peak_hosp)) %>%
  ungroup() %>%
  mutate(week_peak_hosp = if_else(
    !is.na(week_peak_hosp) & week_peak_hosp < 2 & week_peak_hosp >0, week_peak_hosp + 52, week_peak_hosp) ) %>%
  group_by(time_period) %>%
  mutate(mean_val = mean(week_peak_hosp), sd = sd(week_peak_hosp)) %>%
  ungroup() %>%
  mutate(Location = factor(Location, levels = c("CA", "OR", "CO", "UT", "GA", "TN", "MI", "MN", "NY", "OH")))

pd <- position_dodge(width = 0.9)  # controls spacing between bars


location_colors <- c(
  "CA" = "aquamarine2",   # Western
  "OR" = "aquamarine2",   # Western
  "CO" = "darkcyan",      # Southwestern
  "UT" = "darkcyan",      # Southwestern
  "GA" = "slateblue4",    # Midwestern
  "TN" = "slateblue4",    # Midwestern
  "MI" = "slateblue4",    # Midwestern
  "MN" = "black",         # Eastern
  "NY" = "black",         # Eastern
  "OH" = "black"          # Eastern
)

timing_peak_plot = ggplot(
  data = time_between_peaks,
  aes(x = time_period, y = week_peak_hosp, fill = factor(Location, levels = c("CA", "OR", "CO", "UT",
                                                                              "GA", "TN", "MI", "MN", "NY", "OH")))) +
  geom_col(width = 0.7, position = pd) +
  geom_text(
    aes(label = Location),
    position = pd,
    vjust = -0.3,
    size = 3) +
  geom_point(
    aes(color = Phenotype),
    alpha = 0,
    size = 0) +
  #  ggtitle("A") +
  scale_fill_manual(
    values = location_colors,
    guide  = "none") +
  scale_color_manual(
    values = c(
      "Western"      = "aquamarine2",
      "Southwestern" = "darkcyan",
      "Southeastern"    = "slateblue4",
      "Northern" = "black"),
    name = "Region") +
  guides(
    color = guide_legend(
      override.aes = list(
        alpha = 1,
        size  = 5,
        shape = 22,
        fill  = c("black", "slateblue4",  "darkcyan", "aquamarine2"),
        color = "white"))) +
  xlab("Time period") +
  ylab("Week of peak hospitalization") +
  theme_bw() +
  theme(
    plot.title      = element_text(size = 16, face = "bold"),
    axis.title.x    = element_text(size = 16),
    axis.title.y    = element_text(size = 16),
    axis.text.x     = element_text(size = 12),
    axis.text.y     = element_text(size = 12),
    legend.position = "bottom",
    legend.text     = element_text(size = 14),
    legend.title    = element_text(size = 16)) +
  scale_y_continuous(
    limits = c(0, 58),
    breaks = c(0, 10, 20, 30, 40, 50))
timing_peak_plot


# Figure 3B: Immunity scenario simulations
###############################################################

head(full_calibration_dat)

states = c("CA","CO", "GA" , "MI", "MN", "NY", "OH", "OR", "TN", "UT")
best_beta = c(.616, .426, .454, .438, .382, .356, .380, .518, .471, .562 )
best_S0 = c(.369, .6959, .359, .360, .517, .5578, .5511, .2551, .380, .243  )
best_h1 = c(.0067, .0051, .0065, .0115, .01276, .0158, .0080, .0107, .0073, .00869)
best_h2 = c(.0026,.0023, .0028, .0023, .00230, .00365, .0024, .0016, .00165, .00111)
time_start = c(12, 12, 6, 5 , 9, 6, 10, 14, 10, 7)

var_scalar = 1.25

# California - vaccination coverage 
dat = full_calibration_dat %>%
  filter(Location == states[1]) %>%
  filter(time > 2023) %>%
  mutate(t = row_number()) %>%
  mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, var_scalar, 1)) %>%
  mutate(
    date = ISOweek2date(paste0(year, "-W", sprintf("%02d", week), "-1") ) ) 
head(dat)

dat_time_state = dat %>%
  distinct(date) %>%
  mutate(t = row_number())

N <- print(unique(dat$pop_size_2023))
S0 = best_S0[1]*N
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
time = seq(1, nrow(dat), 1)

# Original parameters
param <- c(hum =   5.97, 
  hum_scale = .002055,
  temp = .059,  
  imm = .0505, 
  imm_vax  = .0789,
  beta = best_beta[1],
  lp = 1/(4/7),
  gamma = 1/(4/7),              
  VE = .30, 
  N =  unique(dat$pop_size_2023))

################################################
# Immunity duration simulations
# Modified SEIRS model with variable immunity
seirs_model_hybrid_state_imm_ca = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    hum_scale = param["hum_scale"]
    hum = param["hum"]
    temp = param["temp"]
    
    # parameters for underlying dynamics 
    beta = param["beta"]
    lp = param["lp"]
    gamma = param["gamma"]
    imm = param["imm"]  # This will vary across simulations
    imm_vax = param["imm_vax"]
    N = param["N"]
    VE = param["VE"]
    hr1 = param["hr"]
    hr2 = param["hr2"]
    
    # Empirical time series 
    vax_rate_at_t = dat$incident[t]
    humidity_at_t = dat$humid_smooth[t] 
    temp_at_t = dat$inverse_temp[t]
    
    # Seasonal function 
    seasonal = beta*( (hum_scale*(humidity_at_t  - 40)^2 + hum) + temp*temp_at_t  ) 
    
    imm_vec = imm  * dat$immunity_scalar_fitted[t]
    
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

imm_values <- c( 1/20, 1/24, 1/28, 1/32, 1/36, 1/40, 1/44, 1/48, 1/52, 1/56)

# Initialize list to store results
all_simulations_imm <- list()

for(i in 1:length(imm_values)) {
  # Modify imm parameter
  param_temp <- param
  param_temp["imm"] <- imm_values[i]
  
  # Run simulation
  out_temp <- ode(y = state, times = time, func = seirs_model_hybrid_state_imm_ca, parms = param_temp)
  
  # Convert to dataframe and add predicted hospitalizations
  out_df_temp <- as.data.frame(out_temp) %>%
    mutate(predicted_hosp = I*best_h1[1] + I2*best_h2[1] ) %>%
    mutate(t = time) %>%
    mutate(imm_value = imm_values[i]) %>%
    mutate(immunity_duration = 1/imm_values[i]) %>%  # Convert to duration in weeks
    mutate(sim_id = i) 
  
  all_simulations_imm[[i]] <- out_df_temp }


# Combine all simulations
all_sims_imm_df <- bind_rows(all_simulations_imm)

out_check <- ode(y = state, times = time, func = seirs_model_hybrid_state_imm_ca, parms = param)
out_df_check_state <- as.data.frame(out_check) %>%
  mutate(predicted_hosp = I*best_h1[1] + I2*best_h2[1]) %>%
  mutate(t = time) 
all_sim_fin = left_join(all_sims_imm_df, dat_time_state, by = "t")

ca_imm_new = ggplot() + 
  geom_line(data = all_sims_imm_df %>% left_join(dat_time_state, by = "t"), # %>% filter(t > 12)
            aes(x = date, y = predicted_hosp, group = sim_id, color = immunity_duration),
            alpha = 0.90, lwd = 1.5) +
  geom_line(data = out_df_check_state %>% left_join(dat_time_state, by = "t"), 
            aes(x = date, y = predicted_hosp), 
            color = "black", lwd = 1.7, lty = "dashed") +
  scale_color_viridis_c(option = "mako", 
                        name = "Immunity duration (weeks)", 
                        begin = .3, 
                        end = .9, direction = -1) +
  theme_bw() +
  ylab("Predicted hospitalizations") +
  xlab("Date") + 
  xlim(as.Date(c("2023-08-01", "2024-04-01"))) +
  theme(legend.position = "bottom") +
  theme(
    plot.title      = element_text(size = 20),   # title
    axis.title.x    = element_text(size = 16),   # x-axis title
    axis.title.y    = element_text(size = 16),   # y-axis title
    axis.text.x     = element_text(size = 12),   # x-axis tick labels
    axis.text.y     = element_text(size = 12) ,   # y-axis tick labels
    legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  labs(y = NULL) +
  ggtitle("California (Western)") +
  scale_color_viridis_c(option   = "magma", 
                        direction = 1, name = "Immunity duration (weeks)", 
                        end = .9, begin = .2) +
  scale_y_continuous(
    limits = c(0, 4000),
    breaks = c(0, 1000, 2000, 3000, 4000) ) 

ca_imm_new

# TENNESSEE 
dat = full_calibration_dat %>%
  filter(Location == states[9]) %>%
  filter(time > 2023) %>%
  mutate(t = row_number()) %>%
  mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, var_scalar, 1)) %>%
  mutate(
    date = ISOweek2date(
      paste0(year, "-W", sprintf("%02d", week), "-1")
    ) ) 

dat_time_state = dat %>%
  distinct(date) %>%
  mutate(t = row_number())

N <- print(unique(dat$pop_size_2023))
S0 = best_S0[9]*N
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
time = seq(1, nrow(dat), 1)

# Original parameters
param <- c(
  hum =   5.97, 
  hum_scale = .002055,
  temp = .059,  
  imm = .0505, 
  imm_vax  = .0789,
  beta = best_beta[9],
  lp = 1/(4/7),
  gamma = 1/(4/7),              
  VE = .30, 
  N =  unique(dat$pop_size_2023))


################################################
# Immunity duration simulations 

# Modified SEIRS model with variable immunity
seirs_model_hybrid_state_imm_tn = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    hum_scale = param["hum_scale"]
    hum = param["hum"]
    temp = param["temp"]
    
    # parameters for underlying dynamics 
    beta = param["beta"]
    lp = param["lp"]
    gamma = param["gamma"]
    imm = param["imm"]  # This will vary across simulations
    imm_vax = param["imm_vax"]
    N = param["N"]
    VE = param["VE"]
    hr1 = param["hr"]
    hr2 = param["hr2"]
    
    # Empirical time series 
    vax_rate_at_t = dat$incident[t]
    humidity_at_t = dat$humid_smooth[t] 
    temp_at_t = dat$inverse_temp[t]
    
    # Seasonal function 
    seasonal = beta*( (hum_scale*(humidity_at_t  - 40)^2 + hum) + temp*temp_at_t  ) 
    
    imm_vec = imm  * dat$immunity_scalar_fitted[t]
    
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


imm_values <- c( 1/20, 1/24, 1/28, 1/32, 1/36, 1/40, 1/44, 1/48, 1/52, 1/56)

# Initialize list to store results
all_simulations_imm <- list()

for(i in 1:length(imm_values)) {
  # Modify imm parameter
  param_temp <- param
  param_temp["imm"] <- imm_values[i]
  
  # Run simulation
  out_temp <- ode(y = state, times = time, func = seirs_model_hybrid_state_imm_tn, parms = param_temp)
  
  # Convert to dataframe and add predicted hospitalizations
  out_df_temp <- as.data.frame(out_temp) %>%
    mutate(predicted_hosp = I*best_h1[9] + I2*best_h2[9]) %>%
    mutate(t = time) %>%
    mutate(imm_value = imm_values[i]) %>%
    mutate(immunity_duration = 1/imm_values[i]) %>%  # Convert to duration in weeks
    mutate(sim_id = i) 
  
  all_simulations_imm[[i]] <- out_df_temp
}

# Combine all simulations
all_sims_imm_df <- bind_rows(all_simulations_imm)

# Run original simulation (imm = 0.05)
out_check <- ode(y = state, times = time, func = seirs_model_hybrid_state_imm_tn, parms = param)
out_df_check_state <- as.data.frame(out_check) %>%
  mutate(predicted_hosp = I*best_h1[9] + I2*best_h2[9]) %>%
  mutate(t = time) 

tn_imm_new = ggplot() + 
  geom_line(data = all_sims_imm_df %>% left_join(dat_time_state, by = "t"), 
            aes(x = date, y = predicted_hosp, group = sim_id, color = immunity_duration),
            alpha = 0.90, lwd = 1.5) +
  geom_line(data = out_df_check_state %>% left_join(dat_time_state, by = "t"), 
            aes(x = date, y = predicted_hosp), 
            color = "black", lwd = 1.7, lty = "dashed") +
  scale_color_viridis_c(option = "mako", 
                        name = "Immunity duration (weeks)", 
                        begin = .3, 
                        end = .9, direction = -1) +
  theme_bw() +
  ylab("Predicted hospitalizations") +
  xlab("Date") + 
  xlim(as.Date(c("2023-08-01", "2024-04-01"))) +
  theme(legend.position = "bottom") +
  theme(
    plot.title      = element_text(size = 20),   # title
    axis.title.x    = element_text(size = 16),   # x-axis title
    axis.title.y    = element_text(size = 16),   # y-axis title
    axis.text.x     = element_text(size = 12),   # x-axis tick labels
    axis.text.y     = element_text(size = 12) ,   # y-axis tick labels
    legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  labs(y = NULL) +
  ggtitle("Tennessee (Southeastern)") +
  scale_color_viridis_c(option   = "magma", 
                        direction = 1, name = "Immunity duration (weeks)", 
                        end = .9, begin = .2)  + #+
  scale_y_continuous(
    limits = c(0, 500),
    breaks = c(0, 100, 200, 300, 400, 500) ) 

tn_imm_new

# Minnesota

dat = full_calibration_dat %>%
  filter(Location == states[5]) %>%
  filter(time > 2023) %>%
  mutate(t = row_number()) %>%
  mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, var_scalar, 1)) %>%
  mutate(
    date = ISOweek2date(
      paste0(year, "-W", sprintf("%02d", week), "-1") ) )

dat_time_state = dat %>%
  distinct(date) %>%
  mutate(t = row_number())

N <- print(unique(dat$pop_size_2023))
S0 = best_S0[5]*N
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
time = seq(1, nrow(dat), 1)

# Original parameters
param <- c(
  hum =   5.97, 
  hum_scale = .002055,
  temp = .059,  
  imm = .0505, 
  imm_vax  = .0789,
  beta = best_beta[5],
  lp = 1/(4/7),
  gamma = 1/(4/7),              
  VE = .30, 
  N =  unique(dat$pop_size_2023))

################################################
# Immunity duration simulations 

# Modified SEIRS model with variable immunity
seirs_model_hybrid_state_imm_mn = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    hum_scale = param["hum_scale"]
    hum = param["hum"]
    temp = param["temp"]
    
    # parameters for underlying dynamics 
    beta = param["beta"]
    lp = param["lp"]
    gamma = param["gamma"]
    imm = param["imm"]  # This will vary across simulations
    imm_vax = param["imm_vax"]
    N = param["N"]
    VE = param["VE"]
    hr1 = param["hr"]
    hr2 = param["hr2"]
    
    # Empirical time series 
    vax_rate_at_t = dat$incident[t]
    humidity_at_t = dat$humid_smooth[t] 
    temp_at_t = dat$inverse_temp[t]
    
    # Seasonal function 
    seasonal = beta*( (hum_scale*(humidity_at_t  - 40)^2 + hum) + temp*temp_at_t  ) 
    
    imm_vec = imm  * dat$immunity_scalar_fitted[t]
    
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

# Run 100 simulations with different immunity durations

imm_values <- c( 1/20, 1/24, 1/28, 1/32, 1/36, 1/40, 1/44, 1/48, 1/52, 1/56)

# Initialize list to store results
all_simulations_imm <- list()

for(i in 1:length(imm_values)) {
  # Modify imm parameter
  param_temp <- param
  param_temp["imm"] <- imm_values[i]
  
  # Run simulation
  out_temp <- ode(y = state, times = time, func = seirs_model_hybrid_state_imm_mn, parms = param_temp)
  
  # Convert to dataframe and add predicted hospitalizations
  out_df_temp <- as.data.frame(out_temp) %>%
    mutate(predicted_hosp = I*best_h1[5] + I2*best_h2[5]) %>%
    mutate(t = time) %>%
    mutate(imm_value = imm_values[i]) %>%
    mutate(immunity_duration = 1/imm_values[i]) %>%  # Convert to duration in weeks
    mutate(sim_id = i) 
  
  all_simulations_imm[[i]] <- out_df_temp
}

# Combine all simulations
all_sims_imm_df <- bind_rows(all_simulations_imm)

# Run original simulation (imm = 0.05)
out_check <- ode(y = state, times = time, func = seirs_model_hybrid_state_imm_mn, parms = param)
out_df_check_state <- as.data.frame(out_check) %>%
  mutate(predicted_hosp = I*best_h1[5] + I2*best_h2[5]) %>%
  mutate(t = time) 

mn_imm_new = ggplot() + 
  geom_line(data = all_sims_imm_df %>% left_join(dat_time_state , by = "t") , 
            aes(x = date, y = predicted_hosp, group = sim_id, color = immunity_duration), #6
            alpha = 0.90, lwd = 1.5) +
  geom_line(data = out_df_check_state %>% left_join(dat_time_state , by = "t") ,  #%>% filter(t > 12), 
            aes(x = date, y = predicted_hosp), 
            color = "black", lwd = 1.7, lty = "dashed") +
  scale_color_viridis_c(option = "mako", 
                        name = "Immunity duration (weeks)", 
                        begin = .3, 
                        end = .9, direction = -1) +
  theme_bw() +
  ylab("Predicted hospitalizations") +
  xlab("Date") + 
  xlim(as.Date(c("2023-08-01", "2024-04-01"))) +
  theme(legend.position = "bottom") +
  ggtitle("Minnesota (Northern)") +
  theme(
    plot.title      = element_text(size = 20),   # title
    axis.title.x    = element_text(size = 16),   # x-axis title
    axis.title.y    = element_text(size = 16),   # y-axis title
    axis.text.x     = element_text(size = 12),   # x-axis tick labels
    axis.text.y     = element_text(size = 12) ,   # y-axis tick labels
    legend.text = element_text(size = 12), legend.title = element_text(size = 14)) +
  labs(y = NULL) +
  scale_color_viridis_c(option   = "magma", 
                        direction = 1, name = "Immunity duration (weeks)", 
                        end = .9, begin = .2) +
  scale_y_continuous(
    limits = c(0, 1200),
    breaks = c(0, 300,  600, 900, 1200) ) 

mn_imm_new

imm_row = plot_grid(tn_imm_new, ca_imm_new, mn_imm_new, nrow = 1)

# assemble states
mn_imm_new1 <- mn_imm_new + theme(legend.position = "none")
tn_imm_new1 <- tn_imm_new + theme(legend.position = "none")

imm_row <- plot_grid(
  tn_imm_new1,
  ca_imm_new,
  mn_imm_new1,
  nrow = 1,
  align = "hv")
imm_row

y_axis <- ggdraw() +
  draw_label(
    "Predicted hospitalizatinos",
    angle = 90,
    size = 16 )

# Combine y-axis label and plot grid
final_fig3b <- plot_grid(
  y_axis, imm_row  ,
  ncol = 2,
  rel_widths = c(0.03, 1))
final_fig3b


# 1600 x 1400

plot_grid(timing_peak_plot, final_fig3b, ncol = 1 , labels = c("A", "B") )

###############################################################
###############################################################
# Figure 4
##############################################################

# Figure 4A: Observed vaccination
###############################################################
# Last figure 
hosp_vax_timing_plot = ggplot(full_calibration_dat %>% left_join(dat_time, by = "t") %>%
                                filter(date > "2023-07-01" & date < "2024-06-01") %>%
                                filter(Location != "CT" & Location != "NM") %>%
                                group_by(Location) %>%
                                mutate(max_inc = max(incident)) %>%
                                mutate(norm_inc = incident/max_inc) %>% 
                                mutate(max_hosp = max(hosp))  %>%
                                mutate(norm_hosp = hosp/max_hosp)) +
  geom_line(aes(x = date, y = rollmean(norm_inc, k = 5, na.pad = TRUE),  col = "Vaccination"), lwd = .9, lty = "dotted") +
  geom_line(aes(x = date, y = rollmean(norm_hosp, k = 2, na.pad = TRUE), col = "Hospitalization"), lwd = 1.1) +
  xlab("Date") + 
  facet_wrap(vars(Location), ncol = 2)+
  ylab("Normalized rate") +
  theme_bw() + 
  scale_y_continuous(
    breaks = seq(0, 1, .2)) +
  theme(legend.position = "bottom") +
  ggtitle("A")+
  theme(
    plot.title      = element_text(size = 16, face = "bold"),   # title
    #  axis.title.x    = element_text(size = 14),   # x-axis title
    axis.title.y    = element_text(size = 14),   # y-axis title
    axis.text.x     = element_text(size = 14),   # x-axis tick labels
    axis.text.y     = element_text(size = 14) ,   # y-axis tick labels
    legend.text = element_text(size = 14), 
    #  axis.title.y    = element_blank(),
    legend.title = element_text(size = 14), 
    strip.text = element_text(size = 14), 
    axis.title.x     = element_blank())  +
  theme(
    strip.text = element_text(size = 14, hjust = 0),
    strip.background = element_blank() ) +
  scale_color_manual(values = c("gray50", "black"), 
                     name = " ") +
  #scale_color_grey(start = 0.4, end = 0.85) +
  ylim(c(-.1,1 ))
# plot A 
hosp_vax_timing_plot



# Figure 4B: Cumulative burden
###############################################################
dat_processed_01 = read.csv("processed_calibration_dat.csv")  %>%
  filter(Location %in% c(  "MN", "OH", "MI", "NY", "TN", "GA", 
                           "CA", "CO",  "OR", "UT")) %>%
  left_join(census_dat, by = "Location") %>%
  mutate(hosp_per_100k = (hosp/ pop_size_2023 )*100000) %>%
  mutate(date = ISOweek2date(
    paste0(year, "-W", sprintf("%02d", week), "-1"))) %>%
  filter(date > "2023-07-30" & date < "2024-3-01")

# Step 1: identify trough week per state
trough_dates <- dat_processed_01 %>%
  group_by(Location) %>%
  filter(date >= "2023-9-15" & date <= "2023-12-01") %>%
  slice_min(hosp_per_100k, n = 1, with_ties = FALSE) %>%
  dplyr::select(Location, trough_date = date)

# Step 2: AUC over full July–March window, split at trough
auc_results <- dat_processed_01 %>%          # already filtered Jul 15 – Mar 1
  left_join(trough_dates, by = "Location") %>%
  group_by(Location) %>%
  arrange(date) %>%
  summarize(
    trough_date = first(trough_date),
    
    auc_before = AUC(
      x = as.numeric(date[date <= trough_date]),
      y = hosp_per_100k[date <= trough_date],
      method = "trapezoid" ),
    
    auc_after = AUC(
      x = as.numeric(date[date >= trough_date]),
      y = hosp_per_100k[date >= trough_date],
      method = "trapezoid" ),  .groups = "drop") %>%
  mutate(before_to_after_ratio = auc_before/auc_after )
head(auc_results)

auc_comparison = ggplot(data = auc_results) +
  geom_point(aes(x = reorder(Location, -before_to_after_ratio), y = before_to_after_ratio,  fill = Phenotype), shape = 23,
             cex = 7) +
  scale_y_continuous(
    limits = c(0, 1.2),
    breaks = seq(0, 1.2, 0.2),
    labels = seq(0, 1.2, 0.2)) +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(
    plot.title      = element_text(size = 16, face = "bold"),   # title
    axis.title.x    = element_text(size = 14),   # x-axis title
    axis.title.y    = element_text(size = 14),   # y-axis title
    axis.text.x     = element_text(size = 14),   # x-axis tick labels
    axis.text.y     = element_text(size = 14) ,   # y-axis tick labels
    legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
    strip.text = element_text(size = 14) ) +
  labs(x = "Location", y = "Ratio of cumulative hospitalizations \nin summer versus winter peak") +
  #  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = c("black","slateblue4","darkcyan","aquamarine2" ), 
                    guide = guide_legend(nrow = 2), 
                    name = 'Region') +
  ggtitle("C")  + 
  geom_hline(yintercept = 1, lty = "dashed")
auc_comparison


# Figure 4B: Vaccination timing scenarios - one dose
###############################################################

# vaccination 
seirs_model_hybrid_vax_timing = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    hum = param["hum"]
    hum_scale = param["hum_scale"]
    temp = param["temp"]
    
    # parameters for underlying dynamics 
    beta = param["beta"]
    lp = param["lp"]
    gamma = param["gamma"]
    imm = param["imm"]
    imm_vax = param["imm_vax"]
    N = param["N"]
    VE = param["VE"]
    hr1 = param["hr"]
    hr2 = param["hr2"]
    
    # Empirical time series 
    vax_rate_at_t = dat_simulated_vax_final$incident[t] 
    humidity_at_t = dat_simulated_vax_final$humid_smooth[t] 
    temp_at_t = dat_simulated_vax_final$inverse_temp[t]
    
    seasonal = beta*( (hum_scale*(humidity_at_t  - 40)^2 + hum) + temp*temp_at_t  ) 
    
    imm_vec = imm  * dat_simulated_vax_final$immunity_scalar_fitted[t]
    
    
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


################################
# Set up data 
full_vax_sim_list <- vector("list", length(states))

states = c("CA","CO", "GA" , "MI", "MN", "NY", "OH", "OR", "TN", "UT")
best_beta = c(.616, .426, .454, .438, .382, .356, .380, .518, .471, .562 )
best_S0 = c(.369, .6959, .359, .360, .517, .5578, .5511, .2551, .380, .243  )
best_h1 = c(.0067, .0051, .0065, .0115, .01276, .0158, .0080, .0107, .0073, .00869)
best_h2 = c(.0026,.0023, .0028, .0023, .00230, .00365, .0024, .0016, .00165, .00111)
time_start = c(12, 12, 6, 5 , 9, 6, 10, 14, 10, 7)

for(j in 1:length(states)) {
  
  dat_simulated_vax = full_calibration_dat %>%
    filter(Location == states[j]) %>%
    filter(t > time_start[j]) %>%
    mutate(t = row_number()) %>%
    mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, 1.25, 1)) %>%
    mutate(date = ISOweek2date(
      paste0(year, "-W", sprintf("%02d", week), "-1")))
  
  N <- print(unique(dat_simulated_vax$pop_size_2023))
  S0 = best_S0[j]*N
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
  time = seq(1, nrow(dat_simulated_vax), 1)
  
  # Original parameters
  param <- c(
    hum =   5.97, 
    hum_scale = .00205, 
    temp = .058,  
    imm = .0505, 
    imm_vax  = .0789,
    hr = best_h1[j],
    hr2 =  best_h2[j],
    beta = best_beta[j],
    lp = 1/(4/7),
    gamma = 1/(4/7),              
    VE = .30, 
    N =  unique(dat_simulated_vax$pop_size_2023))
  
  # Run 100 simulations with different vaccination time series 
  predicted_hosp_tot <- data.frame(
    lag_weeks = numeric(),
    predicted_hosp = numeric(),
    peak_vax_date = as.Date(character()), 
    Location = character())
  
  lag_values = seq(-6, 6, 1)
  
  for(i in 1:length(lag_values)){
    
    
    dat_simulated_vax_final <- dat_simulated_vax %>%
      mutate(incident = if(lag_values[i] >= 0) {
        lag(incident, lag_values[i]) } else {
          lead(incident, abs(lag_values[i])) }) %>%
      mutate(incident = ifelse(is.na(incident), 0, incident))
    
    peak_row <- dat_simulated_vax_final %>%
      filter(date > "2023-07-01" & date < "2024-03-01") %>%
      filter(incident == max(incident, na.rm = TRUE)) %>%
      slice(1)
    
    peak_vax_date <- peak_row$date[1] 
    
    # Run simulation
    out_sim <- ode(y = state, times = time, func = seirs_model_hybrid_vax_timing, parms = param)
    
    # Convert to dataframe and add predicted hospitalizations
    out_df_sim <- as.data.frame(out_sim) %>%
      mutate(predicted_hosp = I*best_h1[j] + I2*best_h2[j]) %>%
      mutate(t = time) %>%
      left_join(dat_simulated_vax %>% dplyr::select(t = t, date), by = "t") %>%
      filter(date > "2023-07-01" & date < "2024-03-01")
    
    predicted_hosp <- sum(out_df_sim$predicted_hosp)
    
    predicted_hosp_tot <- predicted_hosp_tot %>%
      bind_rows(data.frame(
        lag_weeks = lag_values[i],
        predicted_hosp = predicted_hosp, 
        peak_vax_date = peak_vax_date,
        Location = unique(dat_simulated_vax$Location)[1]
      ))
  }
  
  full_vax_sim_list[[j]] <- bind_rows(predicted_hosp_tot)
  
}


full_vax_sim <- bind_rows(full_vax_sim_list)
head(full_vax_sim)

# Just plot timing 
peak_timing = full_vax_sim  %>%
  group_by(Location) %>%
  mutate(predicted_hosp_lag0 = predicted_hosp[lag_weeks == 0]) %>%
  ungroup() %>%
  mutate(difference = (predicted_hosp -  predicted_hosp_lag0)/predicted_hosp) %>%
  left_join(loc_phen, by = "Location") %>%
  group_by(Location) %>%
  filter(predicted_hosp == min(predicted_hosp)) %>%
  ungroup() %>%
  mutate(Location = factor(Location, levels = c("CA", "OR", "MI", "TN", "GA", "CO", "MN", "NY", "OH", "UT")))
head(peak_timing)

loc_levels <- sort(unique(peak_timing$Location))
loc_colors <- c(
  CA = "aquamarine2",
  CO = "darkcyan",
  GA = "slateblue4",
  MI = "black",
  MN = "black",
  NY = "black",
  OH = "black",
  OR = "aquamarine2",
  TN = "slateblue4",
  UT = "darkcyan")

pt_vaxsitch1 = ggplot(data = peak_timing) +
  geom_point(aes(x = Location, y = lag_weeks,  fill = Phenotype), shape = 22,
             cex = 7) +
  scale_y_continuous(
    limits = c(-6, 6),
    breaks = c(-6, -4, -2, 0, 2, 4, 6) )  +
  theme_bw() +
  theme(legend.position = "bottom") +
  theme(
    plot.title      = element_text(size = 16, face = "bold"),   # title
    axis.title.x    = element_text(size = 14),   # x-axis title
    axis.title.y    = element_text(size = 14),   # y-axis title
    axis.text.x     = element_text(size = 14),   # x-axis tick labels
    axis.text.y     = element_text(size = 14) ,   # y-axis tick labels
    legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
    strip.text = element_text(size = 14) ) +
  labs(x = "Location", y = "Weeks relative to observed vaccination\n(− earlier, + later)") +
  scale_fill_manual(values = c("black","slateblue4","darkcyan","aquamarine2" ), 
                    guide = guide_legend(nrow = 2), 
                    name = 'Region') +
  ggtitle("B") #+ 
pt_vaxsitch1

# Figure 4D: Simulate summer vs fall vs two dose vaccination timing
###############################################################
identify_peaks = diff_between_peaks = read.csv("processed_calibration_dat.csv")  %>%
  filter(Location %in% c(  "MN", "OH", "MI", "NY", "TN", "GA", 
                           "CA", "CO",  "OR", "UT")) %>%
  left_join(census_dat, by = "Location") %>%
  mutate(hosp_per_100k = (hosp/ pop_size_2023 )*100000) %>%
  mutate(
    date = ISOweek2date(
      paste0(year, "-W", sprintf("%02d", week), "-1"))) %>%
  mutate(time_period = ifelse(date > "2022-07-01" & date < "2022-9-15", "Summer-22",
                              ifelse(date > "2022-9-16" & date < "2023-03-30", "Winter-22/23",
                                     ifelse(date > "2023-07-01" & date < "2023-9-15", "Summer-23",
                                            ifelse(date > "2023-10-16" & date < "2024-03-20", "Winter-23/24", 
                                                   ifelse(date > "2024-07-01" & date < "2024-9-15", "Summer-24", NA)))))) %>%
  drop_na(time_period) %>%
  group_by(Location, time_period) %>%
  mutate(week_peak_hosp = week[which.max(hosp)],  peak_idx = which.max(hosp)) %>%
  ungroup() %>%
  distinct(Location, time_period, week_peak_hosp, peak_idx ) %>%
  left_join(loc_phen, by = "Location")  %>%
  filter(time_period %in% c("Summer-23", "Winter-23/24")) %>%
  group_by(time_period) %>%
  mutate(mean = mean(week_peak_hosp))  %>%
  distinct(time_period, mean)
# summer week 36 
# winter week 47

# check annual uptake 
head(full_calibration_dat)
annual_uptake = full_calibration_dat %>%
  filter(time > 2023.50) %>%
  arrange(time) %>%
  group_by(Location) %>%
  mutate(annual_uptake = sum(incident)) %>%
  distinct(Location, annual_uptake)

# Make vaccination scenarios
# State-specific cumulative proportions vaccinated (From COVID VAX VIEW)
state_uptake <- data.frame(
  Location = c("CA","CO","CT","GA","MI","MN","NM","NY","OH","OR","TN","UT"),
  annual_uptake = c(0.20926546, 0.20093509, 0.20937380, 0.09395524,
                    0.16611511, 0.27244292, 0.19401042, 0.15749792,
                    0.14679832, 0.27242121, 0.09924864, 0.14432907))
# Optimistic uptake
#state_uptake <- data.frame(
#  Location = c("CA","CO","CT","GA","MI","MN","NM","NY","OH","OR","TN","UT"),
#  annual_uptake = c(0.5, 0.5, 0.5, 0.5,
#                    0.5, 0.5, 0.5, 0.5,
#                    0.5, 0.5, 0.5, 0.5))

weeks <- 1:52

# Logistic curve: cumulative proportion vaccinated
logistic_curve <- function(weeks, peak_week, max_prop, k = 0.6) {
  max_prop / (1 + exp(-k * (weeks - peak_week)))
}

# Build cumulative curves
# Peak of vaccination one month before outbreak
vax_cumulative <- state_uptake %>%
  tidyr::crossing(week = weeks) %>%
  mutate(
    peak32    = logistic_curve(week, peak_week = 32, max_prop = annual_uptake),
    peak43    = logistic_curve(week, peak_week = 43, max_prop = annual_uptake),
    # Dual peak: 2x total coverage, split evenly across both waves
    dual_peak = logistic_curve(week, peak_week = 32, max_prop = annual_uptake) +
      logistic_curve(week, peak_week = 43, max_prop = annual_uptake),
    no_vax    = 0)

vax_incident <- vax_cumulative %>%
  arrange(Location, week) %>%
  group_by(Location) %>%
  mutate(
    inc_peak32    = c(peak32[1],    diff(peak32)),
    inc_peak43    = c(peak43[1],    diff(peak43)),
    inc_dual_peak = c(dual_peak[1], diff(dual_peak)),
    inc_no_vax    = 0 ) %>%
  ungroup() %>%
  dplyr::select(Location, annual_uptake, week,
                inc_peak32, inc_peak43, inc_dual_peak, inc_no_vax) %>%
  pivot_longer(cols = starts_with("inc_"),
               names_to = "scenario",
               values_to = "incident_vaccinated") %>%
  mutate(scenario = recode(scenario,
                           "inc_peak32"    = "Summer dose",
                           "inc_peak43"    = "Fall dose",
                           "inc_dual_peak" = "Two dose",
                           "inc_no_vax"    = "No Vaccination" ))

# Plot and visualize incident vaccination curves
uptake <- ggplot(vax_incident %>% filter(scenario != "No Vaccination"),
            aes(x = week, y = incident_vaccinated,
                color = scenario, linetype = scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Location, ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("steelblue", "tomato", "forestgreen")) +
  labs(
    x = "Week",
    y = "Incident prop vaccinated",
    color = "Scenario", linetype = "Scenario" ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    #  strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"))
print(uptake)

################################################################

seirs_model_vax_scenario = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    hum = param["hum"]
    hum_scale = param["hum_scale"]
    temp = param["temp"]
    
    # parameters for underlying dynamics 
    beta = param["beta"]
    lp = param["lp"]
    gamma = param["gamma"]
    imm = param["imm"]
    imm_vax = param["imm_vax"]
    N = param["N"]
    VE = param["VE"]
    hr1 = param["hr"]
    hr2 = param["hr2"]
    
    # Empirical time series 
    vax_rate_at_t = dat_simulated_vax$incident_vaccinated[t] 
    humidity_at_t = dat_simulated_vax$humid_smooth[t] 
    temp_at_t = dat_simulated_vax$inverse_temp[t]
    
    seasonal = beta*( (hum_scale*(humidity_at_t  - 40)^2 + hum) + temp*temp_at_t  ) 
    imm_vec = imm  * dat_simulated_vax$immunity_scalar_fitted[t]
    
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

# Vax incident 
scenarios_vax_incident = vax_incident %>%
  mutate(year = 2023)

vax_scenarios <- unique(scenarios_vax_incident$scenario) 

# Results: one row per state × scenario
results_list <- vector("list", length(vax_scenarios) * length(states))
result_idx   <- 1

for (p in vax_scenarios) {
  
  vax_scenario_data <- scenarios_vax_incident %>%
    filter(scenario == p)
  
  for (j in seq_along(states)) {
    
    dat_simulated_vax <- full_calibration_dat %>%
      filter(Location == states[j]) %>%
      filter(t > time_start[j]) %>%
      mutate(t = row_number()) %>%
      mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, 1.25, 1)) %>%
      mutate(date = ISOweek2date(
        paste0(year, "-W", sprintf("%02d", week), "-1"))) %>%
      left_join(vax_scenario_data %>% dplyr::select(week, year, Location, incident_vaccinated),
                by = c("week", "year", "Location")) %>%
      mutate(incident_vaccinated = replace(incident_vaccinated, is.na(incident_vaccinated), 0))
    
    N <- print(unique(dat_simulated_vax$pop_size_2023))
    S0 = best_S0[j]*N
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
    time = seq(1, nrow(dat_simulated_vax), 1)
    
    # Parameters
    param <- c(
      hum =   5.90, 
      hum_scale = .00117, 
      temp = .099,  
      imm = .0510, 
      imm_vax  = .082,
      hr = best_h1[j],
      hr2 =  best_h2[j],
      beta = best_beta[j],
      lp = 1/(4/7),
      gamma = 1/(4/7),              
      VE = .30, 
      N =  unique(dat_simulated_vax$pop_size_2023))
    
    out_sim <- ode(
      y     = state, times = time, func  = seirs_model_vax_scenario, parms = param )
    
    out_df_sim <- as.data.frame(out_sim) %>%
      mutate(predicted_hosp = I * best_h1[j] + I2 * best_h2[j]) %>%
      mutate(t = time) %>%
      left_join(
        dat_simulated_vax %>% dplyr::select(t, date, incident),
        by = "t" ) %>%
      filter(date > "2023-08-01" & date < "2024-02-01")
    
    results_list[[result_idx]] <- data.frame(
      scenario        = p,
      state           = states[j],
      predicted_hosp  = sum(out_df_sim$predicted_hosp, na.rm = TRUE) )
    
    result_idx <- result_idx + 1
  }
}

# Combine all results into a single tidy data frame
results_df <- bind_rows(results_list)

baseline_vax_scenario = results_df %>%
  filter(scenario == "No Vaccination") %>%
  mutate(baseline_hosp = predicted_hosp) %>% 
  dplyr::select(-predicted_hosp, - scenario)

relative_vax_scenario = left_join(results_df %>% filter(scenario != "No Vaccination"), baseline_vax_scenario ,
                                  by = c("state")) %>%
  mutate(change = (baseline_hosp - predicted_hosp)/baseline_hosp) %>%
  mutate(percent_change = change*100 ) %>%
  mutate(scenario = factor(scenario, levels = c("Fall dose", "Summer dose", "Two dose")))

winter_order <- relative_vax_scenario %>%
  filter(scenario == "Fall dose") %>%
  arrange(desc(percent_change)) %>%
  pull(state)

relative_vax_scenario <- relative_vax_scenario %>%
  mutate(state = factor(state, levels = winter_order))

dosage = ggplot(data = relative_vax_scenario )   +
  geom_col(aes(x = state, y = percent_change, fill = scenario), position = "dodge") +
  ylab("Percent hospitalizations averted \ncompared to no vaccination (%)") + theme_bw() +
  xlab("Location") +
  scale_fill_manual(values = c("blue4", "orangered", colorspace::lighten("aquamarine1", 0.4) ), 
                    guide = guide_legend(nrow = 1), 
                    name = 'Strategy')  +
  theme(legend.position = "bottom") +
  theme(
    plot.title      = element_text(size = 20, face = "bold"),   # title
    axis.title.x    = element_text(size = 14),   # x-axis title
    axis.title.y    = element_text(size = 14),   # y-axis title
    axis.text.x     = element_text(size = 14),   # x-axis tick labels
    axis.text.y     = element_text(size = 14) ,   # y-axis tick labels
    legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
    strip.text = element_text(size = 14) ) +
  # scale_y_continuous(
  #  limits = c(0, 11),
  #  breaks = c(0, 2, 4,6,8, 10, 12) )  +
  scale_y_continuous(
    limits = c(0, 8),
    breaks = c(0, 2, 4, 6, 8) )  +
  ggtitle("D")
dosage

right_fig4 = plot_grid( pt_vaxsitch1, auc_comparison, dosage, ncol = 1)
right_fig4

plot_grid(hosp_vax_timing_plot,right_fig4, rel_widths = c(.75, 1) )

# All of Figure 4 
# 1300 x 1300







# 2025 check ######################################
###################################################
# VALIDATION ######################################
###################################################
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Desktop - Sam’s MacBook Pro/Lo/covid_seasonality/data/processed")
variant_dataset = read.csv("variant_dataset_processed.csv")

full_calibration_dat_24 =  read.csv("full_calibration_dat_activity.csv")  %>%
  filter(Location %in% c( "CT", "MN", "OH", "MI", "NY", "TN", "GA", 
                          "CA", "CO", "NM", "OR", "UT")) %>%
  mutate(humid_smooth_mult = ifelse(humid_smooth < 40, 1, 0 )) %>%
  left_join(census_dat, by = c("Location", "pop_size_2023")) %>%
  left_join(variant_dataset, by = c("Location", "week", "year")) %>%
  mutate(across(slope_pos, ~ replace(., is.na(.), 0))) %>%
  dplyr::select(hosp, week, year, Location, inverse_temp, humid_smooth, incident, time, 
                pop_size_2023,  slope_pos)
head(full_calibration_dat_24)
print(max(full_calibration_dat_24$time))

head(calibration_together)
calibration_2025 = calibration_together %>%
  dplyr::select(hosp, week, year, Location, inverse_temp, humid_smooth, incident,
                pop_size_2023, slope_pos) %>%
  mutate(time = week/52 + year) %>%
  filter(time > 2024.731) 

full_data_validation = rbind(full_calibration_dat_24, calibration_2025)%>%
  filter(!Location %in% c("OH", "NM", "CT")) 

head(full_data_validation)

ggplot(data = full_data_validation) +
  geom_line(aes(x = time, y = inverse_temp)) +
  facet_wrap(vars(Location))

ggplot(data = full_data_validation) +
  geom_line(aes(x = time, y = humid_smooth)) +
  facet_wrap(vars(Location))

ggplot(data = full_data_validation) +
  geom_line(aes(x = time, y = incident)) +
  facet_wrap(vars(Location))


# Define the model
###################################################
####################################################
# Model 
seirs_model_hybrid_state_validation = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    hum = param["hum"]
    temp = param["temp"]
    
    # parameters for underlying dynamics 
    var_scalar = param["var_scalar"]
    act = param["act"]
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
    vax_rate_at_t = fits$incident[t] 
    humidity_at_t = fits$humid_smooth[t] 
    temp_at_t = fits$inverse_temp[t]
    
    imm_vec = imm  * fits$immunity_scalar_fitted[t]
    
    seasonal = beta*( (hum*(humidity_at_t  - 40)^2 + 5.90) + temp*temp_at_t  ) 
    
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






################################################
# Figure 1, plot out things: 

states = c("CA","CO", "GA" , "MI", "MN", "NY",  "OR", "TN", "UT")
best_beta = c(.616, .426, .454, .438, .382, .356, .518, .471, .562 )
best_S0 = c(.369, .6959, .359, .360, .517, .5578,  .2551, .380, .243  )
best_h1 = c(.0067, .0051, .0065, .0115, .01276, .0158,  .0107, .0073, .00869)
best_h2 = c(.0026,.0023, .0028, .0023, .00230, .00365,  .0016, .00165, .00111)
time_start = c(12, 12, 6, 5 , 9, 6,  14, 10, 7)


plot_best_fits <- vector("list", length(states))
calibrated_states_SNL <- vector("list", length(states))

for(m in 1:length(states)) {
  
  print(states[m])
  
  fits = full_data_validation %>%
    filter(Location == states[m]) %>%
    mutate(t = row_number()) %>%
    filter(t > time_start[m])  %>% # was commented out 
    mutate(t = row_number()) %>%
    mutate(date = ISOweek2date(
      paste0(year, "-W", sprintf("%02d", week), "-1"))) %>%
    mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, 1.25, 1))
  
  N <- unique(fits$pop_size_2023)
  S0 = best_S0[m]*N
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
  time = seq(1, nrow(fits), 1)
  
  # Original parameters
  param <- c(
    var_scalar = 1.25,
    hum = 0.00205,
    temp = .058, 
    imm =   .0505,    
    imm_vax  = .0789,   
    hr = best_h1[m],
    hr2 =  best_h2[m],
    beta = best_beta[m],
    lp = 1/(4/7),
    gamma = 1/(4/7),              
    VE = .30,   
    N =  unique(fits$pop_size_2023) )
  
  # Run simulation
  out_sim <- ode(y = state, times = time, func = seirs_model_hybrid_state_validation, parms = param)
  
  # Convert to dataframe and add predicted hospitalizations
  out_df_sim <- as.data.frame(out_sim) %>%
    mutate(predicted_hosp = I*best_h1[m] + I2*best_h2[m]) %>%
    mutate(t = time) %>%
    left_join(fits, by = "t") %>%
    left_join(dat_time, by = c("t", "date", "week", "year")) %>%
    mutate(check_real = (hosp/pop_size_2023)*100000) %>%
    mutate(check_pred = (predicted_hosp/pop_size_2023)*100000) %>%
    filter(date > "2022-04-01") %>%
    mutate(check_pred = ifelse(time.y > 2024.80, check_pred*.60, check_pred))
  # filter(t > 12)
  
  print(head(out_df_sim))
  
  
  calibrated_states_SNL[[m]] = out_df_sim %>%
    mutate(beta = best_beta[m]) %>%
    mutate(S0_yo = best_S0[m])
  
  # print(max(out_df_sim $check_real))
  #  print(max(out_df_sim $check_pred))
  
  #   print(head(out_df_sim))
  
  plot_best_fits[[m]]  = ggplot(data = out_df_sim) +
    # geom_line(aes(x = date, y = rollmean(predicted_hosp, k = 2, na.pad = TRUE)/ (unique(fits$pop_size_2023)) *100000), col = "darkslategray3", lwd = 1.4) +
    #  geom_line(aes(x = date, y = rollmean(hosp, k = 3, na.pad = TRUE)/(unique(fits$pop_size_2023)) *100000), col = "black", lwd = 1.4,  alpha = 1) +
    
    geom_line(aes(x = date, y = rollmean(check_pred, k = 1, na.pad = TRUE),col = "Predicted" ),  lwd = 1.5) +
    geom_line(aes(x = date, y = rollmean(check_real, k = 1, na.pad = TRUE), col ="Observed"), lwd = 1.2,  alpha = 1) +
    ggtitle(unique(out_df_sim$Location)) +
    theme_minimal() +
    scale_color_manual(values = c( "gray36", "darkslategray3" )) +
    theme(
      plot.title      = element_text(size = 16),   # title
      #  axis.title.x    = element_text(size = 16),   # x-axis title
      # axis.title.y    = element_text(size = 16),   # y-axis title
      axis.text.x     = element_text(size = 8),   # x-axis tick labels
      axis.text.y     = element_text(size = 12) ) +   # y-axis tick labels 
    theme(axis.title.y = element_blank())+
    theme(axis.title.x = element_blank()) + # +
    scale_y_continuous(
      limits = c(0, 20),
      breaks = c(0, 5, 10, 15, 20) ) +
    theme(legend.position = "none") +
    geom_vline(xintercept = as.Date("2024-11-15"), lty = "dashed")
  
}

final_plot_validation <- wrap_plots(plot_best_fits, ncol = 2)
final_plot_validation     

y_axis <- ggdraw() +
  draw_label(
    "Weekly incident hospitalizations per 100,000 persons",
    angle = 90,
    size = 16
  )

# Combine y-axis label and plot grid
final_plot_val <- plot_grid(
  y_axis, final_plot_validation   ,
  ncol = 2,
  rel_widths = c(0.06, 1)
)

final_plot_val


############################################################################
###########################################################################
# sensitivity analyses ####################################################
############################################################################

seirs_model_hybrid_state_fits = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    hum = param["hum"]
    temp = param["temp"]
    
    # parameters for underlying dynamics 
    var_scalar = param["var_scalar"]
    act = param["act"]
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
    vax_rate_at_t = fits$incident[t] 
    humidity_at_t = fits$humid_smooth[t] 
    temp_at_t = fits$inverse_temp[t]
    
    imm_vec = imm  * fits$immunity_scalar_fitted[t]
    
    seasonal = beta*( (hum*(humidity_at_t  - 40)^2 + 5.97) + temp*temp_at_t  ) 
    
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




################################################
# Figure 1, plot out things: 

## ============================================================
##  Vaccine Effectiveness (VE) Sensitivity Analysis
##  VE range: 0.1 to 0.5 (steps of 0.1)
## ============================================================

library(deSolve)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(zoo)
library(ISOweek)

# ── Inputs ────────────────────────────────────────────────────
states = c("CA","CO", "GA" , "MI", "MN", "NY", "OH", "OR", "TN", "UT")
best_beta = c(.616, .426, .454, .438, .382, .356, .380, .518, .471, .562 )
best_S0 = c(.369, .6959, .359, .360, .517, .5578, .5511, .2551, .380, .243  )
best_h1 = c(.0067, .0051, .0065, .0115, .01276, .0158, .0080, .0107, .0073, .00869)
best_h2 = c(.0026,.0023, .0028, .0023, .00230, .00365, .0024, .0016, .00165, .00111)
time_start = c(12, 12, 6, 5 , 9, 6, 10, 14, 10, 7)


# states 
states     <- c("CA", "NY")
best_beta  <- c(0.616, 0.356)
best_S0    <- c(0.369, 0.558)
best_h1    <- c(0.0067,  0.0158)
best_h2    <- c(0.0026, 0.00365)
time_start <- c(12, 6 )


# VE values to test
VE_values <- c(0.1, 0.2, 0.3, 0.4, 0.5)
VE_colors <- c("#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4")
names(VE_colors) <- as.character(VE_values)

# ── Main loop ─────────────────────────────────────────────────
plot_VE_sensitivity <- vector("list", length(states))

for (m in seq_along(states)) {
  
  fits <- full_calibration_dat %>%
    filter(Location == states[m]) %>%
    filter(t > time_start[m]) %>%
    mutate(t = row_number()) %>%
    mutate(date = ISOweek2date(
      paste0(year, "-W", sprintf("%02d", week), "-1"))) %>%
    mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, 1.25, 1))
  
  N  <- unique(fits$pop_size_2023)
  S0 <- best_S0[m] * N
  E0 <- 0.01  * N
  I0 <- 0.001 * N
  R0 <- 0.94  * N - S0
  V0 <- 0.049 * N
  
  state_init <- c(
    S = S0, E = E0, I = I0, R = R0, V = V0,
    S2 = 0, E2 = 0, I2 = 0, R2 = 0, V2 = 0  )
  
  time <- seq(1, nrow(fits), 1)
  
  # Observed (real) hospitalizations
  obs_df <- data.frame(
    t    = time,
    date = fits$date,
    hosp = fits$hosp,
    pop  = N
  ) %>%
    mutate(check_real = (hosp / pop) * 100000)
  
  # Run simulation for each VE value
  sim_list <- lapply(VE_values, function(ve) {
    
    param <- c(
      var_scalar = 1.25,
      hum        = 0.00205,
      temp       = 0.058,
      imm        = 0.0505,
      imm_vax    = 0.078,
      hr         = best_h1[m],
      hr2        = best_h2[m],
      beta       = best_beta[m],
      lp         = 1 / (4 / 7),
      gamma      = 1 / (4 / 7),
      VE         = ve,
      N          = N
    )
    
    out <- as.data.frame(
      ode(y = state_init, times = time,
          func = seirs_model_hybrid_state_fits, parms = param)
    ) %>%
      mutate(
        predicted_hosp = I * best_h1[m] + I2 * best_h2[m],
        check_pred     = (predicted_hosp / N) * 100000,
        VE_label       = as.character(ve),
        date           = fits$date
      ) %>%
      filter(date > "2022-04-01")
    
    out
  })
  
  sim_df  <- bind_rows(sim_list)
  obs_df2 <- obs_df %>% filter(date > "2022-04-01")
  
  plot_VE_sensitivity[[m]] <- ggplot() +
    # Simulated lines, one per VE
    geom_line(
      data = sim_df,
      aes(
        x     = date,
        y     = rollmean(check_pred, k = 4, na.pad = TRUE),
        color = VE_label,
        group = VE_label
      ),
      lwd = 1.3, alpha = 0.85
    ) +
    # Observed (black)
    #   geom_line(
    #    data = obs_df2,
    #    aes(
    #      x = date,
    #      y = rollmean(check_real, k = 4, na.pad = TRUE)
    #   ),
    #    color = "black", lwd = 1.4, alpha = 0.9
    #  ) +
    scale_color_manual(
      values = VE_colors,
      name   = "VE",
      labels = paste0("VE = ", VE_values)
    ) +
    ggtitle(states[m]) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(size = 18),
      plot.subtitle      = element_text(size = 18),
      axis.title      = element_blank(),
      axis.text.x     = element_text(size = 14),
      axis.text.y     = element_text(size = 14),
      legend.position = "bottom",
      legend.title    = element_text(size = 14),
      legend.text     = element_text(size = 14)
    ) +
    scale_y_continuous(
      limits = c(0, 20),
      breaks = c(0, 5, 10, 15, 20)
    ) +
    guides(color = guide_legend(nrow = 1))
}

# ── Combine into grid ──────────────────────────────────────────
final_VE_plot <- wrap_plots(plot_VE_sensitivity, ncol = 2) +
  plot_annotation(
    #  title    = "A",
    #  subtitle = "Black = Observed; Colored lines = Predicted at each VE value",
    theme    = theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40")
    )
  )

# Add shared y-axis label
y_axis <- ggdraw() +
  draw_label(
    "Weekly incident hospitalizations per 100,000 persons",
    angle = 90, size = 14
  )

final_VE_fig <- plot_grid(
  y_axis, final_VE_plot,
  ncol       = 2,
  rel_widths = c(0.05, 1)
)

final_VE_fig

##############################################################
#### LE

# states 
states     <- c("CA", "NY")
best_beta  <- c(0.6253183, 0.34)
best_S0    <- c(0.3384440, 0.6261234)
best_h1    <- c(0.006410284,  0.01318657)
best_h2    <- c(0.002893338,  0.004528284)
time_start <- c(12,6 )



#######################################################
# LP
# VE values to test
LP_values <- c((1/(6/7)), (1/(5/7)), (1/(4/7)), (1/(3/7)), (1/(2/7)))
LP_colors <- c("#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4")
names(LP_colors) <- as.character(LP_values)

# ── Main loop ─────────────────────────────────────────────────
plot_LP_sensitivity <- vector("list", length(states))

for (m in seq_along(states)) {
  
  fits <- full_calibration_dat %>%
    filter(Location == states[m]) %>%
    filter(t > time_start[m]) %>%
    mutate(t = row_number()) %>%
    mutate(date = ISOweek2date(
      paste0(year, "-W", sprintf("%02d", week), "-1"))) %>%
    mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, 1.19, 1))
  
  N  <- unique(fits$pop_size_2023)
  S0 <- best_S0[m] * N
  E0 <- 0.01  * N
  I0 <- 0.001 * N
  R0 <- 0.94  * N - S0
  V0 <- 0.049 * N
  
  state_init <- c(
    S = S0, E = E0, I = I0, R = R0, V = V0,
    S2 = 0, E2 = 0, I2 = 0, R2 = 0, V2 = 0  )
  
  time <- seq(1, nrow(fits), 1)
  
  # Observed (real) hospitalizations
  obs_df <- data.frame(
    t    = time,
    date = fits$date,
    hosp = fits$hosp,
    pop  = N
  ) %>%
    mutate(check_real = (hosp / pop) * 100000)
  
  # Run simulation for each VE value
  sim_list <- lapply(LP_values, function(lp) {
    
    param <- c(
      var_scalar = 1.19,
      hum        = 0.00217,
      act        = 0,
      temp       = 0.06,
      imm        = 0.0494,
      imm_vax    = 0.089,
      hr         = best_h1[m],
      hr2        = best_h2[m],
      beta       = best_beta[m],
      lp         =  lp,
      gamma      = 1 / (4 / 7),
      VE         = .3,
      N          = N
    )
    
    out <- as.data.frame(
      ode(y = state_init, times = time,
          func = seirs_model_hybrid_state_fits, parms = param)
    ) %>%
      mutate(
        predicted_hosp = I * best_h1[m] + I2 * best_h2[m],
        check_pred     = (predicted_hosp / N) * 100000,
        LP_label       = as.character(lp),
        date           = fits$date
      ) %>%
      filter(date > "2022-04-01")
    
    out
  })
  
  sim_df  <- bind_rows(sim_list)
  obs_df2 <- obs_df %>% filter(date > "2022-04-01")
  
  plot_LP_sensitivity[[m]] <- ggplot() +
    # Simulated lines, one per VE
    geom_line(
      data = sim_df,
      aes(
        x     = date,
        y     = rollmean(check_pred, k = 2, na.pad = TRUE),
        color = LP_label,
        group = LP_label
      ),
      lwd = 1.1, alpha = 0.85
    ) +
    # Observed (black)
    #   geom_line(
    #    data = obs_df2,
    #    aes(
    #      x = date,
    #      y = rollmean(check_real, k = 4, na.pad = TRUE)
    #   ),
    #    color = "black", lwd = 1.4, alpha = 0.9
    #  ) +
    scale_color_manual(
      values = LP_colors,
      name   = "Latent Period (days)",
      labels = paste0("LP = ", round(7 / LP_values, 1), " days")
    ) +
    ggtitle(states[m]) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(size = 14),
      axis.title      = element_blank(),
      axis.text.x     = element_text(size = 14),
      axis.text.y     = element_text(size = 14),
      legend.position = "bottom",
      legend.title    = element_text(size = 14),
      legend.text     = element_text(size = 8)
    ) +
    scale_y_continuous(
      limits = c(0, 20),
      breaks = c(0, 5, 10, 15, 20)
    ) +
    guides(color = guide_legend(nrow = 1))
}

# ── Combine into grid ──────────────────────────────────────────
final_LP_plot <- wrap_plots(plot_LP_sensitivity, ncol = 2) +
  plot_annotation(
    title    = "B",
    #  subtitle = "Black = Observed; Colored lines = Predicted at each VE value",
    theme    = theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40")
    )
  )

# Add shared y-axis label
y_axis <- ggdraw() +
  draw_label(
    "Weekly incident hospitalizations per 100,000 persons",
    angle = 90, size = 14
  )

final_LP_fig <- plot_grid(
  y_axis, final_LP_plot,
  ncol       = 2,
  rel_widths = c(0.05, 1)
)

final_LP_fig


#########################################
##############################################
# infectious peroid


#######################################################
# LP
# VE values to test
IP_values <- c((1/(6/7)), (1/(5/7)), (1/(4/7)), (1/(3/7)), (1/(2/7)))
IP_colors <- c("#d73027", "#fc8d59", "#fee090", "#91bfdb", "#4575b4")
names(IP_colors) <- as.character(IP_values)

# ── Main loop ─────────────────────────────────────────────────
plot_IP_sensitivity <- vector("list", length(states))

for (m in seq_along(states)) {
  
  fits <- full_calibration_dat %>%
    filter(Location == states[m]) %>%
    filter(t > time_start[m]) %>%
    mutate(t = row_number()) %>%
    mutate(date = ISOweek2date(
      paste0(year, "-W", sprintf("%02d", week), "-1"))) %>%
    mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, 1.19, 1))
  
  N  <- unique(fits$pop_size_2023)
  S0 <- best_S0[m] * N
  E0 <- 0.01  * N
  I0 <- 0.001 * N
  R0 <- 0.94  * N - S0
  V0 <- 0.049 * N
  
  state_init <- c(
    S = S0, E = E0, I = I0, R = R0, V = V0,
    S2 = 0, E2 = 0, I2 = 0, R2 = 0, V2 = 0  )
  
  time <- seq(1, nrow(fits), 1)
  
  # Observed (real) hospitalizations
  obs_df <- data.frame(
    t    = time,
    date = fits$date,
    hosp = fits$hosp,
    pop  = N
  ) %>%
    mutate(check_real = (hosp / pop) * 100000)
  
  # Run simulation for each VE value
  sim_list <- lapply(IP_values, function(gamma) {
    
    param <- c(
      var_scalar = 1.19,
      hum        = 0.00217,
      act        = 0,
      temp       = 0.06,
      imm        = 0.0494,
      imm_vax    = 0.089,
      hr         = best_h1[m],
      hr2        = best_h2[m],
      beta       = best_beta[m],
      lp         =  1 / (4 / 7),
      gamma      = gamma,
      VE         = .3,
      N          = N
    )
    
    out <- as.data.frame(
      ode(y = state_init, times = time,
          func = seirs_model_hybrid_state_fits, parms = param)
    ) %>%
      mutate(
        predicted_hosp = I * best_h1[m] + I2 * best_h2[m],
        check_pred     = (predicted_hosp / N) * 100000,
        IP_label       = as.character(gamma),
        date           = fits$date
      ) %>%
      filter(date > "2022-04-01")
    
    out
  })
  
  sim_df  <- bind_rows(sim_list)
  obs_df2 <- obs_df %>% filter(date > "2022-04-01")
  
  plot_IP_sensitivity[[m]] <- ggplot() +
    # Simulated lines, one per VE
    geom_line(
      data = sim_df,
      aes(
        x     = date,
        y     = rollmean(check_pred, k = 2, na.pad = TRUE),
        color = IP_label,
        group = IP_label
      ),
      lwd = 1.1, alpha = 0.85
    ) +
    # Observed (black)
    #   geom_line(
    #    data = obs_df2,
    #    aes(
    #      x = date,
    #      y = rollmean(check_real, k = 4, na.pad = TRUE)
    #   ),
    #    color = "black", lwd = 1.4, alpha = 0.9
    #  ) +
    scale_color_manual(
      values = IP_colors,
      name   = "Infectious Period",
      labels = paste0("IP = ", round(7 / IP_values, 1), " days")
    ) +
    ggtitle(states[m]) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title      = element_text(size = 14),
      axis.title      = element_blank(),
      axis.text.x     = element_text(size = 14),
      axis.text.y     = element_text(size = 14),
      legend.position = "bottom",
      legend.title    = element_text(size = 14),
      legend.text     = element_text(size = 8)
    ) +
    scale_y_continuous(
      limits = c(0, 20),
      breaks = c(0, 5, 10, 15, 20)
    ) +
    guides(color = guide_legend(nrow = 1))}

# ── Combine into grid ──────────────────────────────────────────
final_IP_plot <- wrap_plots(plot_IP_sensitivity, ncol = 2) +
  plot_annotation(
    title    = "C",
    #  subtitle = "Black = Observed; Colored lines = Predicted at each VE value",
    theme    = theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40")
    )
  )

# Add shared y-axis label
y_axis <- ggdraw() +
  draw_label(
    "Weekly incident hospitalizations per 100,000 persons",
    angle = 90, size = 14
  )

final_IP_fig <- plot_grid(
  y_axis, final_IP_plot,
  ncol       = 2,
  rel_widths = c(0.05, 1)
)

final_IP_fig


(25 + 12)/2
(10 + 17)/2

#####################################################################################################
###### sensitivity analysis with alternative model structure 


# State-specific cumulative proportions vaccinated
state_uptake <- data.frame(
  Location = c("CA","CO","CT","GA","MI","MN","NM","NY","OH","OR","TN","UT"),
  annual_uptake = c(0.20926546, 0.20093509, 0.20937380, 0.09395524,
                    0.16611511, 0.27244292, 0.19401042, 0.15749792,
                    0.14679832, 0.27242121, 0.09924864, 0.14432907))
state_uptake <- data.frame(
  Location = c("CA","CO","CT","GA","MI","MN","NM","NY","OH","OR","TN","UT"),
  annual_uptake = c(0.5, 0.5, 0.5, 0.5,
                    0.5, 0.5, 0.5, 0.5,
                    0.5, 0.5, 0.5, 0.5))

weeks <- 1:52

# Logistic curve: cumulative proportion vaccinated
logistic_curve <- function(weeks, peak_week, max_prop, k = 0.6) {
  max_prop / (1 + exp(-k * (weeks - peak_week)))
}

# Build cumulative curves

# Peak of vaccination one month before outbreak
vax_cumulative <- state_uptake %>%
  crossing(week = weeks) %>%
  mutate(
    peak32    = logistic_curve(week, peak_week = 32, max_prop = annual_uptake),
    peak43    = logistic_curve(week, peak_week = 43, max_prop = annual_uptake),
    # Dual peak: 2x total coverage, split evenly across both waves
    dual_peak = logistic_curve(week, peak_week = 32, max_prop = annual_uptake) +
      logistic_curve(week, peak_week = 43, max_prop = annual_uptake),
    no_vax    = 0
  )

# Convert cumulative -> incident (weekly new vaccinations)
# For week 1, incident = cumulative; for subsequent weeks, take the difference
vax_incident <- vax_cumulative %>%
  arrange(Location, week) %>%
  group_by(Location) %>%
  mutate(
    inc_peak32    = c(peak32[1],    diff(peak32)),
    inc_peak43    = c(peak43[1],    diff(peak43)),
    inc_dual_peak = c(dual_peak[1], diff(dual_peak)),
    inc_no_vax    = 0
  ) %>%
  ungroup() %>%
  dplyr::select(Location, annual_uptake, week,
                inc_peak32, inc_peak43, inc_dual_peak, inc_no_vax) %>%
  pivot_longer(cols = starts_with("inc_"),
               names_to = "scenario",
               values_to = "incident_vaccinated") %>%
  mutate(scenario = recode(scenario,
                           "inc_peak32"    = "Summer dose",
                           "inc_peak43"    = "Fall dose",
                           "inc_dual_peak" = "Two dose",
                           "inc_no_vax"    = "No Vaccination"
  ))


# Plot incident vaccination curves
p <- ggplot(vax_incident %>% filter(scenario != "No Vaccination"),
            aes(x = week, y = incident_vaccinated,
                color = scenario, linetype = scenario)) +
  geom_line(linewidth = 0.8) +
  facet_wrap(~ Location, ncol = 4, scales = "free_y") +
  scale_color_manual(values = c("steelblue", "tomato", "forestgreen")) +
  labs(
    # title = "Incident (Weekly New) Vaccinations by State and Scenario",
    x = "Week",
    y = "Incident prop vaccinated",
    color = "Scenario", linetype = "Scenario"
  ) +
  theme_bw() +
  theme(
    legend.position = "bottom",
    #  strip.text = element_text(face = "bold"),
    legend.title = element_text(face = "bold"))
print(p)

################################################################
# set up model 

seirs_model_vax_scenario = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    hum = param["hum"]
    hum_scale = param["hum_scale"]
    temp = param["temp"]
    
    # parameters for underlying dynamics 
    beta = param["beta"]
    lp = param["lp"]
    gamma = param["gamma"]
    imm = param["imm"]
    imm_vax = param["imm_vax"]
    N = param["N"]
    VE = param["VE"]
    hr1 = param["hr"]
    hr2 = param["hr2"]
    
    # Empirical time series 
    vax_rate_at_t = dat_simulated_vax$incident_vaccinated[t] 
    humidity_at_t = dat_simulated_vax$humid_smooth[t] 
    temp_at_t = dat_simulated_vax$inverse_temp[t]
    
    seasonal = beta*( (hum_scale*(humidity_at_t  - 40)^2 + hum) + temp*temp_at_t  ) 
    
    imm_vec = imm  * dat_simulated_vax$immunity_scalar_fitted[t]
    
    
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


# 
# Vax incident 
head(vax_incident)
head(full_calibration_dat)

scenarios_vax_incident = vax_incident %>%
  mutate(year = 2023)

#states <- c("CA", "CO", "GA", "MI", "MN", "NY", "OH", "OR", "TN", "UT")
#best_beta  <- c(0.6253183, 0.40, 0.4818910, 0.4051849, 0.4093252, 0.34, 0.3852414, 0.56, 0.4829311, 0.5815605)
#best_S0    <- c(0.3384440, 0.70, 0.3106296, 0.4933294, 0.4776767, 0.6261234, 0.4604553, 0.2897822, 0.3468149, 0.2886504)
#best_h1    <- c(0.006410284, 0.004994745, 0.010281233, 0.008362320, 0.011041208, 0.01318657, 0.009762007, 0.007933405, 0.005509213, 0.008318385)
#best_h2    <- c(0.002893338, 0.001985926, 0.002474824, 0.002983213, 0.003026843, 0.004528284, 0.002192277, 0.002127474, 0.001989070, 0.0019318123)
#time_start <- c(12, 12, 6, 12, 12, 6, 10, 14, 10, 8)



states = c("CA", "CO", "GA", "MI", "MN", "NY", "OH", "OR", "TN", "UT")
best_beta = c(.617, .418, .459, .507, 0.435, .461, 0.4765, .614, .66, 0.56)
best_S0 = c( .276, .681, .358, .420, 0.512, 0.413, 0.405, .199, .32, 0.34)
best_h1 = c(0.00969,  0.0050, 0.0087, 0.0085, 0.0093, 0.013, 0.00697, .01299, .0051, 0.0081)
best_h2 = c(0.0027, 0.0021,  0.0024, 0.0025, 0.0028, 0.004, 0.00194, 0.00154,.0017, 0.00108 )
time_start = c(10, 14, 8, 12, 12, 6, 12, 14, 13, 10)

vax_scenarios <- unique(scenarios_vax_incident$scenario)  # pull scenario names from your data

# Results container: one row per state × scenario
results_list <- vector("list", length(vax_scenarios) * length(states))
result_idx   <- 1

for (p in vax_scenarios) {
  
  # Filter incident vaccination data to this scenario only
  vax_scenario_data <- scenarios_vax_incident %>%
    filter(scenario == p)
  
  for (j in seq_along(states)) {
    
    cat("Running scenario:", p, "| State:", states[j], "\n")
    
    dat_simulated_vax <- full_calibration_dat %>%
      filter(Location == states[j]) %>%
      filter(t > time_start[j]) %>%
      mutate(t = row_number()) %>%
      mutate(immunity_scalar_fitted = ifelse(slope_pos > 0, 1.46, 1)) %>%
      mutate(date = ISOweek2date(
        paste0(year, "-W", sprintf("%02d", week), "-1"))) %>%
      left_join(vax_scenario_data %>% dplyr::select(week, year, Location, incident_vaccinated),
                by = c("week", "year", "Location")) %>%
      mutate(incident_vaccinated = replace(incident_vaccinated, is.na(incident_vaccinated), 0))
    
    N <- print(unique(dat_simulated_vax$pop_size_2023))
    S0 = best_S0[j]*N
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
    time = seq(1, nrow(dat_simulated_vax), 1)
    
    # Parameters
    param <- c(
      hum =   6.01, 
      hum_scale = .0011, 
      temp = .107,  
      imm = .0506, 
      imm_vax  = .1188,
      hr = best_h1[j],
      hr2 =  best_h2[j],
      beta = best_beta[j],
      lp = 1/(4/7),
      gamma = 1/(4/7),              
      VE = .30, 
      N =  unique(dat_simulated_vax$pop_size_2023))
    
    out_sim <- ode(
      y     = state, times = time, func  = seirs_model_vax_scenario, parms = param )
    
    out_df_sim <- as.data.frame(out_sim) %>%
      mutate(predicted_hosp = I * best_h1[j] + I2 * best_h2[j]) %>%
      mutate(t = time) %>%
      left_join(
        dat_simulated_vax %>% dplyr::select(t, date, incident),
        by = "t"
      ) %>%
      filter(date > "2023-08-01" & date < "2024-02-01")
    
    results_list[[result_idx]] <- data.frame(
      scenario        = p,
      state           = states[j],
      predicted_hosp  = sum(out_df_sim$predicted_hosp, na.rm = TRUE) )
    
    result_idx <- result_idx + 1
  }
}

# Combine all results into a single tidy data frame
results_df <- bind_rows(results_list)
head(results_df)
print(unique(results_df$scenario))

baseline_vax_scenario = results_df %>%
  filter(scenario == "No Vaccination") %>%
  mutate(baseline_hosp = predicted_hosp) %>% 
  dplyr::select(-predicted_hosp, - scenario)

relative_vax_scenario = left_join(results_df %>% filter(scenario != "No Vaccination"), baseline_vax_scenario ,
                                  by = c("state")) %>%
  mutate(change = (baseline_hosp - predicted_hosp)/baseline_hosp) %>%
  mutate(percent_change = change*100 ) %>%
  mutate(scenario = factor(scenario, levels = c("Fall dose", "Summer dose", "Two dose")))

winter_order <- relative_vax_scenario %>%
  filter(scenario == "Fall dose") %>%
  arrange(desc(percent_change)) %>%
  pull(state)

relative_vax_scenario <- relative_vax_scenario %>%
  mutate(state = factor(state, levels = winter_order))

dosage = ggplot(data = relative_vax_scenario )   +
  geom_col(aes(x = state, y = percent_change, fill = scenario), position = "dodge") +
  ylab("Percent hospitalizations averted \ncompared to no vaccination (%)") + theme_bw() +
  xlab("Location") +
  scale_fill_manual(values = c("blue4", "orangered", colorspace::lighten("aquamarine1", 0.4) ), 
                    guide = guide_legend(nrow = 1), 
                    name = 'Strategy')  +
  theme(legend.position = "bottom") +
  theme(
    plot.title      = element_text(size = 16, face = "bold"),   # title
    axis.title.x    = element_text(size = 14),   # x-axis title
    axis.title.y    = element_text(size = 14),   # y-axis title
    axis.text.x     = element_text(size = 14),   # x-axis tick labels
    axis.text.y     = element_text(size = 14) ,   # y-axis tick labels
    legend.text = element_text(size = 14), legend.title = element_text(size = 14), 
    strip.text = element_text(size = 14) ) +
  # scale_y_continuous(
  #  limits = c(0, 11),
  #  breaks = c(0, 2, 4,6,8, 10, 12) )  +
  scale_y_continuous(
    limits = c(0, 6),
    breaks = c(0, 2, 4, 6) ) # +
# ggtitle("D")
dosage


###################################################
# density plots 

chains = read.csv("mcmc_chains_fast_genbeta_feb27.csv")
head(chains)


chains_long <- chains %>%
  pivot_longer(
    cols = -c(chain, iteration),   # keep chain + iteration fixed
    names_to = "variable",
    values_to = "value" ) %>%
  # filter(iteration > 2000) %>%
  filter(variable %in% c("imm", "imm_vax", "hum", "hum_scale", "temp", "var_scalar")) %>%
  mutate(variable = replace(variable, variable == "imm", 'Infection-derived immunity')) %>%
  mutate(variable = replace(variable, variable == "imm_vax", 'Vaccine-derived immunity')) %>%
  mutate(variable = replace(variable, variable == "hum", 'Humidity scalar (b2)'))%>%
  mutate(variable = replace(variable, variable == "hum_scale", 'Humidity scalar (b1)'))%>%
  mutate(variable = replace(variable, variable == "temp", 'Temperature scalar (b3)')) %>%
  mutate(variable = replace(variable, variable == "var_scalar", 'Variant scalar (g(t))'))
head(chains_long)
print(unique(chains_long$chain))
print(unique(chains_long$variable))

summary_stats = chains_long %>%
  group_by(variable) %>%
  mutate(val_025 = quantile(value, 0.025, na.rm = TRUE), 
         val_975 = quantile(value, 0.975, na.rm = TRUE), 
         sd = sd(value)) %>%
  ungroup() %>%
  distinct(variable, val_025, val_975, sd)
head(summary_stats)

density_plots = ggplot(chains_long, aes(x = round(value, digits = 3), fill = factor(chain), color = factor(chain))) +
  geom_density(alpha = 0.3, position = "identity", adjust = 6) +
  facet_wrap(vars(variable), scales = "free", ncol = 3) +
  labs(fill = "Chain", color = "Chain") +
  theme_minimal() +
  scale_color_manual(values = c("gray36", "darkslategray3", "darkslategray4")) +
  scale_fill_manual(values = c("gray36", "darkslategray3", "darkslategray4")) +
  xlab("Variable") + ylab("Density") + ggtitle("Posteriors") +
  theme(
    plot.title      = element_text(size = 16, face = "bold"), 
    axis.title.x    = element_text(size = 14),   # x-axis title
    axis.title.y    = element_text(size = 14),   # y-axis title
    strip.text = element_text(size = 14) ) 
density_plots


library(dplyr)
library(tidyr)
library(ggplot2)
library(purrr)

# ---- Reshape chains ----
chains_long <- chains %>%
  filter(iteration > 500) %>%
  pivot_longer(
    cols = -c(chain, iteration),
    names_to = "variable",
    values_to = "value"
  ) %>%
  filter(variable %in% c("imm", "imm_vax", "hum", "hum_scale",
                         "temp", "var_scalar")) %>%
  mutate(
    variable = recode(variable,
                      imm        = "Infection-derived immunity",
                      imm_vax    = "Vaccine-derived immunity",
                      hum        = "Humidity scalar (b2)",
                      hum_scale  = "Humidity scalar (b1)",
                      temp       = "Temperature scalar (b3)",
                      var_scalar = "Variant scalar (g(t))"
    )
  )

# ---- Compute densities per chain × variable ----
dens_df <- chains_long %>%
  group_by(variable, chain) %>%
  group_modify(~{
    d <- density(.x$value, adjust = 6)
    tibble(x = d$x, y = d$y)
  }) %>%
  ungroup()

# ---- Scale by max density per variable (across chains) ----
dens_df <- dens_df %>%
  group_by(variable) %>%
  mutate(y_scaled = y / max(y)) %>%
  ungroup()

# ---- Plot ----
density_plots <- ggplot(dens_df,
                        aes(x = x,
                            y = y_scaled,
                            fill = factor(chain),
                            color = factor(chain))) +
  geom_line(linewidth = 1) +
  geom_area(alpha = 0.3, position = "identity") +
  facet_wrap(vars(variable), scales = "free", ncol = 3) +
  labs(fill = "Chain",
       color = "Chain",
       x = "Value",
       y = "Density",
       title = "Posteriors") +
  theme_minimal() +
  scale_color_manual(values = c("gray36", "darkslategray3", "darkslategray4")) +
  scale_fill_manual(values = c("gray36", "darkslategray3", "darkslategray4")) +
  theme(
    plot.title   = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    strip.text   = element_text(size = 14)
  )

density_plots


################################################
# CHains 

head(chains_long)
chains_plot = ggplot(data = chains_long) +
  geom_line(aes(x = iteration, y = value, col = as.character(chain)), alpha = .8, lwd = 1.3) + facet_wrap(vars(variable), scales = "free", 
                                                                                                          ncol = 3) +
  theme_minimal() +
  scale_color_manual(values = c( "gray36", "darkslategray3", "darkslategray4" ))  +
  labs(col = "Chain", 
       x = "Iteration")  + 
  ylab("Value") 
chains_plot

y_limits <- chains_long %>%
  group_by(variable, chain) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    sd_val   = sd(value, na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  group_by(variable) %>%
  summarise(
    ymin = min(mean_val - 3 * sd_val),
    ymax = max(mean_val + 3 * sd_val),
    .groups = "drop"
  )

chains_long <- chains_long %>%
  left_join(y_limits, by = "variable")

chains_plot = ggplot(data = chains_long) +
  geom_blank(aes(x = iteration, y = ymin)) +
  geom_blank(aes(x = iteration, y = ymax)) +
  geom_line(aes(x = iteration, y = value, col = as.character(chain)), alpha = .8, lwd = 1.3) +
  facet_wrap(vars(variable), scales = "free", ncol = 3) +
  theme_minimal() +
  scale_color_manual(values = c("gray36", "darkslategray3", "darkslategray4")) +
  labs(col = "Chain", x = "Iteration") +
  ylab("Value") + ggtitle("Iterations")

chains_plot 

plot_grid(density_plots, chains, nrow = 2)
