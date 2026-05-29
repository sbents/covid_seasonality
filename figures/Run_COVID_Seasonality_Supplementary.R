# Run Supplementary Analyses

library(deSolve)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)
library(zoo)
library(ISOweek)

# VALIDATION on 2025 
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

# Model 
seirs_model_hybrid_state_validation = function(t, state, param) {
  with(as.list(state), {
    
    # parameters for seasonality multipliers 
    hum = param["hum"]
    temp = param["temp"]
    
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
  
  plot_best_fits[[m]]  = ggplot(data = out_df_sim) +
    geom_line(aes(x = date, y = rollmean(check_pred, k = 1, na.pad = TRUE),col = "Predicted" ),  lwd = 1.5) +
    geom_line(aes(x = date, y = rollmean(check_real, k = 1, na.pad = TRUE), col ="Observed"), lwd = 1.2,  alpha = 1) +
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
    theme(legend.position = "none") +
    geom_vline(xintercept = as.Date("2024-11-15"), lty = "dashed")}

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

# VE
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
    pop  = N ) %>%
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
      N          = N )
    
    out <- as.data.frame(
      ode(y = state_init, times = time,
          func = seirs_model_hybrid_state_fits, parms = param)) %>%
      mutate(
        predicted_hosp = I * best_h1[m] + I2 * best_h2[m],
        check_pred     = (predicted_hosp / N) * 100000,
        VE_label       = as.character(ve),
        date           = fits$date ) %>%
      filter(date > "2022-04-01")
      out
  })
  
  sim_df  <- bind_rows(sim_list)
  obs_df2 <- obs_df %>% filter(date > "2022-04-01")
  
  plot_VE_sensitivity[[m]] <- ggplot() +
    geom_line(
      data = sim_df,
      aes(
        x     = date,
        y     = rollmean(check_pred, k = 4, na.pad = TRUE),
        color = VE_label,
        group = VE_label
      ),
      lwd = 1.3, alpha = 0.85 ) +
    scale_color_manual(
      values = VE_colors,
      name   = "VE",
      labels = paste0("VE = ", VE_values)) +
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
      legend.text     = element_text(size = 14) ) +
    scale_y_continuous(
      limits = c(0, 20),
      breaks = c(0, 5, 10, 15, 20)  ) +
    guides(color = guide_legend(nrow = 1))
}

# ── Combine into grid ──────────────────────────────────────────
final_VE_plot <- wrap_plots(plot_VE_sensitivity, ncol = 2) +
  plot_annotation(
    theme    = theme(
      plot.title    = element_text(size = 16, face = "bold"),
      plot.subtitle = element_text(size = 12, color = "gray40") ))

# Add shared y-axis label
y_axis <- ggdraw() +
  draw_label(
    "Weekly incident hospitalizations per 100,000 persons",
    angle = 90, size = 14 )

final_VE_fig <- plot_grid(
  y_axis, final_VE_plot,
  ncol       = 2,
  rel_widths = c(0.05, 1))

final_VE_fig

# states 
states     <- c("CA", "NY")
best_beta  <- c(0.6253183, 0.34)
best_S0    <- c(0.3384440, 0.6261234)
best_h1    <- c(0.006410284,  0.01318657)
best_h2    <- c(0.002893338,  0.004528284)
time_start <- c(12,  6 )


#######################################################
# Test latent period sensitivity 
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

  obs_df <- data.frame(
    t    = time,
    date = fits$date,
    hosp = fits$hosp,
    pop  = N
  ) %>%
    mutate(check_real = (hosp / pop) * 100000)

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
      N          = N  )
    
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
    geom_line(
      data = sim_df,
      aes(
        x     = date,
        y     = rollmean(check_pred, k = 2, na.pad = TRUE),
        color = LP_label,
        group = LP_label ),
      lwd = 1.1, alpha = 0.85 ) +
    scale_color_manual(
      values = LP_colors,
      name   = "Latent Period (days)",
      labels = paste0("LP = ", round(7 / LP_values, 1), " days") ) +
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


# Sensitivity analysis with alternative model structure 
#########################################################

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
  tidyr::crossing(week = weeks) %>%
  mutate(
    peak32    = logistic_curve(week, peak_week = 32, max_prop = annual_uptake),
    peak43    = logistic_curve(week, peak_week = 43, max_prop = annual_uptake),
    # Dual peak: 2x total coverage, split evenly across both waves
    dual_peak = logistic_curve(week, peak_week = 32, max_prop = annual_uptake) +
      logistic_curve(week, peak_week = 43, max_prop = annual_uptake),
    no_vax    = 0)

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
    strip.text = element_text(size = 14) ) 
dosage

# Density plots from MCMC fits 
###################################################
chains = read.csv("mcmc_chains_fast_genbeta_feb27.csv")

chains_long <- chains %>%
  pivot_longer(
    cols = -c(chain, iteration),   # keep chain + iteration fixed
    names_to = "variable",
    values_to = "value" ) %>%
  filter(variable %in% c("imm", "imm_vax", "hum", "hum_scale", "temp", "var_scalar")) %>%
  mutate(variable = replace(variable, variable == "imm", 'Infection-derived immunity')) %>%
  mutate(variable = replace(variable, variable == "imm_vax", 'Vaccine-derived immunity')) %>%
  mutate(variable = replace(variable, variable == "hum", 'Humidity scalar (b2)'))%>%
  mutate(variable = replace(variable, variable == "hum_scale", 'Humidity scalar (b1)'))%>%
  mutate(variable = replace(variable, variable == "temp", 'Temperature scalar (b3)')) %>%
  mutate(variable = replace(variable, variable == "var_scalar", 'Variant scalar (g(t))'))

summary_stats = chains_long %>%
  group_by(variable) %>%
  mutate(val_025 = quantile(value, 0.025, na.rm = TRUE), 
         val_975 = quantile(value, 0.975, na.rm = TRUE), 
         sd = sd(value)) %>%
  ungroup() %>%
  distinct(variable, val_025, val_975, sd)

chains_long <- chains %>%
  filter(iteration > 500) %>%
  pivot_longer(
    cols = -c(chain, iteration),
    names_to = "variable",
    values_to = "value" ) %>%
  filter(variable %in% c("imm", "imm_vax", "hum", "hum_scale",
                         "temp", "var_scalar")) %>%
  mutate(
    variable = recode(variable,
                      imm        = "Infection-derived immunity",
                      imm_vax    = "Vaccine-derived immunity",
                      hum        = "Humidity scalar (b2)",
                      hum_scale  = "Humidity scalar (b1)",
                      temp       = "Temperature scalar (b3)",
                      var_scalar = "Variant scalar (g(t))" ) )

dens_df <- chains_long %>%
  group_by(variable, chain) %>%
  group_modify(~{
    d <- density(.x$value, adjust = 6)
    tibble(x = d$x, y = d$y)}) %>%
  ungroup()

dens_df <- dens_df %>%
  group_by(variable) %>%
  mutate(y_scaled = y / max(y)) %>%
  ungroup()

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
    strip.text   = element_text(size = 14) )
density_plots

