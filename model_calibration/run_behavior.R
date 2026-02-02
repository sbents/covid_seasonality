#setwd("/Users/sambents/Desktop/Lo/covid_seasonality/data/processed")

# load packages 
library(deSolve)
library(ggplot2) 
library(lubridate)
library(dplyr)
library(tidyr)
library(zoo)
library(lhs)
library(stringr)
library(readr)
library(tidyverse)
library(tibble)
library(reshape2)
library(parallel)
library(future)
library(coda)
library(gridExtra)
library(cowplot)

#################################################################
# Comment out data processing 
# Data 
#Location = c( "CT", "GA", "MI",  "MN", "NY",  "OH", "TN", "CA", "CO", "NM", "OR", "UT")
#pop_size_2023 = c( 3643023,  11064432, 10083356, 5753048, 19737367, 11824034,  7148304, 
#                   39198693,   5901339,  2121164, 4253653, 3443222)
#census_dat = data.frame(Location, pop_size_2023) # %>%

# Load in data and check that hosp rate per 100k is proper! 
#full_calibration_dat = read.csv("processed_calibration_dat.csv")  %>%
#  filter(Location %in% c( "CT", "MN", "OH", "MI", "NY", "TN", "GA", 
#                          "CA", "CO", "NM", "OR", "UT")) %>%
#  mutate(humid_smooth_mult = ifelse(humid_smooth < 40, 1, 0 )) %>%
#  left_join(census_dat, by = "Location") %>%
#  mutate(hosp_per_100k = (hosp/ pop_size_2023 )*100000)
#head(full_calibration_dat )

# load in behavior data
#url <- "https://media.githubusercontent.com/media/bansallab/indoor_outdoor/refs/heads/main/indoor_activity_data/indoor_activity_2018_2020.csv"
#behavior_dat <- read.csv(url)
#head(behavior_dat)

#county_list = c(6001, 6013, 6075, 8001, 8005, 8031, 8035, 8059, 9009, 9007, 13097, 13067, 13089, 
#         13097, 13121, 13135, 13217, 13243, 26049, 26051, 26049, 26065, 26161, 27003, 
#         27019, 27037, 27053, 27123, 27139, 27163, 35001, 35005, 35013, 35021, 35029, 
#         35045, 35049, 36001, 36021, 36039, 36041, 36049, 36055, 36053, 36073, 36075, 
#         36083, 36091, 36093, 36095, 36119, 36135, 39041, 39045, 39049, 39075, 39089, 
#         39097, 39129, 39143, 39127, 39165, 41005, 41051, 41067, 47017, 47037, 47043, 
#         47137, 47157, 47165, 47187, 47189, 49035)

#county = c(6001, 6013, 6075, 8001, 8005, 8031, 8035, 8059, 9009, 9007, 13097, 13067, 13089, 
#                13097, 13121, 13135, 13217, 13243, 26049, 26051, 26049, 26065, 26161, 27003, 
#                27019, 27037, 27053, 27123, 27139, 27163, 35001, 35005, 35013, 35021, 35029, 
#                35045, 35049, 36001, 36021, 36039, 36041, 36049, 36055, 36053, 36073, 36075, 
#                36083, 36091, 36093, 36095, 36119, 36135, 39041, 39045, 39049, 39075, 39089, 
#                39097, 39129, 39143, 39127, 39165, 41005, 41051, 41067, 47017, 47037, 47043, 
#                47137, 47157, 47165, 47187, 47189, 49035)
#Location = c(rep("CA", 3), rep("CO", 5), rep("CT", 2), rep("GA", 8 ), rep("MI", 5), rep("MN", 7),
#          rep("NM", 7), rep("NY", 15), rep("OH", 10), rep("OR", 3), rep("TN", 8), rep("UT", 1))
#fips_county = data.frame(county, Location)
#head(fips_county)

#behavior_dat_clean = behavior_dat %>%
#  mutate(county = as.numeric(county)) %>%
#  filter(county %in% county_list) %>%
#  left_join(fips_county, by = c("county"), relationship = "many-to-many") %>%
#  group_by(Location, date) %>%
#  mutate(mean_indoor = mean(indoor_activity)) %>%
#  mutate(date = as.Date(date)) %>%
#  mutate(year = year(date), week = week(date)) %>%
#  filter(year == 2019) 
  
#head(behavior_dat_clean)

# plot out indoor activity 
#ggplot(data = behavior_dat_clean) +
#  geom_line(aes(x = as.Date(date), y = mean_indoor, col = Location)) +
#  theme_bw() +
#  facet_wrap(vars(Location))

#full_calibration_dat_activity = left_join(full_calibration_dat, behavior_dat_clean %>% distinct(Location, week, mean_indoor), 
#                                          by = c("week", "Location")) %>%
#  mutate(mean_indoor = if_else(is.na(mean_indoor), 1, mean_indoor)) %>%
#  mutate(mean_indoor = mean_indoor*100)
#head(full_calibration_dat_activity)

#write.csv(full_calibration_dat_activity, "full_calibration_dat_activity.csv")


full_calibration_dat_activity = read.csv("full_calibration_dat_activity.csv")

# MCMC ###################################################################################

run_hierarchical_mcmc_calibration <- function(full_calibration_dat_activity, iterations = 18000, burn_in = 8000, thin = 2, n_chains = 3) {
  
  # check states to be fit in the calibration
  states <- unique(full_calibration_dat_activity$Location)
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
    "CA" =  39198693,
    "CO" = 5901339, 
    "NM" = 2121164, 
    "OR" = 4253653, 
    "UT" = 3443222
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
    act = parms[["act"]]
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
    state_data <- parms[["state_data"]] # %>%
     # mutate(humid_smooth = ifelse(humid_smooth < 40, humid_smooth* hum_scale, humid_smooth))
    
    # time, specify it must be less than the number of rows in calibration data 
    t_idx <- round(t)
    if(t_idx < 1) t_idx <- 1
    if(t_idx > nrow(state_data)) t_idx <- nrow(state_data)
    
    # pull in time series influencing beta 
    vax_rate_at_t <- state_data$incident[t_idx]
   # humidity_at_t <- state_data$humid_smooth[t_idx]
    act_at_t <- state_data$mean_indoor[t_idx]
    
    
    # seasonal function
    seasonal <- beta * ( act* act_at_t )
    
    # make seasonal value something small and positive if negative- this should not occur to start with 
    if(!is.finite(seasonal) || seasonal < 0) seasonal <- 0.01
    
    #  diff eq's 
    dS <- -seasonal * S * (I + I2) / N - VE * vax_rate_at_t * S + imm_vax * V
    dE <- seasonal * S * (I + I2) / N - lp * E
    dI <- lp * E - gamma * I
    dR <- I * gamma - R * imm
    dV <- VE * vax_rate_at_t * S - imm_vax * V
    
    dS2 <- -seasonal * S2 * (I + I2) / N + imm * R + imm * R2 + imm_vax * V2 - VE * vax_rate_at_t * S2
    dE2 <- seasonal * S2 * (I + I2) / N - lp * E2
    dI2 <- lp * E2 - gamma * I2
    dR2 <- I2 * gamma - R2 * imm
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
      state_data <- full_calibration_dat_activity[full_calibration_dat_activity$Location == state_name, ]
      state_data <- state_data[order(state_data$t), ]
      
      # state-specific parameters
      #  beta_param <- paste0("beta_", state_name)
      #  beta_value <- params[[beta_param]]
      beta_value <- params[["beta"]]
      
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
      #  hum_scale = params[["hum_scale"]] , 
      #  hum = params["hum"], 
        act = params["act"], 
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
      
      predicted_hosp <- round((((out_df$I * all_params$hr1 + out_df$I2 * all_params$hr2)/N_val)*100000), digits = 0)
      
      # observed data 
      observed_hosp <- round((((state_data$hosp)/N_val)*100000), digits = 0)
      
    #  predicted_hosp <- out_df$I * all_params$hr1 + out_df$I2 * all_params$hr2
      
      # observed data 
    #  observed_hosp <- state_data$hosp
      
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
  
  
  
  # prior function 
  log_prior <- function(params) {
    
    # global parameters 
    global_bounds <- list(
      # var = c(0.05, .09),
      act = c(.05, 200),
      imm = c(1/50, 1/19),
      imm_vax = c(1/20, 1/8),
      hr1 = c(0.0001, 0.05),
      hr2 = c(0.0001, 0.05), 
      S0 = c(0.30, 0.80), 
      beta = c(.2, .7)
    )
    
    # state-specific parameters 
    #  state_bounds <- list(
    #     beta = c(0.1, 0.70)
    #  )
    
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
    #  for(state_name in states) {
    #    # Check beta
    #    beta_param <- paste0("beta_", state_name)
    #    if(beta_param %in% names(params)) {
    #      bounds <- state_bounds[["beta"]]
    #      if(params[beta_param] < bounds[1] || params[beta_param] > bounds[2]) {
    #        return(-Inf)
    #      }
    #    }
    
    #  }
    
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
    # var = .07,
  #  hum_scale = 1.73, 
   # hum = 0.27,
    act = .89, 
    imm = .044,
    imm_vax = .050,
    hr1 = 0.047,
    hr2 = 0.002, 
    S0 = 0.567, 
    beta = .344
  )
  
  
  # state-specific initial values (only fitted parameters)
  # state_initial <- list()
  # for(state_name in states) {
  #   state_initial[[paste0("beta_", state_name)]] <- runif(1, 0.2, .45)
  # }
  
  # combine initial parameters 
  # initial_params <- c(global_initial, state_initial)
  initial_params <- c(global_initial)
  param_names <- names(initial_params)
  
  cat("=== PARAMETER SUMMARY ===\n")
  cat("Total FITTED parameters:", length(param_names), "\n")
  cat("  Global fitted parameters:", length(global_initial), "\n")
  cat("    - beta (transmission rate): GLOBAL (shared across all states)\n")
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
    
    fixed_scales <- abs(current_params) * 0.16  # FIXED 16% proposal
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
  
  # return(list(
  #    samples = all_samples,
  #    mcmc_list = mcmc_list,
  #    acceptance_rates = acceptance_rates,
  #    gelman_diag = gelman_diag,
  #    summary = summary_stats,
  #    states = states,
  #    param_names = param_names,
  #    final_proposal_covs = lapply(chain_results, function(x) x$final_proposal_cov)
  #  ))
  
  return(list(
    samples = all_samples,
    mcmc_list = mcmc_list,
    acceptance_rates = acceptance_rates,
    gelman_diag = gelman_diag,
    summary = summary_stats,
    states = states,
    param_names = param_names,
    final_proposal_covs = lapply(chain_results, function(x) x$final_proposal_cov),
    likelihood_func = likelihood_func,      # ADD THIS LINE
    log_posterior = log_posterior,          # ADD THIS LINE
    full_calibration_dat_activity = full_calibration_dat_activity  # ADD THIS LINE (needed by likelihood_func)
  ))
}



#############################################################
### Run the hierarchical MCMC calibration 
results <- run_hierarchical_mcmc_calibration(
  full_calibration_dat_activity = full_calibration_dat_activity,  # your data
  iterations = 14000,     # start small for testing
  burn_in = 0,         # burn-in period
  thin = 1,              # thinning interval
  n_chains = 3 )          # number of MCMC chains


#############################################################
### Post-processing MCMC results 

# Combine all chains into a single dataframe with chain ID
all_chains_df <- do.call(rbind, lapply(1:length(results$samples), function(chain_id) {
  chain_data <- as.data.frame(results$samples[[chain_id]])
  chain_data$chain <- chain_id
  chain_data$iteration <- 1:nrow(chain_data)
  return(chain_data)
}))

# Save to CSV
write.csv(all_chains_df, "mcmc_chains_immact_genbeta_jan15.csv", row.names = FALSE)



sink("mcmc_summary_immact_genbeta_jan15.txt")

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

# Calculate final likelihood from posterior samples
cat("\n=== FINAL LIKELIHOOD ===\n")

# Combine all samples from all chains
all_samples_combined <- do.call(rbind, results$samples)

# Option 1: Use mean of all posterior samples
cat("\n--- Mean of all posterior samples ---\n")
mean_params <- colMeans(all_samples_combined)
names(mean_params) <- results$param_names

cat("Mean parameter values across all chains:\n")
print(mean_params)

mean_likelihood <- results$likelihood_func(mean_params)
cat("\nLog-likelihood at mean parameters:", mean_likelihood, "\n")

# Option 2: Get maximum a posteriori (MAP) estimate
#cat("\n--- Maximum a posteriori (MAP) estimate ---\n")
# Calculate log posterior for each sample
#log_posts <- apply(all_samples_combined, 1, function(params) {
#  names(params) <- results$param_names
#  results$log_posterior(params)
#})

#map_index <- which.max(log_posts)
#map_params <- all_samples_combined[map_index, ]
#names(map_params) <- results$param_names

#cat("MAP parameter values:\n")
#print(map_params)

#map_likelihood <- likelihood_func(map_params)
#cat("\nLog-likelihood at MAP parameters:", map_likelihood, "\n")
#cat("Log-posterior at MAP parameters:", log_posts[map_index], "\n")

cat("=== END LIKELIHOOD SUMMARY ===\n")

# Close the sink
sink()








