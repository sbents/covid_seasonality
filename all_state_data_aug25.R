# ALL STATE CALIBRATION

###############################################################
# Load the necessary library
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
library(cowplot)

###############################################################
# Set working directory 
setwd("/Users/sambents/Desktop/Lo/covid_seasonality")

# DATA ########################################################
# Hospitalization data  #######################################
Location = c("CA", "CO", "CT", "GA", "MD",
             "MI", "MN", "NM", "NY", "OH", 
             "OR", "TN", "UT")
Site = c("California", "Colorado", "Connecticut", 
         "Georgia", "Maryland", "Michigan", "Minnesota", 
         "New Mexico", "New York", "Ohio", "Oregon",
         "Tennessee", "Utah")
location_site = data.frame(Location, Site)

net = read.csv("data/hosp/RSV_net.csv")  %>%
  mutate(week = as.Date(Week.Ending.Date)) %>%
  filter(Surveillance.Network == "COVID-NET") %>%
  left_join(location_site, by = "Site")
head(net)

# Visualize 
ggplot(data = net %>% filter(Site != "Overall") %>% filter(as.Date(Week.Ending.Date) > "2022-03-01")
       %>% filter(Site != "Iowa")  %>% filter(Site != "North Carolina") ) +
  geom_line(aes(x = as.Date(Week.Ending.Date), y = Weekly.Rate), lwd = 1.2) +
  theme_bw() + 
  facet_wrap(vars(Location), scales = "free_y") +
  ylab("Weeky hospitalization rate per 100k") +
  ggtitle("COVID-NET") # +
 
# geom_vline(xintercept = as.Date("2022-08-01"), col = "red")+
 # geom_vline(xintercept = as.Date("2023-01-01"), col = "red")

# Vaccination data  ###########################################################
# Filter for COVID-Net states 
# https://data.cdc.gov/Vaccinations/COVID-19-Vaccinations-in-the-United-States-Jurisdi/unsk-b7fc/about_data
covid_vax = read.csv("data/vaccine/covid_vax.csv") %>%
  mutate(date = as.Date(Date, format = "%m/%d/%Y")) %>% 
  filter(Location %in% c("CA", "CO", "CT", "GA", "MD",
                         "MI", "MN", "NM", "NY", "OH", "OR", "TN", "UT")) %>%
  dplyr::select(date, MMWR_week, Location, Administered, Admin_Per_100K, Distributed) %>%
  arrange(Location, date) %>%
  mutate(vax_admin_incident = Administered - lag(Administered))
head(covid_vax)

# Population Size, use 2023 estimate from census 
# https://www.census.gov/data/tables/time-series/demo/popest/2020s-state-total.html
Location = c("CA", "CO", "CT", "GA", "MD",
             "MI", "MN", "NM", "NY", "OH", 
             "OR", "TN", "UT")
pop_size_2023 = c(39198693, 5901339, 3643023, 11064432, 6217062,
                  10083356, 5753048, 2121164, 19737367, 11824034,
                  4253653, 7148304, 3443222)
census_dat = data.frame(Location, pop_size_2023) 
head(census_dat)

covid_vax_rate = left_join(covid_vax, census_dat, by = "Location") %>%
  mutate(vax_rate = (vax_admin_incident/pop_size_2023)) %>%
  filter(date > "2021-11-15") %>%
  mutate(month = month(date), year = year(date), week = week(date)) %>%
  group_by(week, year, Location) %>%
  mutate(week_vax_rate = sum(vax_rate, na.rm = TRUE)) %>%
  ungroup() %>%
  distinct(Location, week_vax_rate, year, week, week_vax_rate) %>%
  mutate(smooth_vax_rate = rollmean(as.numeric(week_vax_rate), k = 3, fill = NA)) %>%
  mutate(smooth_vax_rate = ifelse(smooth_vax_rate < 0, 0, smooth_vax_rate)) 

vaccine_check1 = covid_vax_rate %>%
  distinct(year, Location, week, week_vax_rate) %>%
  group_by(year, Location) %>%
  mutate(annual_vax_rate = sum(week_vax_rate)) %>%
  distinct(year, Location, annual_vax_rate)

ggplot(data = vaccine_check1, aes(x = as.character(year), y = annual_vax_rate, col = Location)) +
  geom_point() + xlab("Year") + theme_bw()

# Visualize weekly incident 
ggplot(data = covid_vax_rate) +
  geom_line(aes(x= week/52 + year, y = smooth_vax_rate , col = Location), lwd = 2) +
  theme_bw() +
  facet_wrap(vars(Location)) +
  scale_color_viridis_d() +
  ylab("Weekly vaccination rate (%)") +
  ggtitle("CDC Vaccination data: 11/15/2021- 5/10/2023")  
head(covid_vax_rate)
  
# Visualize monthly incident 
#ggplot(data = covid_vax_rate) +
#  geom_smooth(aes(x= date, y = month_vax_rate, col = Location), span = .19, lwd = 2, se = FALSE) +
#  theme_bw() +
#  scale_color_viridis_d() +
#  ylab("Monthly vaccination rate (%)") +
#  ggtitle("CDC Vaccination data: 11/15/2021- 5/10/2023")  

# More vaccination data 
# https://data.cdc.gov/Vaccinations/Monthly-Cumulative-Number-and-Percent-of-Persons-W/3bmy-cyyd/about_data

month_code <- c(
  JAN = "01", FEB = "02", MAR = "03", APR = "04",
  MAY = "05", JUN = "06", JUL = "07", AUG = "08",
  SEP = "09", OCT = "10", NOV = "11", DEC = "12"
)

covid_vax_2425 = read.csv("data/vaccine/monthly_covid_vax.csv") %>%
  filter(Jurisdiction %in% c("Colorado", "California", "Minnesota" ,
                             "Oregon", "Michigan", "Maryland", "Ohio",
                             "New York State", "Georgia", "Connecticut",
                             "New Mexico", "Tennessee" , "Utah")) %>%
  filter(Age_group_label == "Overall") %>%
  mutate(
    start_year = as.numeric(substr(current_season, 1, 4)),
    end_year   = as.numeric(substr(current_season, 6, 7)) + 2000,
    month_num  = month_code[Month],
    year       = ifelse(Month %in% c("AUG","SEP","OCT","NOV","DEC"),
                        start_year, end_year),
    year_month = paste0(year, "-", month_num, "-01")
  ) %>%
  mutate(year_month = as.Date(year_month)) %>%
  arrange(Jurisdiction, year_month) %>%
  mutate(incident = Estimate - lag(Estimate)) %>%
  mutate(incident = ifelse(incident < 0, 0, incident)) %>%
  mutate(incident = ifelse(is.na(incident), 0, incident)) %>%
  filter(Jurisdiction != "Maryland")
head(covid_vax_2425)

vaccine_check2 = covid_vax_2425 %>%
  distinct(Jurisdiction, year_month, incident, year) %>%
  group_by(year, Jurisdiction) %>%
  mutate(annual_vax_rate = sum(incident)) %>%
  distinct(year,Jurisdiction, annual_vax_rate)

# Visualize 
ggplot(data = covid_vax_2425, aes(x = as.Date(year_month), y = incident*100, col = Jurisdiction)) +
  geom_line(lwd = 2) +
  theme_bw() +
  scale_color_viridis_d() +
  ylab("Monthly percent of pop vaccinated (%)") +
  ggtitle("CDC Vaccination data: 07/01/2023 -7/01/2025") +
  xlab("Date") #+ 
 # facet_wrap(vars(Jurisdiction))

ggplot(data = vaccine_check2, aes(x = as.character(year), y = annual_vax_rate, col = Jurisdiction)) +
  geom_point() + xlab("Year") + theme_bw()

# Temperature data  ###########################################################
# states to read in
Location = c("CA", "CO", "CT", "GA", "MD",
             "MI", "MN", "NM", "NY", "OH", 
             "OR", "TN", "UT")

# Aggregate all temperature data 
all_temp <- do.call(rbind, lapply(Location, function(st) {
  f <- file.path("data", "temperature", paste0(st, "_temp.csv"))
  df <- read.csv(f)
  df$Location <- st
  return(df)
}))

state_temp = all_temp %>%
  mutate(date = as.Date(date)) %>%
  group_by(Location, date) %>%
  mutate(mean = mean(mean)) %>%
  distinct(Location, date, mean) %>%
  ungroup() %>%
  mutate(temp= mean) %>%
  drop_na() %>%
  group_by(Location) %>%
  mutate(inverse_temp = - temp) %>%
  mutate(inverse_temp = inverse_temp - min(inverse_temp))  %>%
 # mutate(inverse_temp = ifelse(inverse_temp < 15, 15, inverse_temp)) %>%
  mutate(temp_smooth = rollmean(as.numeric(temp ), k = 30, fill = NA)) %>%
  mutate(temp_inverse_smooth = rollmean(as.numeric(inverse_temp ), k = 30, fill = NA)) %>%
  ungroup()

ggplot(data = state_temp) +
  geom_line(aes(x = date, y = temp_inverse_smooth, col = Location), lwd = 1.3) + 
  theme_bw() +
  scale_color_viridis_d() +
  ylab("Mean temperature in Celsius (degrees)") + 
  xlab("Date") + facet_wrap(vars(Location))


# Humidity data  ###########################################################

# states to read in
Location = c("CA", "CO", "CT", "GA", "MD",
             "MI", "MN", "NM", "NY", "OH", 
             "OR", "TN", "UT")

# Aggregate all temperature data 
all_hum <- do.call(rbind, lapply(Location, function(st) {
  f <- file.path("data", "humidity", paste0(st, "_hum.csv"))
  df <- read.csv(f)
  df$Location <- st
  return(df)
}))

state_hum = all_hum %>%
  group_by(Location, date) %>%
  mutate(mean = mean(mean)) %>%
  distinct(date, mean) %>%
  ungroup() %>%
  mutate(humid = mean) %>%
  group_by(Location) %>%
  mutate(humid_smooth = rollmean(as.numeric(humid), k = 45, fill = NA))  %>%
  mutate(date = as.Date(date)) %>%
  ungroup()
head(state_hum)

ggplot(data = state_hum) +
  geom_line(aes(x = date, y = humid_smooth, col = Location), lwd = 1.3) + 
  theme_bw() +
  scale_color_viridis_d() +
  ylab("Relative humidity (%)") + 
  xlab("Date") # + facet_wrap(vars(Location))

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
  geom_line() +
  geom_line(aes(x = week, y = rescaled_share))

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
  mutate(week= as.Date(week))

# Visualize 
ggplot(data = frequency, aes(x = week, y =rescaled_lof, col = usa_or_hhsregion )) +
  geom_line(lwd = 1.5) +
  theme_bw() +
  ylab("Frequency odds") + xlab("Date") +
  scale_color_viridis_d() +
  ggtitle('CDC variant date')


##################################################################
# Transform into weekly data and make into one large calibration dataset 

# Hosp #######################################################
all_states_hosp <- net %>%
  left_join(census_dat, by = 'Location') %>%
  drop_na() %>%
  mutate(hosp = Weekly.Rate*(pop_size_2023/100000)) %>%
  mutate(hosp = round(hosp, digits = 0)) %>%
  mutate(date = week) %>% dplyr::select(-week) %>%
  mutate(week = week(date), year = year(date)) %>%
  dplyr::select(hosp, week, year, Location)
head(all_states_hosp) 

ggplot(data = all_states_hosp, aes(x = week, y = hosp)) +
  geom_line() +
  facet_wrap(vars(Location), scales = "free_y") +
  ylab("Weekly hospitalizations")

# Temperature  #######################################################
weekly_temp = state_temp %>%
  mutate(week = week(date), year = year(date)) %>%
  group_by(year, week, Location) %>%
  dplyr::summarise(across(temp:inverse_temp, mean, na.rm = TRUE), .groups = "drop")
head(weekly_temp)

# Humidity #######################################################
weekly_hum = state_hum %>%
  mutate(week = week(date), year = year(date)) %>%
  group_by(year, week, Location) %>%
  dplyr::summarise(across(humid_smooth, mean, na.rm = TRUE), .groups = "drop")
head(weekly_hum)

# Variants #######################################################
Location = c("CA", "CO", "CT", "GA", "MD",
             "MI", "MN", "NM", "NY", "OH", 
             "OR", "TN", "UT")
usa_or_hhsregion = c(9, 8, 1, 4, 3, 5, 5, 
                    6, 2, 5, 10, 4, 8)
location_hhs = data.frame(Location, usa_or_hhsregion) %>%
  mutate(usa_or_hhsregion = as.character(usa_or_hhsregion))

weekly_variants = frequency %>%
  mutate(date = week) %>%
  dplyr::select(-week) %>%
  mutate(week = week(date), year = year(date) ) %>%
  left_join(location_hhs, by = "usa_or_hhsregion", relationship = "many-to-many") %>%
  dplyr::select(Location, week, year, odds_frequency)
head(weekly_variants)

# Vaccination #######################################################
# Up to 2023
head(covid_vax_rate)

#weekly_vax_rate_2223 = covid_vax_rate %>%
#  mutate(week = week(date)) %>%
#  mutate(
#    month_start = as.Date(paste(year, month, "01", sep = "-"))) 
#head(weekly_vax_rate_2223)

# 2. Generate all weekly dates per location
#weekly_grid <- weekly_vax_rate_2223 %>%
#  group_by(Location) %>%
#  dplyr::summarise(
#    week_start = seq.Date(
#      from = min(month_start),
#      to   = max(month_start) + days(days_in_month(max(month_start)) - 1),
#      by   = "week"),
#    .groups = "drop")
#head(weekly_grid)

#weekly_vax_2223 <- weekly_grid %>%
#  group_by(Location) %>%
#  do({
#    df_loc <- filter(weekly_vax_rate_2223, Location == .$Location[1])
#    data.frame(
#      week_start = .$week_start,
#      incident   = approx(x = df_loc$month_start, y = df_loc$month_vax_rate, xout = .$week_start, rule = 2)$y
#    )
#  }) %>%
#  ungroup() %>%
#  mutate(week = week(week_start), year = year(week_start)) %>%
#  mutate(incident = incident/100) %>%
#  dplyr::select(-week_start)
#head(weekly_vax_2223)

# 2023-2024 
head(covid_vax_rate)
weekly_vax_2223 = covid_vax_rate %>%
  dplyr::select(week, year, smooth_vax_rate, Location) %>%
  drop_na() %>%
  mutate(incident = smooth_vax_rate) %>%
  dplyr::select(-smooth_vax_rate)
head(weekly_vax_2223)


# 2024-2025
# Join with common vaccination data 
Location = c("CA", "CO", "CT", "GA", "MD",
             "MI", "MN", "NM", "NY", "OH", 
             "OR", "TN", "UT")
Jurisdiction = c("California", "Colorado", "Connecticut", 
                 "Georgia", "Maryland", "Michigan", "Minnesota", 
                 "New Mexico", "New York State", "Ohio", "Oregon",
                 "Tennessee", "Utah")
location_jur = data.frame(Location, Jurisdiction)

monthly_vax <- covid_vax_2425 %>%
  mutate(
    month_start = as.Date(paste(year, month_num, "01", sep = "-")))
head(monthly_vax)

monthly_vax <- monthly_vax %>%
  left_join(location_jur, by = "Jurisdiction") %>%
  mutate(month_start = as.Date(paste(year, month_num, "01", sep = "-"))) %>%
  arrange(Location, month_start)

weekly_grid <- monthly_vax %>%
  group_by(Location) %>%
  reframe(
    week_start = seq.Date(
      from = min(month_start),
      to = max(month_start) + days(days_in_month(max(month_start)) - 1),
      by = "week"
    )
  )

#  smooth interpolation, then adjust to preserve totals
weekly_vax_2425 <- weekly_grid %>%
  group_by(Location) %>%
  arrange(week_start) %>%
  do({
    location_data <- filter(monthly_vax, Location == .$Location[1])
    week_dates <- .$week_start
    
    # Step 1: Create smooth interpolation using spline
    smooth_values <- spline(
      x = as.numeric(location_data$month_start),
      y = location_data$incident,
      xout = as.numeric(week_dates),
      method = "natural"
    )$y
    
    # Step 2: Adjust to preserve monthly totals
    weekly_data <- data.frame(
      week_start = week_dates,
      smooth_incident = smooth_values,
      month_start = floor_date(week_dates, "month")
    )
    
    # Calculate adjustment factors for each month
    monthly_adjustments <- weekly_data %>%
      group_by(month_start) %>%
      dplyr::summarise(
        smooth_sum = sum(smooth_incident),
        .groups = "drop"
      ) %>%
      left_join(
        location_data %>% dplyr::select(month_start, incident),
        by = "month_start"
      ) %>%
      mutate(
        adjustment_factor = incident / smooth_sum
      )
    
    # Apply adjustments to preserve totals
    final_data <- weekly_data %>%
      left_join(monthly_adjustments %>% dplyr::select(month_start, adjustment_factor),
                by = "month_start") %>%
      mutate(
        incident = smooth_incident * adjustment_factor
      ) %>%
      dplyr::select(week_start, incident)
    
    final_data
  }) %>%
  ungroup() %>%
  mutate(
    week = week(week_start),
    year = year(week_start)
  ) %>%
  dplyr::select(Location, incident, week, year) %>%
  filter(!is.na(incident))
head(weekly_vax_2425)

check_again = weekly_vax_2425 %>%
  group_by(year, Location) %>%
  mutate(annual = sum(incident)) %>%
  distinct(year, Location, annual)

ggplot(data = check_again, aes(x = as.character(year), y = annual, col = Location)) +
  geom_point() + xlab("Year") + theme_bw()

# Bind together 
weekly_vaccination_dat = rbind(weekly_vax_2223, weekly_vax_2425 ) %>%
  ungroup()
head(weekly_vaccination_dat)


##############################################################
# Join all covariates and filter by time 

full_calibration_dat = left_join(all_states_hosp, weekly_temp, by = c("Location", "year", "week") ) %>%
  left_join(weekly_hum, by = c("Location", "year", "week")) %>%
  left_join(weekly_variants,  by = c("Location", "year", "week")) %>%
  left_join(weekly_vaccination_dat, by = c("Location", "year", "week")) %>%
  mutate(time = week/52 + year) %>%
  filter(time > 2022 & time < 2024.75) %>%
  arrange(Location, time) %>%
  group_by(Location) %>%
  mutate(incident = na.approx(incident, x = time, na.rm = FALSE) ) %>%
  mutate(humid_smooth = ifelse(is.na(humid_smooth), 78, humid_smooth)) %>%
  ungroup() %>%
  filter(Location != "MD") %>%
  group_by(Location) %>%
  mutate( t = seq(1, 143, 1)) %>%
  ungroup() 
  
write.csv(full_calibration_dat,  "processed_calibration_dat_humidsmooth.csv")
getwd()

# 2 week lag 
full_calibration_dat = left_join(all_states_hosp, weekly_temp, by = c("Location", "year", "week") ) %>%
  left_join(weekly_hum, by = c("Location", "year", "week")) %>%
  left_join(weekly_variants,  by = c("Location", "year", "week")) %>%
  left_join(weekly_vaccination_dat, by = c("Location", "year", "week")) %>%
  mutate(time = week/52 + year) %>%
  filter(time > 2022 & time < 2024.75) %>%
  arrange(Location, time) %>%
  group_by(Location) %>%
  mutate(incident = na.approx(incident, x = time, na.rm = FALSE) ) %>%
  mutate(humid_smooth = ifelse(is.na(humid_smooth), 78, humid_smooth)) %>%
  ungroup() %>%
  filter(Location != "MD") %>%
  group_by(Location) %>%
  ungroup() %>%
  group_by(Location) %>%
  mutate(inverse_temp = lag(inverse_temp, 2)) %>%
  mutate(humid_smooth = lag(humid_smooth, 2)) %>%
  drop_na() %>%
  mutate( t = seq(1, 141, 1)) %>% ungroup()
head(full_calibration_dat) 

write.csv(full_calibration_dat,  "processed_calibration_dat_lag2weeks.csv")

# Plot out all covariates to assess 
# hosp
ggplot(data = full_calibration_dat, aes( x = time, y = hosp)) +
  geom_line(aes(col = Location), lwd = 1.3 ) + facet_wrap(vars(Location), scales = "free_y") +
  ylab("Weekly hospitalizations") + theme_bw()  + scale_color_viridis_d()
# humid
ggplot(data = full_calibration_dat, aes( x = time, y = humid_smooth)) +
  geom_line(aes(col = Location), lwd = 1.3 ) + facet_wrap(vars(Location)) +
  ylab("Weekly humidity") + theme_bw()  + scale_color_viridis_d()
# temp
ggplot(data = full_calibration_dat, aes( x = time, y = inverse_temp)) +
  geom_line(aes(col = Location), lwd = 1.3) + facet_wrap(vars(Location)) +
  ylab("Weekly temp")+ theme_bw()  + scale_color_viridis_d()
# incident 
ggplot(data = full_calibration_dat, aes( x = time, y = incident)) +
  geom_line(aes(col = Location), lwd = 1.3) + facet_wrap(vars(Location)) + theme_bw() +
  ylab("Weekly vaccination (%)") + theme_bw()  + scale_color_viridis_d()

# check vax again 
check_vax_fin = full_calibration_dat %>%
  group_by(year, Location) %>%
  mutate(sum_vax = sum(incident)) %>%
  distinct(year, Location, sum_vax)

# variants
ggplot(data = full_calibration_dat, aes( x = time, y = odds_frequency)) +
  geom_line(aes(col = Location), lwd = 1.3) + facet_wrap(vars(Location))+
  ylab("Weekly variant odds") + theme_bw()  + scale_color_viridis_d()

head(full_calibration_dat)

##############################################################
##############################################################
# Model set up 
seirs_model_hybrid = function(t, state, param) {
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
    vax_rate_at_t = calibration_dataset$incident[t] 
    humidity_at_t = calibration_dataset$humid_smooth[t] 
    temp_at_t = calibration_dataset$inverse_temp[t]
    variant_at_t = calibration_dataset$odds_frequency[t] 
    
    
    # Seasonal function 
    seasonal = beta*( var*variant_at_t + hum*humidity_at_t + temp*temp_at_t)
    
    dS = -seasonal*S*(I + I2)/N - VE*vax_rate_at_t*S + imm_vax*V
    dE =  seasonal*S*(I+ I2)/N - lp*E
    dI =  lp*E - gamma*I
    dR = I*gamma -  R*imm   
    dV = VE*vax_rate_at_t*S  - imm_vax*V   
    
    dS2 = -seasonal*S2*(I + I2)/N  + imm*R  + imm*R2 + imm_vax*V2 - VE*vax_rate_at_t*S2
    dE2 =  seasonal*S2*(I + I2)/N - lp*E2
    dI2 =  lp*E2 - gamma*I2
    dR2 = I2*gamma -  R2*imm
    dV2 = VE*vax_rate_at_t*S2  - imm_vax*V2  
    
    # Return the rates of change
    list(c(dS, dE, dI, dR, dV, dS2, dE2, dI2, dR2, dV2))
  })
}


##############################################################
##############################################################
# MCMC fitting function 

run_hierarchical_mcmc_calibration <- function(full_calibration_dat, iterations = 2000, burn_in = 300, thin = 2, n_chains = 3) {
  
  # states to be fit in the calibration
  states <- unique(full_calibration_dat$Location)
  n_states <- length(states)
  cat("Fitting", n_states, "states:", paste(states, collapse = ", "), "\n")
  
  # Define the model for a single state
  seirs_model_hybrid <- function(t, state, param, state_data) {
    
    # Extract state variables
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
    
    # Global parameters
    var <- param[["var"]]
    hum <- param[["hum"]] 
    temp <- param[["temp"]]
    lp <- param[["lp"]]
    gamma <- param[["gamma"]]
    imm <- param[["imm"]]
    imm_vax <- param[["imm_vax"]]
    VE <- param[["VE"]]
    hr1 <- param[["hr1"]]
    hr2 <- param[["hr2"]]
    
    # State-specific parameters
    beta <- param[["beta"]]
    N <- param[["N"]]
    
    # Get time series data for this time point
    t_idx <- round(t)  # Ensure integer indexing
    if(t_idx >= 1 && t_idx <= nrow(state_data)) {
      vax_rate_at_t <- state_data$incident[t_idx]
      humidity_at_t <- state_data$humid[t_idx]
      temp_at_t <- state_data$inverse_temp[t_idx]
      variant_at_t <- state_data$odds_frequency[t_idx]
    } else {
      # Use last available values if t exceeds data
      vax_rate_at_t <- tail(state_data$incident, 1)
      humidity_at_t <- tail(state_data$humid, 1)
      temp_at_t <- tail(state_data$inverse_temp, 1)
      variant_at_t <- tail(state_data$odds_frequency, 1)
    }
    

    if(is.na(vax_rate_at_t)) vax_rate_at_t <- 0
    if(is.na(humidity_at_t)) humidity_at_t <- 0.5
    if(is.na(temp_at_t)) temp_at_t <- 0.05
    if(is.na(variant_at_t)) variant_at_t <- 0.5
    
    # Seasonal function
    seasonal <- beta * (var * variant_at_t + hum * humidity_at_t + temp * temp_at_t)
    
    # Ensure seasonal is positive and finite
    if(!is.finite(seasonal) || seasonal < 0) seasonal <- 0.01
    
    # Differential equations
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
  
  # Hierarchical likelihood function
  likelihood_func <- function(params) {
    
    # Fixed parameters
    fixed_params <- c(
      lp = 1/4,        # latent period 
      gamma = 1/4,     # recovery rate
      VE = 0.30        # vaccine efficacy
    )
    
    # global parameters to be fit (shared across states)
    global_params <- c(
      var = params["var"],          # variant multiplier
      hum = params["hum"],          # humidity multiplier
      temp = params["temp"],        # temp multiplier
      imm = params["imm"],          # immunity
      imm_vax = params["imm_vax"],  # vaccine immunity
      hr1 = params["hr1"],          # hosp rate 1 
      hr2 = params["hr2"]           # hosp rate 2 
    )
    
    total_log_lik <- 0
    
    # Loop through each state
    for(i in 1:n_states) {
      
      # pull out state and subset to only that state's data 
      state_name <- states[i]
      state_data <- full_calibration_dat[full_calibration_dat$Location == state_name, ]
      state_data <- state_data[order(state_data$t), ]  # Ensure proper time ordering
      
    #  print(head(state_data))
      
      # State-specific parameters
      beta_param <- paste0("beta_", state_name)
      S0_param <- paste0("S0_", state_name)
      
      state_params <- c(
        beta = params[beta_param],
        N = state_populations[[state_name]]  # Use fixed population from list
      )
      
      # Combine all parameters
      all_params <- c(global_params, state_params, fixed_params)
      
      # Initial conditions for this state
      N_val <- state_params["N"]  # This comes from fixed population list
      S0_frac <- params[S0_param]
      
      S0 <- S0_frac * N_val
      E0 <- 0.01 * N_val
      I0 <- 0.001 * N_val
    #  R0 <- 0.30 * N_val  # Adjusted based on S0
      R0 = (1 - S0_frac - 0.01 - 0.001 - 0.049) * N_val
      V0 <- 0.049 * N_val
      
      S20 <- 0.00 * N_val
      E20 <- 0.00 * N_val
      I20 <- 0.00 * N_val
      R20 <- 0.00 * N_val
      V20 <- 0.00 * N_val
      
      state_init <- c(S = S0, E = E0, I = I0, R = R0, V = V0, 
                      S2 = S20, E2 = E20, I2 = I20, R2 = R20, V2 = V20)
      
      # Time vector for this state
      time_vec <- 1:nrow(state_data)
      
      
      # Solve ODE for this state
      tryCatch({
        out <- ode(y = state_init, times = time_vec, func = seirs_model_hybrid, 
                   parms = all_params, state_data = state_data)
        
        out_df <- as.data.frame(out)
        
        # ---- DEBUGGING PRINTS: check alignment before likelihood ----
        cat(sprintf("State %s: start_idx=%d, end_idx=%d, length(predicted)=%d, length(observed)=%d\n",
                    state_name, start_idx, end_idx, length(predicted_hosp), length(state_data$hosp)))
        cat("Observed weeks head:", paste(head(state_data$week), collapse=", "), "\n")
        cat("ODE time_vec head:", paste(head(time_vec), collapse=", "), "\n")
        cat("ODE time_vec tail:", paste(tail(time_vec), collapse=", "), "\n\n")
        
        # Calculate predicted hospitalizations
        predicted_hosp <- (out_df$I * all_params["hr1"] + out_df$I2 * all_params["hr2"])
        
        # Observed hospitalizations for this state
        observed_hosp <- state_data$hosp
        
        
        # Skip first few time points to avoid initialization effects
        start_idx <- min(10, length(observed_hosp) - 50)
        end_idx <- length(observed_hosp)
        
        # check
        
        if(start_idx < end_idx && length(predicted_hosp) >= end_idx) {
          
          # Calculate log-likelihood for this state using poisson likelihood
  
          state_log_lik <- sum(dpois(observed_hosp[start_idx:end_idx], 
                                     lambda = pmax(predicted_hosp[start_idx:end_idx], 1e-6), 
                                     log = TRUE), na.rm = TRUE)

          cat(sprintf("State %s: logLik = %.2f\n", state_name, state_log_lik))
          
          # sum of log-likelihoods across all states 
          total_log_lik <- total_log_lik + state_log_lik
          
        } else {
          return(-1e6)  # Penalty for invalid dimensions
        }
        
      }, error = function(e) {
        return(-1e6)  # Penalty for ODE solver failure
      })
    }
    
    return(total_log_lik)
  }
  
  # Hierarchical prior function
  # Set bounds on fitted parameters 
  log_prior <- function(params) {
    
    # Global parameter bounds
    global_bounds <- list(
      var = c(0.01, 2.0),
      hum = c(0.01, 2.0), 
      temp = c(0.01, 1.5),
      imm = c(1/300, 1/90),
      imm_vax = c(1/120, 1/30),
      hr1 = c(0.0001, 0.2),
      hr2 = c(0.0001, 0.2)
    )
    
    # State-specific parameter bounds
    state_bounds <- list(
      beta = c(0.05, 1.50),
      S0 = c(0.30, 0.80)
    )
    
    # Population sizes (you'll need to adjust these based on your data)
    state_populations <- list(
      CA = 39198693, CO = 5901339, CT = 3643023, GA = 11064432,
      MI = 10083356, MN = 5753048, NM = 2121164, NY = 19737367,
      OH = 11824034, OR = 4253653, TN = 7148304, VA = 3443222
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
      
      # Check S0
      S0_param <- paste0("S0_", state_name)
      if(S0_param %in% names(params)) {
        bounds <- state_bounds[["S0"]]
        if(params[S0_param] < bounds[1] || params[S0_param] > bounds[2]) {
          return(-Inf)
        }
      }
      
      # Population parameters are fixed (not checked here)
    }
    
    return(0)  # Uniform priors within bounds
  }
  
  # Posterior function
  # posterior = prior + likelihood 
  log_posterior <- function(params) {
    lp <- log_prior(params)
    
    if(is.infinite(lp)) {
      return(-Inf)
    }
    
    ll <- likelihood_func(params)
    return(lp + ll)
  }
  
  # Initial parameter values
  global_initial <- list(
    var = 0.25,
    hum = 0.30,
    temp = 0.15,
    imm = 1/180,
    imm_vax = 1/75,
    hr1 = 0.008,
    hr2 = 0.003
  )
  
  # State-specific initial values (only fitted parameters)
  state_initial <- list()
  for(state_name in states) {
    state_initial[[paste0("beta_", state_name)]] <- runif(1, 0.02, 1.25)
    state_initial[[paste0("S0_", state_name)]] <- runif(1, 0.45, 0.65)
  }
  
  # Combine initial parameters (only the ones we're fitting)
  initial_params <- c(global_initial, state_initial)
  param_names <- names(initial_params)
  
  # Fixed population sizes (not fitted)
  state_populations <- list(
    CA = 39198693, CO = 5901339, CT = 3643023, GA = 11064432,
    MI = 10083356, MN = 5753048, NM = 2121164, NY = 19737367,
    OH = 11824034, OR = 4253653, TN = 7148304, VA = 3443222
  )
  
  cat("Total parameters to fit:", length(param_names), "\n")
  cat("Global parameters:", length(global_initial), "\n")
  cat("State-specific fitted parameters per state:", (length(initial_params) - length(global_initial)) / n_states, "\n")
  cat("Population parameters (fixed):", length(states), "\n")
  
  # Rest of the MCMC code remains similar but with updated parameter handling...
  # [Include the rest of the MCMC chain running code from your original function]
  
  # Function to run a single MCMC chain
  # component wise Metropolis MCMC w adaptive proposal variance 
  # Accepts or rejects based on log posterior ratio 
  run_chain <- function(chain_id, iterations, burn_in, thin) {
    
    current_params <- unlist(initial_params)
    names(current_params) <- param_names
    
    # Perturb initial values for different chains (only fitted parameters)
    if(chain_id > 1) {
      current_params <- current_params * runif(length(current_params), 0.90, 1.10)
    }
    
    samples <- matrix(NA, nrow = ceiling((iterations - burn_in) / thin), 
                      ncol = length(param_names))
    colnames(samples) <- param_names
    
    # Check initial log posterior
    current_log_post <- log_posterior(current_params)
    cat(sprintf("Chain %d: Initial log posterior = %.2f\n", chain_id, current_log_post))
    
    if(!is.finite(current_log_post)) {
      cat(sprintf("Warning: Chain %d has non-finite initial log posterior. Adjusting initial values.\n", chain_id))
      # Try to find better initial values
      attempts <- 0
      while(!is.finite(current_log_post) && attempts < 10) {
        fitted_params <- !grepl("^N_", names(current_params))
        current_params[fitted_params] <- current_params[fitted_params] * runif(sum(fitted_params), 0.8, 1.2)
        current_log_post <- log_posterior(current_params)
        attempts <- attempts + 1
      }
      
      if(!is.finite(current_log_post)) {
        stop(sprintf("Chain %d: Cannot find valid initial parameters after 10 attempts", chain_id))
      }
      cat(sprintf("Chain %d: Found valid initial log posterior = %.2f after %d attempts\n", 
                  chain_id, current_log_post, attempts))
    }
    
    # Adaptive proposal standard deviations
    proposal_sd <- abs(current_params) * 0.10
    
    # Ensure minimum step sizes
    proposal_sd[proposal_sd < 1e-4] <- 1e-3 
    
    # Adaptation parameters
    adapt_start <- 100
    adapt_interval <- 50
    accept_count <- numeric(length(param_names))
    proposal_attempts <- numeric(length(param_names))
    target_accept_rate <- 0.234
    
    sample_idx <- 1
    
    for(i in 1:iterations) {
      
      # Component-wise updates
      for(j in 1:length(param_names)) {
        param_name <- param_names[j]
        
        # Skip if proposal_sd is zero or very small
        if(proposal_sd[j] <= 1e-10) next
        
        proposed_params <- current_params
        proposed_params[j] <- current_params[j] + rnorm(1, 0, proposal_sd[j])
        
        proposed_log_post <- log_posterior(proposed_params)
        
        if (i %% 100 == 0 && j == 1) {  # every 100 iterations, only once per sweep
          cat(sprintf("Iter %d: Current log post = %.2f\n", i, current_log_post))
        }
        
        if(!is.na(proposed_log_post) && is.finite(proposed_log_post)) {
          log_accept_ratio <- proposed_log_post - current_log_post
          proposal_attempts[j] <- proposal_attempts[j] + 1
          
          # Ensure log_accept_ratio is finite
          if(is.finite(log_accept_ratio) && log(runif(1)) < log_accept_ratio) {
            current_params <- proposed_params
            current_log_post <- proposed_log_post
            accept_count[j] <- accept_count[j] + 1
          }
        }
        
        # Adapt proposal standard deviation
        if(i > adapt_start && i %% adapt_interval == 0 && i <= burn_in) {
          if(proposal_attempts[j] > 0) {
            accept_rate_j <- accept_count[j] / proposal_attempts[j]
            
            if(accept_rate_j < target_accept_rate) {
              proposal_sd[j] <- max(proposal_sd[j] * 0.9, 1e-8)
            } else {
              proposal_sd[j] <- min(proposal_sd[j] * 1.1, 1)
            }
          }
        }
      }
      
      # Store samples after burn-in and thinning
      if(i > burn_in && (i - burn_in) %% thin == 0) {
        samples[sample_idx, ] <- current_params
        sample_idx <- sample_idx + 1
      }
      
      # Print progress
      if(i %% 500 == 0) {
        overall_accept <- sum(accept_count[proposal_attempts > 0]) / sum(proposal_attempts[proposal_attempts > 0])
        cat(sprintf("Chain %d, Iteration %d: Overall acceptance rate = %.3f\n", 
                    chain_id, i, overall_accept))
      }
    }
    
    return(list(
      samples = samples,
      acceptance_rate = sum(accept_count) / sum(proposal_attempts),
      param_accept_rates = accept_count / pmax(proposal_attempts, 1)
    ))
  }
  
  # Run multiple chains
  cat("Starting hierarchical MCMC with", n_chains, "chains,", iterations, "iterations\n")
  
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
  
  # Combine results
  all_samples <- list()
  acceptance_rates <- numeric(n_chains)
  param_accept_rates <- matrix(0, nrow = n_chains, ncol = length(param_names))
  colnames(param_accept_rates) <- param_names
  
  for(i in 1:n_chains) {
    all_samples[[i]] <- chain_results[[i]]$samples
    acceptance_rates[i] <- chain_results[[i]]$acceptance_rate
    param_accept_rates[i, ] <- chain_results[[i]]$param_accept_rates
  }
  
  # Create MCMC objects
  mcmc_objects <- lapply(all_samples, mcmc)
  mcmc_list <- mcmc.list(mcmc_objects)
  
  # Diagnostics - with error handling
  gelman_diag <- NULL
  if(n_chains > 1) {
    tryCatch({
      # Check if chains have converged enough for Gelman diagnostic
      # Remove any columns with zero variance or infinite values
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
        cat("This often indicates poor chain mixing or convergence issues\n")
      }
    }, error = function(e) {
      cat("Warning: Could not calculate Gelman-Rubin diagnostics:", e$message, "\n")
      cat("This is often due to poor chain convergence or numerical instability\n")
    })
  }
  
  summary_stats <- summary(mcmc_list)
  
  return(list(
    samples = all_samples,
    mcmc_list = mcmc_list,
    acceptance_rates = acceptance_rates,
    param_accept_rates = param_accept_rates,
    gelman_diag = gelman_diag,
    summary = summary_stats,
    states = states,
    param_names = param_names
  ))
}


#############################################################
### Run the hierarchical MCMC calibration 

results <- run_hierarchical_mcmc_calibration(
  full_calibration_dat = full_calibration_dat,  # your data
  iterations = 1000,     # start small for testing
  burn_in = 100,         # burn-in period
  thin = 1,              # thinning interval
  n_chains = 1           # number of MCMC chains
)



#############################################################
### Post-processing MCMC results 

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

# Extract specific parameter types
global_params <- results$summary$statistics[
  rownames(results$summary$statistics) %in% c("var", "hum", "temp", "imm", "imm_vax", "hr1", "hr2"), 
]
print("Global parameter estimates:")
print(global_params)

# Extract state-specific beta parameters
beta_params <- results$summary$statistics[
  grepl("beta_", rownames(results$summary$statistics)), 
]
print("State-specific beta estimates:")
print(beta_params)

# Extract state-specific S0 parameters
S0_params <- results$summary$statistics[
  grepl("S0_", rownames(results$summary$statistics)), 
]
print("State-specific S0 estimates:")
print(S0_params)

# Plot trace plots for key parameters (optional)
# You can plot traces to check convergence visually
if(requireNamespace("coda", quietly = TRUE)) {
  # Plot traces for global parameters
  plot(results$mcmc_list[, "var"], main = "Trace plot for 'var' parameter")
  plot(results$mcmc_list[, "imm"], main = "Trace plot for 'imm' parameter")
  
  # Plot trace for a state-specific parameter
  first_state <- results$states[2]
  beta_param_name <- paste0("beta_", first_state)
  if(beta_param_name %in% colnames(results$mcmc_list[[1]])) {
    plot(results$mcmc_list[, beta_param_name], 
         main = paste("Trace plot for", beta_param_name))
  }
}





###########################################
# Aug 28 try to debug 

run_hierarchical_mcmc_calibration <- function(full_calibration_dat, iterations = 2000, burn_in = 300, thin = 2, n_chains = 3) {
  
  # Load required library
  if (!requireNamespace("deSolve", quietly = TRUE)) {
    stop("deSolve package is required but not installed. Please install it with: install.packages('deSolve')")
  }
  library(deSolve)
  
  # states to be fit in the calibration
  states <- unique(full_calibration_dat$Location)
  n_states <- length(states)
  cat("Fitting", n_states, "states:", paste(states, collapse = ", "), "\n")
  
  # Define state_populations at the top level so it's available throughout
  state_populations <- list(
    "CA" = 39198693, "CO" = 5901339, "CT" = 3643023, "GA" = 11064432,
    "MI" = 10083356, "MN" = 5753048, "NM" = 2121164, "NY" = 19737367,
    "OH" = 11824034, "OR" = 4253653, "TN" = 7148304, "UT" = 3380800
  )
  
  # Define the model for a single state - CLEAN VERSION
  seirs_model_hybrid <- function(t, state, parms) {
    
    # Extract state variables
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
    
    # Extract parameters
    var <- parms[["var"]]
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
    
    # Get the state_data from the parameter list
    state_data <- parms[["state_data"]]
    
    # Get time series data for this time point
    t_idx <- round(t)
    if(t_idx < 1) t_idx <- 1
    if(t_idx > nrow(state_data)) t_idx <- nrow(state_data)
    
    # Access data with bounds checking
    vax_rate_at_t <- state_data$incident[t_idx]
    humidity_at_t <- state_data$humid[t_idx]
    temp_at_t <- state_data$inverse_temp[t_idx]
    variant_at_t <- state_data$odds_frequency[t_idx]
    
    # Handle NA values
    if(is.na(vax_rate_at_t)) vax_rate_at_t <- 0
    if(is.na(humidity_at_t)) humidity_at_t <- 0.5
    if(is.na(temp_at_t)) temp_at_t <- 0.05
    if(is.na(variant_at_t)) variant_at_t <- 0.5
    
    # Seasonal function
    seasonal <- beta * (var * variant_at_t + hum * humidity_at_t + temp * temp_at_t)
    
    # Ensure seasonal is positive and finite
    if(!is.finite(seasonal) || seasonal < 0) seasonal <- 0.01
    
    # Differential equations
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
  
  # Hierarchical likelihood function
  likelihood_func <- function(params) {
    
    # Fixed parameters
    fixed_params <- c(
      lp = 1/4,        # latent period 
      gamma = 1/4,     # recovery rate
      VE = 0.30        # vaccine efficacy
    )
    
    # global parameters to be fit (shared across states)
    global_params <- c(
      var = params["var"],          # variant multiplier
      hum = params["hum"],          # humidity multiplier
      temp = params["temp"],        # temp multiplier
      imm = params["imm"],          # immunity
      imm_vax = params["imm_vax"],  # vaccine immunity
      hr1 = params["hr1"],          # hosp rate 1 
      hr2 = params["hr2"]           # hosp rate 2 
    )
    
    total_log_lik <- 0
    
    # Loop through each state
    for(i in 1:n_states) {
      
      # pull out state and subset to only that state's data 
      state_name <- states[i]
      state_data <- full_calibration_dat[full_calibration_dat$Location == state_name, ]
      state_data <- state_data[order(state_data$t), ]  # Ensure proper time ordering
      
      # Check if state_data is empty
      if(nrow(state_data) == 0) {
        cat("Warning: No data for state", state_name, "\n")
        return(-1e6)
      }
      
      # Debug: Check data structure
      cat(sprintf("State %s: nrows=%d, ncols=%d\n", state_name, nrow(state_data), ncol(state_data)))
      
      # Check for required columns
      required_cols <- c("incident", "humid", "inverse_temp", "odds_frequency", "hosp")
      missing_cols <- setdiff(required_cols, colnames(state_data))
      if(length(missing_cols) > 0) {
        cat(sprintf("Warning: Missing columns for state %s: %s\n", state_name, paste(missing_cols, collapse=", ")))
        return(-1e6)
      }
      
      # State-specific parameters: beta is fitted, N is fixed
      beta_param <- paste0("beta_", state_name)
      S0_param <- paste0("S0_", state_name)
      
      # Check if state exists in population list
      if(!state_name %in% names(state_populations)) {
        cat("Warning: Population not found for state", state_name, "\n")
        return(-1e6)
      }
      
      # Check if fitted parameters exist in params vector
      if(!beta_param %in% names(params)) {
        cat(sprintf("Error: Parameter '%s' not found in params vector\n", beta_param))
        return(-1e6)
      }
      
      if(!S0_param %in% names(params)) {
        cat(sprintf("Error: Parameter '%s' not found in params vector\n", S0_param))
        return(-1e6)
      }
      
      # Access fitted parameters
      beta_value <- params[[beta_param]]
      S0_value <- params[[S0_param]]
      
      # Check for missing fitted parameters
      if(is.na(beta_value) || is.na(S0_value)) {
        cat("Warning: Missing fitted parameters for state", state_name, "\n")
        return(-1e6)
      }
      
      # Create parameter list for ODE: fitted + fixed + state_data
      all_params <- list(
        # Global fitted parameters
        var = global_params["var"],
        hum = global_params["hum"],
        temp = global_params["temp"],
        imm = global_params["imm"],
        imm_vax = global_params["imm_vax"],
        hr1 = global_params["hr1"],
        hr2 = global_params["hr2"],
        # State-specific fitted parameters
        beta = beta_value,
        # Fixed parameters
        N = state_populations[[state_name]],
        lp = fixed_params["lp"],
        gamma = fixed_params["gamma"],
        VE = fixed_params["VE"],
        # Data for this state
        state_data = state_data
      )
      
      # Initial conditions for this state
      N_val <- state_populations[[state_name]]
      S0_frac <- S0_value
      
      S0 <- S0_frac * N_val
      E0 <- 0.01 * N_val
      I0 <- 0.001 * N_val
      R0 <- (1 - S0_frac - 0.01 - 0.001 - 0.049) * N_val
      V0 <- 0.049 * N_val
      
      S20 <- 0.00 * N_val
      E20 <- 0.00 * N_val
      I20 <- 0.00 * N_val
      R20 <- 0.00 * N_val
      V20 <- 0.00 * N_val
      
      state_init <- c(S = S0, E = E0, I = I0, R = R0, V = V0, 
                      S2 = S20, E2 = E20, I2 = I20, R2 = R20, V2 = V20)
      
      # Time vector for this state
      time_vec <- 1:nrow(state_data)
      
      # Solve ODE for this state - CLEAN CALL
      tryCatch({
        cat(sprintf("Solving ODE for state %s...\n", state_name))
        out <- ode(y = state_init, times = time_vec, func = seirs_model_hybrid, parms = all_params)
        cat(sprintf("ODE solved successfully for state %s\n", state_name))
        
        out_df <- as.data.frame(out)
        
        # Calculate predicted hospitalizations
        predicted_hosp <- (out_df$I * all_params[["hr1"]] + out_df$I2 * all_params[["hr2"]])
        
        # Observed hospitalizations for this state
        observed_hosp <- state_data$hosp
        
        # Skip first few time points to avoid initialization effects
        start_idx <- min(10, length(observed_hosp) - 50)
        end_idx <- length(observed_hosp)
        
        # Check for alignment issues
        if(start_idx < end_idx && length(predicted_hosp) >= end_idx && length(observed_hosp) >= end_idx) {
          
          # Extract the relevant portions
          pred_subset <- predicted_hosp[start_idx:end_idx]
          obs_subset <- observed_hosp[start_idx:end_idx]
          
          # Check for valid data
          if(any(is.na(pred_subset)) || any(is.na(obs_subset))) {
            cat("Warning: NA values found for state", state_name, "\n")
            return(-1e6)
          }
          
          # Ensure predictions are positive
          pred_subset <- pmax(pred_subset, 1e-6)
          
          # Calculate log-likelihood for this state using poisson likelihood
          state_log_lik <- sum(dpois(obs_subset, lambda = pred_subset, log = TRUE), na.rm = TRUE)
          
          cat(sprintf("State %s: logLik = %.2f (mean pred=%.2f, mean obs=%.2f)\n", 
                      state_name, state_log_lik, mean(pred_subset), mean(obs_subset)))
          
          # Check if log-likelihood is finite
          if(!is.finite(state_log_lik)) {
            cat("Warning: Non-finite log-likelihood for state", state_name, "\n")
            return(-1e6)
          }
          
          # sum of log-likelihoods across all states 
          total_log_lik <- total_log_lik + state_log_lik
          
        } else {
          cat(sprintf("Error: Dimension mismatch for state %s\n", state_name))
          return(-1e6)
        }
        
      }, error = function(e) {
        cat("ODE solver error for state", state_name, ":", e$message, "\n")
        return(-1e6)
      })
    }
    
    cat(sprintf("Total log-likelihood across all states: %.2f\n", total_log_lik))
    return(total_log_lik)
  }
  
  # Hierarchical prior function
  log_prior <- function(params) {
    
    # Global parameter bounds
    global_bounds <- list(
      var = c(0.01, 2.0),
      hum = c(0.01, 2.0), 
      temp = c(0.01, 1.5),
      imm = c(1/300, 1/90),
      imm_vax = c(1/120, 1/30),
      hr1 = c(0.0001, 0.2),
      hr2 = c(0.0001, 0.2)
    )
    
    # State-specific parameter bounds
    state_bounds <- list(
      beta = c(0.05, 1.50),
      S0 = c(0.30, 0.80)
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
      
      # Check S0
      S0_param <- paste0("S0_", state_name)
      if(S0_param %in% names(params)) {
        bounds <- state_bounds[["S0"]]
        if(params[S0_param] < bounds[1] || params[S0_param] > bounds[2]) {
          return(-Inf)
        }
      }
    }
    
    return(0)  # Uniform priors within bounds
  }
  
  # Posterior function
  log_posterior <- function(params) {
    lp <- log_prior(params)
    
    if(is.infinite(lp)) {
      return(-Inf)
    }
    
    ll <- likelihood_func(params)
    return(lp + ll)
  }
  
  # Initial parameter values
  global_initial <- list(
    var = 0.25,
    hum = 0.30,
    temp = 0.15,
    imm = 1/180,
    imm_vax = 1/75,
    hr1 = 0.008,
    hr2 = 0.003
  )
  
  # State-specific initial values (only fitted parameters)
  state_initial <- list()
  for(state_name in states) {
    state_initial[[paste0("beta_", state_name)]] <- runif(1, 0.2, 0.8)
    state_initial[[paste0("S0_", state_name)]] <- runif(1, 0.45, 0.65)
  }
  
  # Combine initial parameters (only the ones we're fitting)
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
  
  # Function to run a single MCMC chain
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
    
    # Check initial log posterior
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
    
    # Adaptive proposal standard deviations
    proposal_sd <- abs(current_params) * 0.10
    proposal_sd[proposal_sd < 1e-4] <- 1e-3 
    
    # Adaptation parameters
    adapt_start <- 100
    adapt_interval <- 50
    accept_count <- numeric(length(param_names))
    proposal_attempts <- numeric(length(param_names))
    target_accept_rate <- 0.234
    
    sample_idx <- 1
    
    for(i in 1:iterations) {
      
      # Component-wise updates
      for(j in 1:length(param_names)) {
        param_name <- param_names[j]
        
        if(proposal_sd[j] <= 1e-10) next
        
        proposed_params <- current_params
        proposed_params[j] <- current_params[j] + rnorm(1, 0, proposal_sd[j])
        
        proposed_log_post <- log_posterior(proposed_params)
        
        if (i %% 100 == 0 && j == 1) {
          cat(sprintf("Iter %d: Current log post = %.2f\n", i, current_log_post))
        }
        
        if(!is.na(proposed_log_post) && is.finite(proposed_log_post)) {
          log_accept_ratio <- proposed_log_post - current_log_post
          proposal_attempts[j] <- proposal_attempts[j] + 1
          
          if(is.finite(log_accept_ratio) && log(runif(1)) < log_accept_ratio) {
            current_params <- proposed_params
            current_log_post <- proposed_log_post
            accept_count[j] <- accept_count[j] + 1
          }
        }
        
        # Adapt proposal standard deviation
        if(i > adapt_start && i %% adapt_interval == 0 && i <= burn_in) {
          if(proposal_attempts[j] > 0) {
            accept_rate_j <- accept_count[j] / proposal_attempts[j]
            
            if(accept_rate_j < target_accept_rate) {
              proposal_sd[j] <- max(proposal_sd[j] * 0.9, 1e-8)
            } else {
              proposal_sd[j] <- min(proposal_sd[j] * 1.1, 1)
            }
          }
        }
      }
      
      # Store samples after burn-in and thinning
      if(i > burn_in && (i - burn_in) %% thin == 0) {
        samples[sample_idx, ] <- current_params
        sample_idx <- sample_idx + 1
      }
      
      # Print progress
      if(i %% 500 == 0) {
        overall_accept <- sum(accept_count[proposal_attempts > 0]) / sum(proposal_attempts[proposal_attempts > 0])
        cat(sprintf("Chain %d, Iteration %d: Overall acceptance rate = %.3f\n", 
                    chain_id, i, overall_accept))
      }
    }
    
    return(list(
      samples = samples,
      acceptance_rate = sum(accept_count) / sum(proposal_attempts),
      param_accept_rates = accept_count / pmax(proposal_attempts, 1)
    ))
  }
  
  # Run multiple chains
  cat("Starting hierarchical MCMC with", n_chains, "chains,", iterations, "iterations\n")
  
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
  
  # Combine results
  all_samples <- list()
  acceptance_rates <- numeric(n_chains)
  param_accept_rates <- matrix(0, nrow = n_chains, ncol = length(param_names))
  colnames(param_accept_rates) <- param_names
  
  for(i in 1:n_chains) {
    all_samples[[i]] <- chain_results[[i]]$samples
    acceptance_rates[i] <- chain_results[[i]]$acceptance_rate
    param_accept_rates[i, ] <- chain_results[[i]]$param_accept_rates
  }
  
  # Create MCMC objects
  mcmc_objects <- lapply(all_samples, mcmc)
  mcmc_list <- mcmc.list(mcmc_objects)
  
  # Diagnostics - with error handling
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
    param_accept_rates = param_accept_rates,
    gelman_diag = gelman_diag,
    summary = summary_stats,
    states = states,
    param_names = param_names
  ))
}


results <- run_hierarchical_mcmc_calibration(
  full_calibration_dat = full_calibration_dat,  # your data
  iterations = 100,     # start small for testing
  burn_in = 20,         # burn-in period
  thin = 1,              # thinning interval
  n_chains = 1           # number of MCMC chains
)









#####################################
# DEBUG AGAIN 
# This currently seems to be producing functioning likelihoods

run_hierarchical_mcmc_calibration <- function(full_calibration_dat, iterations = 2000, burn_in = 300, thin = 2, n_chains = 3) {
  
  # Load required library
  if (!requireNamespace("deSolve", quietly = TRUE)) {
    stop("deSolve package is required but not installed. Please install it with: install.packages('deSolve')")
  }
  library(deSolve)
  
  # states to be fit in the calibration
  states <- unique(full_calibration_dat$Location)
  n_states <- length(states)
  cat("Fitting", n_states, "states:", paste(states, collapse = ", "), "\n")
  
  # Define state_populations at the top level so it's available throughout
  state_populations <- list(
    "CA" = 39198693, "CO" = 5901339, "CT" = 3643023, "GA" = 11064432,
    "MI" = 10083356, "MN" = 5753048, "NM" = 2121164, "NY" = 19737367,
    "OH" = 11824034, "OR" = 4253653, "TN" = 7148304, "UT" = 3380800
  )
  
  # Define the model for a single state - CLEAN VERSION
  seirs_model_hybrid <- function(t, state, parms) {
    
    # Extract state variables
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
    
    # Extract parameters
    var <- parms[["var"]]
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
    
    # Get the state_data from the parameter list
    state_data <- parms[["state_data"]]
    
    # Get time series data for this time point
    t_idx <- round(t)
    if(t_idx < 1) t_idx <- 1
    if(t_idx > nrow(state_data)) t_idx <- nrow(state_data)
    
    # Access data with bounds checking
    vax_rate_at_t <- state_data$incident[t_idx]
    humidity_at_t <- state_data$humid[t_idx]
    temp_at_t <- state_data$inverse_temp[t_idx]
    variant_at_t <- state_data$odds_frequency[t_idx]
    
    # Handle NA values
    if(is.na(vax_rate_at_t)) vax_rate_at_t <- 0
    if(is.na(humidity_at_t)) humidity_at_t <- 0.5
    if(is.na(temp_at_t)) temp_at_t <- 0.05
    if(is.na(variant_at_t)) variant_at_t <- 0.5
    
    # Seasonal function
    seasonal <- beta * (var * variant_at_t + hum * humidity_at_t + temp * temp_at_t)
    
    # Ensure seasonal is positive and finite
    if(!is.finite(seasonal) || seasonal < 0) seasonal <- 0.01
    
    # Differential equations
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
  
  likelihood_func <- function(params) {
    
    total_log_lik <- 0
    
    for(i in 1:n_states) {
      
      state_name <- states[i]
      state_data <- full_calibration_dat[full_calibration_dat$Location == state_name, ]
      state_data <- state_data[order(state_data$t), ]
      
#      cat("\n=============================\n")
#      cat("Processing state:", state_name, "\n")
   #   cat("Data rows:", nrow(state_data), "columns:", ncol(state_data), "\n")
  #    cat("Head of state data:\n")
  #    print(head(state_data))
      
      # State-specific parameters
      beta_param <- paste0("beta_", state_name)
      S0_param <- paste0("S0_", state_name)
      
      beta_value <- params[[beta_param]]
      S0_value <- params[[S0_param]]
      
      N_val <- state_populations[[state_name]]
      
      state_init <- c(
        S = S0_value * N_val,
        E = 0.01 * N_val,
        I = 0.001 * N_val,
        R = (1 - S0_value - 0.01 - 0.001 - 0.049) * N_val,
        V = 0.049 * N_val,
        S2 = 0, E2 = 0, I2 = 0, R2 = 0, V2 = 0
      )
      
      all_params <- list(
        var = params["var"], hum = params["hum"], temp = params["temp"],
        imm = params["imm"], imm_vax = params["imm_vax"],
        hr1 = params["hr1"], hr2 = params["hr2"],
        beta = beta_value, N = N_val,
        lp = 1/4, gamma = 1/4, VE = 0.30,
        state_data = state_data
      )
      
      time_vec <- 1:nrow(state_data)
      
#      cat("About to solve ODE with initial conditions:\n")
#      print(state_init)
#      cat("Time vector length:", length(time_vec), "\n")
      
      # Solve ODE with debug
      out <- tryCatch({
        ode(y = state_init, times = time_vec, func = seirs_model_hybrid, parms = all_params)
      }, error = function(e) {
#        cat("ODE solver ERROR for state", state_name, ":\n")
#       cat("Parameters:\n")
#        print(all_params)
#        cat("Initial state:\n")
#        print(state_init)
#        cat("Error message:", e$message, "\n")
        return(NULL)
      })
      
      if(is.null(out)) {
 #       cat("Skipping state due to ODE failure.\n")
        return(-1e6)
      }
      
#      cat("ODE solved successfully. Head of output:\n")
#      print(tail(out))
      
      out_df <- as.data.frame(out)
      predicted_hosp <- out_df$I * all_params$hr1 + out_df$I2 * all_params$hr2
      observed_hosp <- state_data$hosp
      
 #     cat("Predicted hosp (head):", paste(head(predicted_hosp), collapse=", "), "\n")
#      cat("Observed hosp (head):", paste(head(observed_hosp), collapse=", "), "\n")
      
      # Select subset for likelihood
      start_idx <- min(10, length(observed_hosp) - 50)
      end_idx <- length(observed_hosp)
      
      if(start_idx >= end_idx || length(predicted_hosp) < end_idx || length(observed_hosp) < end_idx) {
#        cat("Dimension mismatch: skipping likelihood calculation for state", state_name, "\n")
        return(-1e6)
      }
      
      pred_subset <- pmax(predicted_hosp[start_idx:end_idx], 1e-6)
      obs_subset <- observed_hosp[start_idx:end_idx]
      
      if(any(is.na(pred_subset)) || any(is.na(obs_subset))) {
#        cat("NA values detected in predictions or observations!\n")
        return(-1e6)
      }
      
      state_log_lik <- sum(dpois(obs_subset, lambda = pred_subset, log = TRUE))
 #     cat(sprintf("State %s log-likelihood = %.2f (mean pred=%.2f, mean obs=%.2f)\n", 
#                  state_name, state_log_lik, mean(pred_subset), mean(obs_subset)))
      
      total_log_lik <- total_log_lik + state_log_lik
    }
    
 #   cat("Total log-likelihood across all states:", total_log_lik, "\n")
    print(total_log_lik)
    return(total_log_lik)
  }
  
  
  
  # Hierarchical prior function
  log_prior <- function(params) {
    
    # Global parameter bounds
    global_bounds <- list(
      var = c(0.01, 2.0),
      hum = c(0.01, 2.0), 
      temp = c(0.01, 1.5),
      imm = c(1/300, 1/90),
      imm_vax = c(1/120, 1/30),
      hr1 = c(0.0001, 0.2),
      hr2 = c(0.0001, 0.2)
    )
    
    # State-specific parameter bounds
    state_bounds <- list(
      beta = c(0.05, 1.50),
      S0 = c(0.30, 0.80)
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
      
      # Check S0
      S0_param <- paste0("S0_", state_name)
      if(S0_param %in% names(params)) {
        bounds <- state_bounds[["S0"]]
        if(params[S0_param] < bounds[1] || params[S0_param] > bounds[2]) {
          return(-Inf)
        }
      }
    }
    
    return(0)  # Uniform priors within bounds
  }
  
  # Posterior function
  log_posterior <- function(params) {
    lp <- log_prior(params)
    
    if(is.infinite(lp)) {
      return(-Inf)
    }
    
    ll <- likelihood_func(params)
    return(lp + ll)
  }
  
  # Initial parameter values
  global_initial <- list(
    var = 0.25,
    hum = 0.30,
    temp = 0.15,
    imm = 1/180,
    imm_vax = 1/75,
    hr1 = 0.008,
    hr2 = 0.003
  )
  
  # State-specific initial values (only fitted parameters)
  state_initial <- list()
  for(state_name in states) {
    state_initial[[paste0("beta_", state_name)]] <- runif(1, 0.2, 0.8)
    state_initial[[paste0("S0_", state_name)]] <- runif(1, 0.45, 0.65)
  }
  
  # Combine initial parameters (only the ones we're fitting)
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
  
  # Function to run a single MCMC chain
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
    
    # Check initial log posterior
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
    
    # Adaptive proposal standard deviations
    proposal_sd <- abs(current_params) * 0.10
    proposal_sd[proposal_sd < 1e-4] <- 1e-3 
    
    # Adaptation parameters
    adapt_start <- 100
    adapt_interval <- 50
    accept_count <- numeric(length(param_names))
    proposal_attempts <- numeric(length(param_names))
    target_accept_rate <- 0.234
    
    sample_idx <- 1
    
    for(i in 1:iterations) {
      
      # Component-wise updates
      for(j in 1:length(param_names)) {
        param_name <- param_names[j]
        
        if(proposal_sd[j] <= 1e-10) next
        
        proposed_params <- current_params
        proposed_params[j] <- current_params[j] + rnorm(1, 0, proposal_sd[j])
        
        proposed_log_post <- log_posterior(proposed_params)
        
        if (i %% 100 == 0 && j == 1) {
          cat(sprintf("Iter %d: Current log post = %.2f\n", i, current_log_post))
        }
        
        if(!is.na(proposed_log_post) && is.finite(proposed_log_post)) {
          log_accept_ratio <- proposed_log_post - current_log_post
          proposal_attempts[j] <- proposal_attempts[j] + 1
          
          if(is.finite(log_accept_ratio) && log(runif(1)) < log_accept_ratio) {
            current_params <- proposed_params
            current_log_post <- proposed_log_post
            accept_count[j] <- accept_count[j] + 1
          }
        }
        
        # Adapt proposal standard deviation
        if(i > adapt_start && i %% adapt_interval == 0 && i <= burn_in) {
          if(proposal_attempts[j] > 0) {
            accept_rate_j <- accept_count[j] / proposal_attempts[j]
            
            if(accept_rate_j < target_accept_rate) {
              proposal_sd[j] <- max(proposal_sd[j] * 0.9, 1e-8)
            } else {
              proposal_sd[j] <- min(proposal_sd[j] * 1.1, 1)
            }
          }
        }
      }
      
      # Store samples after burn-in and thinning
      if(i > burn_in && (i - burn_in) %% thin == 0) {
        samples[sample_idx, ] <- current_params
        sample_idx <- sample_idx + 1
      }
      
      # Print progress
      if(i %% 500 == 0) {
        overall_accept <- sum(accept_count[proposal_attempts > 0]) / sum(proposal_attempts[proposal_attempts > 0])
        cat(sprintf("Chain %d, Iteration %d: Overall acceptance rate = %.3f\n", 
                    chain_id, i, overall_accept))
      }
    }
    
    return(list(
      samples = samples,
      acceptance_rate = sum(accept_count) / sum(proposal_attempts),
      param_accept_rates = accept_count / pmax(proposal_attempts, 1)
    ))
  }
  
  # Run multiple chains
  cat("Starting hierarchical MCMC with", n_chains, "chains,", iterations, "iterations\n")
  
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
  
  # Combine results
  all_samples <- list()
  acceptance_rates <- numeric(n_chains)
  param_accept_rates <- matrix(0, nrow = n_chains, ncol = length(param_names))
  colnames(param_accept_rates) <- param_names
  
  for(i in 1:n_chains) {
    all_samples[[i]] <- chain_results[[i]]$samples
    acceptance_rates[i] <- chain_results[[i]]$acceptance_rate
    param_accept_rates[i, ] <- chain_results[[i]]$param_accept_rates
  }
  
  # Create MCMC objects
  mcmc_objects <- lapply(all_samples, mcmc)
  mcmc_list <- mcmc.list(mcmc_objects)
  
  # Diagnostics - with error handling
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
    param_accept_rates = param_accept_rates,
    gelman_diag = gelman_diag,
    summary = summary_stats,
    states = states,
    param_names = param_names
  ))
}


##########################################################################################
##########################################################################################
# Sep 3 2025

run_hierarchical_mcmc_calibration <- function(full_calibration_dat, iterations = 2000, burn_in = 300, thin = 2, n_chains = 3) {
  
  # Load required library
  if (!requireNamespace("deSolve", quietly = TRUE)) {
    stop("deSolve package is required but not installed. Please install it with: install.packages('deSolve')")
  }
  library(deSolve)
  
  # states to be fit in the calibration
  states <- unique(full_calibration_dat$Location)
  # number of states being fit 
  n_states <- length(states)
  
  cat("Fitting", n_states, "states:", paste(states, collapse = ", "), "\n")
  
  # state populations as of the 2023 census 
  state_populations <- list(
    "CA" = 39198693, "CO" = 5901339, "CT" = 3643023, "GA" = 11064432,
    "MI" = 10083356, "MN" = 5753048, "NM" = 2121164, "NY" = 19737367,
    "OH" = 11824034, "OR" = 4253653, "TN" = 7148304, "UT" = 3380800
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
    var <- parms[["var"]]
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
    
    # pull in time series being fit 
    vax_rate_at_t <- state_data$incident[t_idx]
    humidity_at_t <- state_data$humid_smooth[t_idx]
    temp_at_t <- state_data$inverse_temp[t_idx]
    variant_at_t <- state_data$odds_frequency[t_idx]
    
   # Handle NA values
    if(is.na(vax_rate_at_t)) vax_rate_at_t <- 0
     if(is.na(humidity_at_t)) humidity_at_t <- 0.5
     if(is.na(temp_at_t)) temp_at_t <- 0.05
     if(is.na(variant_at_t)) variant_at_t <- 0.5
    
    # seasonal function
    seasonal <- beta * (var * variant_at_t + hum * humidity_at_t + temp * temp_at_t)
    
    # make seasonal value something small and positive if negative- this should not occur to start with 
    if(!is.finite(seasonal) || seasonal < 0) seasonal <- 0.01
    
    #  diffy q's 
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
      state_data <- full_calibration_dat[full_calibration_dat$Location == state_name, ]
      state_data <- state_data[order(state_data$t), ]
      
      # state-specific parameters
      beta_param <- paste0("beta_", state_name)
      S0_param <- paste0("S0_", state_name)
      
      beta_value <- params[[beta_param]]
      S0_value <- params[[S0_param]]
      
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
        var = params["var"], hum = params["hum"], temp = params["temp"],
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
      
      # discard first few samples in case of initial dynamics 
      start_idx <- min(10, length(observed_hosp) - 50)
      end_idx <- length(observed_hosp)
      
      if(start_idx >= end_idx || length(predicted_hosp) < end_idx || length(observed_hosp) < end_idx) {
        return(-1e6)
      }
      
      pred_subset <- pmax(predicted_hosp[start_idx:end_idx], 1e-6)
      obs_subset <- observed_hosp[start_idx:end_idx]
      
      if(any(is.na(pred_subset)) || any(is.na(obs_subset))) {
        return(-1e6)
      }
      
      # calculate likelihood using poission distribution 
      state_log_lik <- sum(dpois(obs_subset, lambda = pred_subset, log = TRUE))
      
      # iteratively add likeliloods from across states 
      total_log_lik <- total_log_lik + state_log_lik
    }
    
    print(total_log_lik)
    return(total_log_lik)
  }
  
  
  # prior function 
  log_prior <- function(params) {
    
    # global parameters 
    global_bounds <- list(
      var = c(0.01, 2.0),
      hum = c(0.01, 2.0), 
      temp = c(0.01, 1.5),
      imm = c(1/50, 1/12),
      imm_vax = c(1/20, 1/6),
      hr1 = c(0.0001, 0.2),
      hr2 = c(0.0001, 0.2)
    )
    
    # state-specific parameters 
    state_bounds <- list(
      beta = c(0.05, 0.80),
      S0 = c(0.30, 0.80)
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
      
      # Check S0
      S0_param <- paste0("S0_", state_name)
      if(S0_param %in% names(params)) {
        bounds <- state_bounds[["S0"]]
        if(params[S0_param] < bounds[1] || params[S0_param] > bounds[2]) {
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
    var = 0.45,
    hum = 0.30,
    temp = 0.15,
    imm = 1/24,
    imm_vax = 1/8,
    hr1 = 0.008,
    hr2 = 0.003
  )
  
  # state-specific initial values (only fitted parameters)
  state_initial <- list()
  for(state_name in states) {
    state_initial[[paste0("beta_", state_name)]] <- runif(1, 0.2, .6)
    state_initial[[paste0("S0_", state_name)]] <- runif(1, 0.25, 0.75)
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
  
  # Function to run a single MCMC chain with JOINT PROPOSALS
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
    
    # Initialize proposal covariance matrix
    n_params <- length(param_names)
    
    # Start with diagonal covariance based on parameter scales
    # Builds a diagonal covariance matrix for the multivariate normal proposal distribution.
    # Each diagonal entry is the variance of proposals for that parameter (scale^2).
    # Off-diagonal terms are zero → no correlation between parameters in the initial proposal.
    # tells you how related param are to each other 
    
    initial_scales <- abs(current_params) * 0.05  # Start with smaller proposals for joint updates
    initial_scales[initial_scales < 1e-5] <- 1e-4
    proposal_cov <- diag(initial_scales^2)
    rownames(proposal_cov) <- param_names
    colnames(proposal_cov) <- param_names
    
    # Adaptation parameters for covariance
    adapt_start <- max(100, n_params * 2)  # Start adaptation after enough samples
    adapt_interval <- 50
    accept_count <- 0
    total_proposals <- 0
    target_accept_rate <- 0.234  # Optimal for multivariate normal
    
    # For adaptive covariance estimation
    sample_history <- matrix(NA, nrow = 0, ncol = n_params)
    colnames(sample_history) <- param_names
    
    sample_idx <- 1
    
    # Cholesky decomposition for multivariate normal sampling
    # essentially, if two parameters are correlated, propose parmeters guesses 
    # along the ridge, instead of proposing parameters zig zagging across space 
    chol_cov <- chol(proposal_cov)
    
    for(i in 1:iterations) {
      
      # JOINT PROPOSAL: Propose all parameters at once
      tryCatch({
        # Generate multivariate normal proposal
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
        
        # Store current parameters for covariance adaptation
        if(i > adapt_start && i <= burn_in) {
          sample_history <- rbind(sample_history, current_params)
          
          # Adapt covariance matrix periodically
          if(i %% adapt_interval == 0 && nrow(sample_history) > n_params) {
            
            # Calculate empirical covariance
            empirical_cov <- cov(sample_history)
            
            # Scale by 2.38^2 / d (optimal scaling)
            scaling_factor <- (2.38^2) / n_params
            proposal_cov <- scaling_factor * empirical_cov
            
            # Add small diagonal component for numerical stability
            proposal_cov <- proposal_cov + diag(1e-6, n_params)
            
            # Check if matrix is positive definite
             # eigenvalues of the covariance matrix must be positive (or not multivariate sampling)
            eigenvals <- eigen(proposal_cov, only.values = TRUE)$values
            if(any(eigenvals <= 0)) {
              # Fall back to diagonal matrix if not positive definite
              proposal_cov <- diag(diag(empirical_cov)) * scaling_factor
              proposal_cov <- proposal_cov + diag(1e-6, n_params)
            }
            
            # Update Cholesky decomposition
            tryCatch({
              chol_cov <- chol(proposal_cov)
            }, error = function(e) {
              # If Cholesky fails, use diagonal
              proposal_cov <<- diag(diag(proposal_cov))
              chol_cov <<- chol(proposal_cov)
            })
            
            # Current acceptance rate
            current_accept_rate <- accept_count / total_proposals
            
            # Additional scaling based on acceptance rate
            if(current_accept_rate < 0.15) {
              # Too low acceptance, make proposals smaller
              chol_cov <- chol_cov * 0.9
            } else if(current_accept_rate > 0.35) {
              # Too high acceptance, make proposals larger  
              chol_cov <- chol_cov * 1.1
            }
            
            cat(sprintf("Chain %d, Iter %d: Adapted covariance, accept rate = %.3f\n", 
                        chain_id, i, current_accept_rate))
          }
        }
        
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
        cat(sprintf("Chain %d, Iteration %d: Acceptance rate = %.3f\n", 
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
  cat("Starting hierarchical MCMC with JOINT PROPOSALS,", n_chains, "chains,", iterations, "iterations\n")
  
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
  iterations = 2000,     # start small for testing
  burn_in = 400,         # burn-in period
  thin = 1,              # thinning interval
  n_chains = 2           # number of MCMC chains
)

#############################################################
### Post-processing MCMC results 
head(results)

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

# Extract specific parameter types
global_params <- results$summary$statistics[
  rownames(results$summary$statistics) %in% c("var", "hum", "temp", "imm", "imm_vax", "hr1", "hr2"), 
]
print("Global parameter estimates:")
print(global_params)

# Extract state-specific beta parameters
beta_params <- results$summary$statistics[
  grepl("beta_", rownames(results$summary$statistics)), 
]
print("State-specific beta estimates:")
print(beta_params)

# Extract state-specific S0 parameters
S0_params <- results$summary$statistics[
  grepl("S0_", rownames(results$summary$statistics)), 
]
print("State-specific S0 estimates:")
print(S0_params)

# Plot trace plots for key parameters (optional)
# You can plot traces to check convergence visually
if(requireNamespace("coda", quietly = TRUE)) {
  # Plot traces for global parameters
  plot(results$mcmc_list[, "var"], main = "Trace plot for 'var' parameter")
  plot(results$mcmc_list[, "imm"], main = "Trace plot for 'imm' parameter")
  
  # Plot trace for a state-specific parameter
  first_state <- results$states[10]
  beta_param_name <- paste0("beta_", first_state)
  if(beta_param_name %in% colnames(results$mcmc_list[[1]])) {
    plot(results$mcmc_list[, beta_param_name], 
         main = paste("Trace plot for", beta_param_name))
  }
}


# check a fit 
head(calibration_dataset)

seirs_model_hybrid = function(t, state, param) {
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
    vax_rate_at_t = calibration_dataset$incident[t] 
    humidity_at_t = calibration_dataset$humid_smooth[t] 
    temp_at_t = calibration_dataset$inverse_temp[t]
    variant_at_t = calibration_dataset$odds_frequency[t] 
    
    
    # Seasonal function 
    seasonal = beta*( var*variant_at_t + hum*humidity_at_t + temp*temp_at_t)
    
    print(seasonal)
    
    dS = -seasonal*S*(I + I2)/N - VE*vax_rate_at_t*S + imm_vax*V
    dE =  seasonal*S*(I+ I2)/N - lp*E
    dI =  lp*E - gamma*I
    dR = I*gamma -  R*imm   
    dV = VE*vax_rate_at_t*S  - imm_vax*V   
    
    dS2 = -seasonal*S2*(I + I2)/N  + imm*R  + imm*R2 + imm_vax*V2 - VE*vax_rate_at_t*S2
    dE2 =  seasonal*S2*(I + I2)/N - lp*E2
    dI2 =  lp*E2 - gamma*I2
    dR2 = I2*gamma -  R2*imm
    dV2 = VE*vax_rate_at_t*S2  - imm_vax*V2  
    
    # Return the rates of change
    list(c(dS, dE, dI, dR, dV, dS2, dE2, dI2, dR2, dV2))
  })
}


##############################
# loop to plot all outputs 
model_outputs = left_join(full_calibration_dat, census_dat, by = "Location")
head(model_outputs)

tab = data.frame(results$summary$statistics)
head(tab)

# initalize plots
plot_list <- list()

state = print(unique(model_outputs$Location))

for(i in state){
  
  dat = model_outputs %>%
    filter(Location == i)
  
  tab_state <- results$summary$statistics[grepl(paste0(i, "$"), rownames(results$summary$statistics)), ]
  
  # states 
  N <- print(unique(dat$pop_size_2023))
    
  #S0 <- .13*N
  # S0 <- (tab_state[2][i]) **N
  S0 <- as.numeric(tab_state[grepl("S0_", rownames(tab_state)), "Mean"]) * N
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
    var = results$summary$statistics[1],   
    hum = results$summary$statistics[2], 
    temp = results$summary$statistics[3],  
    imm = results$summary$statistics[4], 
    imm_vax  = results$summary$statistics[5],
    hr =  results$summary$statistics[6],
    hr2 = results$summary$statistics[7],
 #   beta =  tab_state1[1][i],
    beta = as.numeric(tab_state[grepl("beta_", rownames(tab_state)), "Mean"]),
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
      variant_at_t = dat$odds_frequency[t] 
      
      # Seasonal function 
      seasonal = beta*( var*variant_at_t + hum*humidity_at_t + temp*temp_at_t)
      
      dS = -seasonal*S*(I + I2)/N - VE*vax_rate_at_t*S + imm_vax*V
      dE =  seasonal*S*(I+ I2)/N - lp*E
      dI =  lp*E - gamma*I
      dR = I*gamma -  R*imm   
      dV = VE*vax_rate_at_t*S  - imm_vax*V   
      
      dS2 = -seasonal*S2*(I + I2)/N  + imm*R  + imm*R2 + imm_vax*V2 - VE*vax_rate_at_t*S2
      dE2 =  seasonal*S2*(I + I2)/N - lp*E2
      dI2 =  lp*E2 - gamma*I2
      dR2 = I2*gamma -  R2*imm
      dV2 = VE*vax_rate_at_t*S2  - imm_vax*V2  
      
      # Return the rates of change
      list(c(dS, dE, dI, dR, dV, dS2, dE2, dI2, dR2, dV2))
    })
  }
  
  # Solve the ODE system
  out_check <- ode(y = state, times = time, func = seirs_model_hybrid_state, parms = param)
  # hop
  out_df_check <- as.data.frame(out_check) %>%
    mutate(predicted_hosp = I* results$summary$statistics[6] + I2* results$summary$statistics[7] )
  
  # plot 
  plot =  ggplot() +
    geom_line(data = out_df_check %>% filter(time > 10), aes(x = time, y = predicted_hosp), col = "red", lwd = 1.4) +
    geom_line(data = dat %>% filter(t > 10), aes(x = t, y = hosp), lwd = 1.4) +
    ggtitle(paste("State:", i)) +
    ylab("Hospitalizations") +
    theme_bw()
  
  plot_list[[i]] <- plot
  
}

wrap_plots(plotlist = plot_list, ncol = 4)



