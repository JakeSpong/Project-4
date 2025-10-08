# Load libraries
library(tidyverse)
library(here)
library(ggplot2)
library(lubridate)
library(ggpubr)
library(multcompView) #for significant difference letters
library(scales)
library(dplyr)
library(rlang)


#### Analyse flux data based on time segments- 29-09-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/29-09-2025/gas-flux-29-09-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
#ensure the data is viewed as number, not character
flux_data$`[CO2]d_ppm` <- as.numeric(as.character(flux_data$`[CO2]d_ppm`))
flux_data$`[CH4]d_ppm` <- as.numeric(as.character(flux_data$`[CH4]d_ppm`))

# Combine Date and Time safely (for dd/mm/yyyy format with fractional seconds)
flux_data <- flux_data %>%
  mutate(DateTime = as.POSIXct(paste(Time), format = "%d/%m/%Y %H:%M:%OS", tz = "UTC"))

#remove all rows with time 01/01/2002
flux_data <- flux_data[as.Date(flux_data$DateTime) != as.Date("2002-01-01"), ]


#ggplot(flux_data, aes(x = Time, y = `[CO2]d_ppm`, group = 1)) +
# geom_line(color = "blue") +
#geom_point(color = "red", size = 1.5) +
#labs(
# title = "[CO2]d_ppm over Time",
#x = "Time",
#y = "[CO2]d_ppm"
#) 

# --- Load time ranges with sample IDs ---
time_ranges <- readr::read_delim(
  here::here("Data", "Gas Flux Measurements", "measurement-times-29-09-2025.csv"))

# Combine Date + start_time and Date + end_time into POSIXct
# (Assuming the Date column matches the start/end times)
time_ranges <- time_ranges %>%
  mutate(
    start_datetime = as.POSIXct(paste(Date, Start), format = "%d/%m/%Y %H:%M:%OS", tz = "UTC"),
    end_datetime   = as.POSIXct(paste(Date, End),   format = "%d/%m/%Y %H:%M:%OS", tz = "UTC")
  )

#add one minute to the start time so see if that improves r2
time_ranges$start_datetime <- time_ranges$start_datetime + 25

# Initialize list to store results
results_list <- list()

# --- Loop through each time range ---
for (i in 1:nrow(time_ranges)) {
  sample_id <- time_ranges$Mesocosm[i]
  start_dt <- time_ranges$start_datetime[i]
  end_dt <- time_ranges$end_datetime[i]
  
  # Subset flux data
  data_chunk <- flux_data %>%
    filter(DateTime >= start_dt & DateTime <= end_dt) %>%
    mutate(Time_sec = as.numeric(DateTime - min(DateTime)))
  
  # Initialize result row
  result_row <- tibble(
    sample_id = sample_id,
    start_time = start_dt,
    end_time = end_dt
  )
  
  if (nrow(data_chunk) >= 2) {
    
    for (gas in c("[CO2]d_ppm", "[CH4]d_ppm")) {
      
      gas_label <- ifelse(gas == "[CO2]d_ppm", "CO2", "CH4")
      
      # Fit model
      model <- lm(as.formula(paste0("`", gas, "` ~ Time_sec")), data = data_chunk)
      model_summary <- summary(model)
      
      slope <- coef(model)[["Time_sec"]]
      r_squared <- model_summary$r.squared
      adj_r_squared <- model_summary$adj.r.squared
      p_value <- model_summary$coefficients["Time_sec", "Pr(>|t|)"]
      
      # Plot if R² < 0.9
      if (r_squared < 0.9) {
        p <- ggplot(data_chunk, aes(x = Time_sec, y = .data[[gas]])) +
          geom_point(color = "blue") +
          geom_smooth(method = "lm", se = FALSE, color = "red", na.rm = TRUE) +
          labs(
            title = paste("Low R²:", round(r_squared, 3), "| Sample:", sample_id, "| Gas:", gas_label),
            x = "Time (sec)",
            y = paste0(gas_label, " (ppm)")
          ) +
          scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
          theme_minimal() +
          theme(axis.text.y = element_text(angle = 0, hjust = 1))
        
        print(p)
      }
      
      # Add results to row under appropriate column names
      result_row[[paste0(gas_label, "_slope")]] <- slope
      result_row[[paste0(gas_label, "_r2")]] <- r_squared
      result_row[[paste0(gas_label, "_adj_r2")]] <- adj_r_squared
      result_row[[paste0(gas_label, "_p_value")]] <- p_value
    }
  }
  
  # Append the result row
  results_list[[i]] <- result_row
}

# --- Final combined results in wide format ---
results_df <- bind_rows(results_list)




# --- View or export ---
print(results_df)

#see spread of r2 values
hist(results_df$CO2_r2)
hist(results_df$CH4_r2)

# load factors and numbers needed for gas flux calculations
factors <- readr::read_delim(
  here::here("Data", "Gas Flux Measurements", "mesocosm_id_factors.csv"))

# Left join to keep all rows from results_df and add matching Habitat and Vegetation
results_df <- results_df %>%
  left_join(factors, by = "sample_id")

#calculate gas flux - co2
results_df$`CO2 flux (micromol CO2 per s per m2)` <- results_df$CO2_slope*((results_df$`p (atm)`*results_df$`V (L)`)/(results_df$`A (m2)`*results_df$`R (L atm mol-1 K-1)`*results_df$`T (K)`))
#calculate gas flux - ch4
results_df$`CH4 flux (nmol CH4 per s per m2)` <- 1000*(results_df$CH4_slope*((results_df$`p (atm)`*results_df$`V (L)`)/(results_df$`A (m2)`*results_df$`R (L atm mol-1 K-1)`*results_df$`T (K)`)))


#reorder the sites so they show up on the plot from west (LHS) to east (RHS)
results_df$Bracken <- factor(results_df$Bracken, levels = c("Present", "Absent"))
#reorder the habitats so they appear in the specified order
results_df$Habitat <- factor(results_df$Habitat, levels = c("Grassland", "Heathland"))

#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
co2_bxp <- ggboxplot(results_df, x = "Habitat", aes(y = `CO2 flux (micromol CO2 per s per m2)`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  expression("CO"[2] * " rate (" * mu * "mol" ~ s^{-1} ~ m^{-2} * ")")
  ) + theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
  ) 

show(co2_bxp)  
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_29-09-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$CO2_slope ~ results_df$Bracken*results_df$Habitat)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(d$pH ~ d$Vegetation*d$Site)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residuals to test for normality
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
#print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)


#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
ch4_bxp <- ggboxplot(results_df, x = "Habitat", aes(y = `CH4 flux (nmol CH4 per s per m2)`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  expression("CH"[4] * " rate (" * "n" * "mol" ~ s^{-1} ~ m^{-2} * ")")
  ) + theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
  ) 

show(ch4_bxp)  
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_29-09-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$CH4_slope ~ results_df$Bracken*results_df$Habitat)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(d$pH ~ d$Vegetation*d$Site)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residuals to test for normality
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
#print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)



# Optional: write to CSV
write_csv(results_df, "Processed Data/gas_flux_29-09-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 01-10-2025 ----


# --- Load flux data ---
flux_data1 <- read_csv("Data/Gas Flux Measurements/01-10-2025/gas-flux-01-10-2025_pt1.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
flux_data2 <- read_csv("Data/Gas Flux Measurements/01-10-2025/gas-flux-01-10-2025_pt2.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)

#combine the datasets
flux_data <- rbind(flux_data1, flux_data2)
#ensure the data is viewed as number, not character
flux_data$`[CO2]d_ppm` <- as.numeric(as.character(flux_data$`[CO2]d_ppm`))
flux_data$`[CH4]d_ppm` <- as.numeric(as.character(flux_data$`[CH4]d_ppm`))

# Combine Date and Time safely (for dd/mm/yyyy format with fractional seconds)
flux_data <- flux_data %>%
  mutate(DateTime = as.POSIXct(paste(Time), format = "%d/%m/%Y %H:%M:%OS", tz = "UTC"))

#remove all rows with time 01/01/2002
flux_data <- flux_data[as.Date(flux_data$DateTime) != as.Date("2002-01-01"), ]


#ggplot(flux_data, aes(x = Time, y = `[CO2]d_ppm`, group = 1)) +
# geom_line(color = "blue") +
#geom_point(color = "red", size = 1.5) +
#labs(
# title = "[CO2]d_ppm over Time",
#x = "Time",
#y = "[CO2]d_ppm"
#) 

# --- Load time ranges with sample IDs ---
time_ranges <- readr::read_delim(
  here::here("Data", "Gas Flux Measurements", "measurement-times-01-10-2025.csv"))

# Combine Date + start_time and Date + end_time into POSIXct
# (Assuming the Date column matches the start/end times)
time_ranges <- time_ranges %>%
  mutate(
    start_datetime = as.POSIXct(paste(Date, Start), format = "%d/%m/%Y %H:%M:%OS", tz = "UTC"),
    end_datetime   = as.POSIXct(paste(Date, End),   format = "%d/%m/%Y %H:%M:%OS", tz = "UTC")
  )

#add one minute to the start time so see if that improves r2
time_ranges$start_datetime <- time_ranges$start_datetime + 25

# Initialize list to store results
results_list <- list()

# --- Loop through each time range ---
for (i in 1:nrow(time_ranges)) {
  sample_id <- time_ranges$Mesocosm[i]
  start_dt <- time_ranges$start_datetime[i]
  end_dt <- time_ranges$end_datetime[i]
  
  # Subset flux data
  data_chunk <- flux_data %>%
    filter(DateTime >= start_dt & DateTime <= end_dt) %>%
    mutate(Time_sec = as.numeric(DateTime - min(DateTime)))
  
  # Initialize result row
  result_row <- tibble(
    sample_id = sample_id,
    start_time = start_dt,
    end_time = end_dt
  )
  
  if (nrow(data_chunk) >= 2) {
    
    for (gas in c("[CO2]d_ppm", "[CH4]d_ppm")) {
      
      gas_label <- ifelse(gas == "[CO2]d_ppm", "CO2", "CH4")
      
      # Fit model
      model <- lm(as.formula(paste0("`", gas, "` ~ Time_sec")), data = data_chunk)
      model_summary <- summary(model)
      
      slope <- coef(model)[["Time_sec"]]
      r_squared <- model_summary$r.squared
      adj_r_squared <- model_summary$adj.r.squared
      p_value <- model_summary$coefficients["Time_sec", "Pr(>|t|)"]
      
      # Plot if R² < 0.9
      if (r_squared < 0.9) {
        p <- ggplot(data_chunk, aes(x = Time_sec, y = .data[[gas]])) +
          geom_point(color = "blue") +
          geom_smooth(method = "lm", se = FALSE, color = "red", na.rm = TRUE) +
          labs(
            title = paste("Low R²:", round(r_squared, 3), "| Sample:", sample_id, "| Gas:", gas_label),
            x = "Time (sec)",
            y = paste0(gas_label, " (ppm)")
          ) +
          scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
          theme_minimal() +
          theme(axis.text.y = element_text(angle = 0, hjust = 1))
        
        print(p)
      }
      
      # Add results to row under appropriate column names
      result_row[[paste0(gas_label, "_slope")]] <- slope
      result_row[[paste0(gas_label, "_r2")]] <- r_squared
      result_row[[paste0(gas_label, "_adj_r2")]] <- adj_r_squared
      result_row[[paste0(gas_label, "_p_value")]] <- p_value
    }
  }
  
  # Append the result row
  results_list[[i]] <- result_row
}

# --- Final combined results in wide format ---
results_df <- bind_rows(results_list)




# --- View or export ---
print(results_df)

#see spread of r2 values
hist(results_df$CO2_r2)
hist(results_df$CH4_r2)

# load factors and numbers needed for gas flux calculations
factors <- readr::read_delim(
  here::here("Data", "Gas Flux Measurements", "mesocosm_id_factors.csv"))

# Left join to keep all rows from results_df and add matching Habitat and Vegetation
results_df <- results_df %>%
  left_join(factors, by = "sample_id")

#calculate gas flux - co2
results_df$`CO2 flux (micromol CO2 per s per m2)` <- results_df$CO2_slope*((results_df$`p (atm)`*results_df$`V (L)`)/(results_df$`A (m2)`*results_df$`R (L atm mol-1 K-1)`*results_df$`T (K)`))
#calculate gas flux - ch4
results_df$`CH4 flux (nmol CH4 per s per m2)` <- 1000*(results_df$CH4_slope*((results_df$`p (atm)`*results_df$`V (L)`)/(results_df$`A (m2)`*results_df$`R (L atm mol-1 K-1)`*results_df$`T (K)`)))


#reorder the sites so they show up on the plot from west (LHS) to east (RHS)
results_df$Bracken <- factor(results_df$Bracken, levels = c("Present", "Absent"))
#reorder the habitats so they appear in the specified order
results_df$Habitat <- factor(results_df$Habitat, levels = c("Grassland", "Heathland"))

#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
co2_bxp <- ggboxplot(results_df, x = "Habitat", aes(y = `CO2 flux (micromol CO2 per s per m2)`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  expression("CO"[2] * " rate (" * mu * "mol" ~ s^{-1} ~ m^{-2} * ")")
  ) + theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
  ) 

show(co2_bxp)  
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_01-10-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$CO2_slope ~ results_df$Bracken*results_df$Habitat)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(d$pH ~ d$Vegetation*d$Site)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residuals to test for normality
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
#print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)


#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
ch4_bxp <- ggboxplot(results_df, x = "Habitat", aes(y = `CH4 flux (nmol CH4 per s per m2)`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  expression("CH"[4] * " rate (" * "n" * "mol" ~ s^{-1} ~ m^{-2} * ")")
  ) + theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
  ) 

show(ch4_bxp)  
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_01-10-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$CH4_slope ~ results_df$Bracken*results_df$Habitat)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(d$pH ~ d$Vegetation*d$Site)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residuals to test for normality
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
#print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)



# Optional: write to CSV
write_csv(results_df, "Processed Data/gas_flux_01-10-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 03-10-2025 ----


# --- Load flux data ---
flux_data1 <- read_csv("Data/Gas Flux Measurements/03-10-2025/gas-flux-03-10-2025_pt1.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
flux_data2 <- read_csv("Data/Gas Flux Measurements/03-10-2025/gas-flux-03-10-2025_pt2.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)

#combine the datasets
flux_data <- rbind(flux_data1, flux_data2)
#ensure the data is viewed as number, not character
flux_data$`[CO2]d_ppm` <- as.numeric(as.character(flux_data$`[CO2]d_ppm`))
flux_data$`[CH4]d_ppm` <- as.numeric(as.character(flux_data$`[CH4]d_ppm`))

# Combine Date and Time safely (for dd/mm/yyyy format with fractional seconds)
flux_data <- flux_data %>%
  mutate(DateTime = as.POSIXct(paste(Time), format = "%d/%m/%Y %H:%M:%OS", tz = "UTC"))

#remove all rows with time 01/01/2002
flux_data <- flux_data[as.Date(flux_data$DateTime) != as.Date("2002-01-01"), ]


#ggplot(flux_data, aes(x = Time, y = `[CO2]d_ppm`, group = 1)) +
# geom_line(color = "blue") +
#geom_point(color = "red", size = 1.5) +
#labs(
# title = "[CO2]d_ppm over Time",
#x = "Time",
#y = "[CO2]d_ppm"
#) 

# --- Load time ranges with sample IDs ---
time_ranges <- readr::read_delim(
  here::here("Data", "Gas Flux Measurements", "measurement-times-03-10-2025.csv"))

# Combine Date + start_time and Date + end_time into POSIXct
# (Assuming the Date column matches the start/end times)
time_ranges <- time_ranges %>%
  mutate(
    start_datetime = as.POSIXct(paste(Date, Start), format = "%d/%m/%Y %H:%M:%OS", tz = "UTC"),
    end_datetime   = as.POSIXct(paste(Date, End),   format = "%d/%m/%Y %H:%M:%OS", tz = "UTC")
  )

#add one minute to the start time so see if that improves r2
time_ranges$start_datetime <- time_ranges$start_datetime + 25

# Initialize list to store results
results_list <- list()

# --- Loop through each time range ---
for (i in 1:nrow(time_ranges)) {
  sample_id <- time_ranges$Mesocosm[i]
  start_dt <- time_ranges$start_datetime[i]
  end_dt <- time_ranges$end_datetime[i]
  
  # Subset flux data
  data_chunk <- flux_data %>%
    filter(DateTime >= start_dt & DateTime <= end_dt) %>%
    mutate(Time_sec = as.numeric(DateTime - min(DateTime)))
  
  # Initialize result row
  result_row <- tibble(
    sample_id = sample_id,
    start_time = start_dt,
    end_time = end_dt
  )
  
  if (nrow(data_chunk) >= 2) {
    
    for (gas in c("[CO2]d_ppm", "[CH4]d_ppm")) {
      
      gas_label <- ifelse(gas == "[CO2]d_ppm", "CO2", "CH4")
      
      # Fit model
      model <- lm(as.formula(paste0("`", gas, "` ~ Time_sec")), data = data_chunk)
      model_summary <- summary(model)
      
      slope <- coef(model)[["Time_sec"]]
      r_squared <- model_summary$r.squared
      adj_r_squared <- model_summary$adj.r.squared
      p_value <- model_summary$coefficients["Time_sec", "Pr(>|t|)"]
      
      # Plot if R² < 0.9
      if (r_squared < 0.9) {
        p <- ggplot(data_chunk, aes(x = Time_sec, y = .data[[gas]])) +
          geom_point(color = "blue") +
          geom_smooth(method = "lm", se = FALSE, color = "red", na.rm = TRUE) +
          labs(
            title = paste("Low R²:", round(r_squared, 3), "| Sample:", sample_id, "| Gas:", gas_label),
            x = "Time (sec)",
            y = paste0(gas_label, " (ppm)")
          ) +
          scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
          theme_minimal() +
          theme(axis.text.y = element_text(angle = 0, hjust = 1))
        
        print(p)
      }
      
      # Add results to row under appropriate column names
      result_row[[paste0(gas_label, "_slope")]] <- slope
      result_row[[paste0(gas_label, "_r2")]] <- r_squared
      result_row[[paste0(gas_label, "_adj_r2")]] <- adj_r_squared
      result_row[[paste0(gas_label, "_p_value")]] <- p_value
    }
  }
  
  # Append the result row
  results_list[[i]] <- result_row
}

# --- Final combined results in wide format ---
results_df <- bind_rows(results_list)




# --- View or export ---
print(results_df)

#see spread of r2 values
hist(results_df$CO2_r2)
hist(results_df$CH4_r2)

# load factors and numbers needed for gas flux calculations
factors <- readr::read_delim(
  here::here("Data", "Gas Flux Measurements", "mesocosm_id_factors.csv"))

# Left join to keep all rows from results_df and add matching Habitat and Vegetation
results_df <- results_df %>%
  left_join(factors, by = "sample_id")

#calculate gas flux - co2
results_df$`CO2 flux (micromol CO2 per s per m2)` <- results_df$CO2_slope*((results_df$`p (atm)`*results_df$`V (L)`)/(results_df$`A (m2)`*results_df$`R (L atm mol-1 K-1)`*results_df$`T (K)`))
#calculate gas flux - ch4
results_df$`CH4 flux (nmol CH4 per s per m2)` <- 1000*(results_df$CH4_slope*((results_df$`p (atm)`*results_df$`V (L)`)/(results_df$`A (m2)`*results_df$`R (L atm mol-1 K-1)`*results_df$`T (K)`)))


#reorder the sites so they show up on the plot from west (LHS) to east (RHS)
results_df$Bracken <- factor(results_df$Bracken, levels = c("Present", "Absent"))
#reorder the habitats so they appear in the specified order
results_df$Habitat <- factor(results_df$Habitat, levels = c("Grassland", "Heathland"))

#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
co2_bxp <- ggboxplot(results_df, x = "Habitat", aes(y = `CO2 flux (micromol CO2 per s per m2)`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  expression("CO"[2] * " rate (" * mu * "mol" ~ s^{-1} ~ m^{-2} * ")")
  ) + theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
  ) 

show(co2_bxp)  
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_03-10-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$CO2_slope ~ results_df$Bracken*results_df$Habitat)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(d$pH ~ d$Vegetation*d$Site)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residuals to test for normality
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
#print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)


#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
ch4_bxp <- ggboxplot(results_df, x = "Habitat", aes(y = `CH4 flux (nmol CH4 per s per m2)`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  expression("CH"[4] * " rate (" * "n" * "mol" ~ s^{-1} ~ m^{-2} * ")")
  ) + theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
  ) 

show(ch4_bxp)  
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_03-10-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$CH4_slope ~ results_df$Bracken*results_df$Habitat)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(d$pH ~ d$Vegetation*d$Site)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residuals to test for normality
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
#print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

# Optional: write to CSV
write_csv(results_df, "Processed Data/gas_flux_03-10-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 06-10-2025 ----


# --- Load flux data ---
flux_data1 <- read_csv("Data/Gas Flux Measurements/06-10-2025/gas-flux-06-10-2025_pt1.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
flux_data2 <- read_csv("Data/Gas Flux Measurements/06-10-2025/gas-flux-06-10-2025_pt2.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
flux_data3 <- read_csv("Data/Gas Flux Measurements/06-10-2025/gas-flux-06-10-2025_pt3.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)

#combine the datasets
flux_data <- rbind(flux_data1, flux_data2, flux_data3)
#ensure the data is viewed as number, not character
flux_data$`[CO2]d_ppm` <- as.numeric(as.character(flux_data$`[CO2]d_ppm`))
flux_data$`[CH4]d_ppm` <- as.numeric(as.character(flux_data$`[CH4]d_ppm`))

# Combine Date and Time safely (for dd/mm/yyyy format with fractional seconds)
flux_data <- flux_data %>%
  mutate(DateTime = as.POSIXct(paste(Time), format = "%d/%m/%Y %H:%M:%OS", tz = "UTC"))

#remove all rows with time 01/01/2002
flux_data <- flux_data[as.Date(flux_data$DateTime) != as.Date("2002-01-01"), ]


#ggplot(flux_data, aes(x = Time, y = `[CO2]d_ppm`, group = 1)) +
# geom_line(color = "blue") +
#geom_point(color = "red", size = 1.5) +
#labs(
# title = "[CO2]d_ppm over Time",
#x = "Time",
#y = "[CO2]d_ppm"
#) 

# --- Load time ranges with sample IDs ---
time_ranges <- readr::read_delim(
  here::here("Data", "Gas Flux Measurements", "measurement-times-06-10-2025.csv"))

# Combine Date + start_time and Date + end_time into POSIXct
# (Assuming the Date column matches the start/end times)
time_ranges <- time_ranges %>%
  mutate(
    start_datetime = as.POSIXct(paste(Date, Start), format = "%d/%m/%Y %H:%M:%OS", tz = "UTC"),
    end_datetime   = as.POSIXct(paste(Date, End),   format = "%d/%m/%Y %H:%M:%OS", tz = "UTC")
  )

#add one minute to the start time so see if that improves r2
time_ranges$start_datetime <- time_ranges$start_datetime + 25

# Initialize list to store results
results_list <- list()

# --- Loop through each time range ---
for (i in 1:nrow(time_ranges)) {
  sample_id <- time_ranges$Mesocosm[i]
  start_dt <- time_ranges$start_datetime[i]
  end_dt <- time_ranges$end_datetime[i]
  
  # Subset flux data
  data_chunk <- flux_data %>%
    filter(DateTime >= start_dt & DateTime <= end_dt) %>%
    mutate(Time_sec = as.numeric(DateTime - min(DateTime)))
  
  # Initialize result row
  result_row <- tibble(
    sample_id = sample_id,
    start_time = start_dt,
    end_time = end_dt
  )
  
  if (nrow(data_chunk) >= 2) {
    
    for (gas in c("[CO2]d_ppm", "[CH4]d_ppm")) {
      
      gas_label <- ifelse(gas == "[CO2]d_ppm", "CO2", "CH4")
      
      # Fit model
      model <- lm(as.formula(paste0("`", gas, "` ~ Time_sec")), data = data_chunk)
      model_summary <- summary(model)
      
      slope <- coef(model)[["Time_sec"]]
      r_squared <- model_summary$r.squared
      adj_r_squared <- model_summary$adj.r.squared
      p_value <- model_summary$coefficients["Time_sec", "Pr(>|t|)"]
      
      # Plot if R² < 0.9
      if (r_squared < 0.9) {
        p <- ggplot(data_chunk, aes(x = Time_sec, y = .data[[gas]])) +
          geom_point(color = "blue") +
          geom_smooth(method = "lm", se = FALSE, color = "red", na.rm = TRUE) +
          labs(
            title = paste("Low R²:", round(r_squared, 3), "| Sample:", sample_id, "| Gas:", gas_label),
            x = "Time (sec)",
            y = paste0(gas_label, " (ppm)")
          ) +
          scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
          theme_minimal() +
          theme(axis.text.y = element_text(angle = 0, hjust = 1))
        
        print(p)
      }
      
      # Add results to row under appropriate column names
      result_row[[paste0(gas_label, "_slope")]] <- slope
      result_row[[paste0(gas_label, "_r2")]] <- r_squared
      result_row[[paste0(gas_label, "_adj_r2")]] <- adj_r_squared
      result_row[[paste0(gas_label, "_p_value")]] <- p_value
    }
  }
  
  # Append the result row
  results_list[[i]] <- result_row
}

# --- Final combined results in wide format ---
results_df <- bind_rows(results_list)




# --- View or export ---
print(results_df)

#see spread of r2 values
hist(results_df$CO2_r2)
hist(results_df$CH4_r2)

# load factors and numbers needed for gas flux calculations
factors <- readr::read_delim(
  here::here("Data", "Gas Flux Measurements", "mesocosm_id_factors.csv"))

# Left join to keep all rows from results_df and add matching Habitat and Vegetation
results_df <- results_df %>%
  left_join(factors, by = "sample_id")

#calculate gas flux - co2
results_df$`CO2 flux (micromol CO2 per s per m2)` <- results_df$CO2_slope*((results_df$`p (atm)`*results_df$`V (L)`)/(results_df$`A (m2)`*results_df$`R (L atm mol-1 K-1)`*results_df$`T (K)`))
#calculate gas flux - ch4
results_df$`CH4 flux (nmol CH4 per s per m2)` <- 1000*(results_df$CH4_slope*((results_df$`p (atm)`*results_df$`V (L)`)/(results_df$`A (m2)`*results_df$`R (L atm mol-1 K-1)`*results_df$`T (K)`)))


#reorder the sites so they show up on the plot from west (LHS) to east (RHS)
results_df$Bracken <- factor(results_df$Bracken, levels = c("Present", "Absent"))
#reorder the habitats so they appear in the specified order
results_df$Habitat <- factor(results_df$Habitat, levels = c("Grassland", "Heathland"))

#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
co2_bxp <- ggboxplot(results_df, x = "Habitat", aes(y = `CO2 flux (micromol CO2 per s per m2)`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  expression("CO"[2] * " rate (" * mu * "mol" ~ s^{-1} ~ m^{-2} * ")")
  ) + theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
  ) 

show(co2_bxp)  
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_06-10-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$CO2_slope ~ results_df$Bracken*results_df$Habitat)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(d$pH ~ d$Vegetation*d$Site)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residuals to test for normality
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
#print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)


#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
ch4_bxp <- ggboxplot(results_df, x = "Habitat", aes(y = `CH4 flux (nmol CH4 per s per m2)`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  expression("CH"[4] * " rate (" * "n" * "mol" ~ s^{-1} ~ m^{-2} * ")")
  ) + theme(
    # Remove panel border
    panel.border = element_blank(),  
    # Remove panel grid lines
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Remove panel background
    panel.background = element_blank(),
    # Add axis line
    axis.line = element_line(colour = "black", linewidth = 0.5),
    #change colour and thickness of axis ticks
    axis.ticks = element_line(colour = "black", linewidth = 0.5),
    #change axis labels colour
    axis.title.x = element_text(colour = "black"),
    axis.title.y = element_text(colour = "black"),
    #change tick labels colour
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
  ) 

show(ch4_bxp)  
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_06-10-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$CH4_slope ~ results_df$Bracken*results_df$Habitat)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(d$pH ~ d$Vegetation*d$Site)
#check normality.  
plot(anova, 2)
#conduct shapiro-wilk test on ANOVA residuals to test for normality
#extract the residuals
aov_residuals <- residuals(object = anova)
#run shapiro-wilk test.  if p > 0.05 the data is normal
shapiro.test(x = aov_residuals)

summary(anova)
#tukey's test to identify significant interactions
tukey <- TukeyHSD(anova)
#print(tukey)
#compact letter display
cld <- multcompLetters4(anova, tukey)
print(cld)

# Optional: write to CSV
write_csv(results_df, "Processed Data/gas_flux_06-10-2025_slopes_output.csv")












#### read in processed gas flux data, filter for data with r2 > 0.9, and do a master plot ----