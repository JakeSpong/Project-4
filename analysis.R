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

#### Lysate pH ----
pH_lys <- read_csv("Data/Lysate pH.csv")

#reorder the sites so they show up on the plot from west (LHS) to east (RHS)
pH_lys$Bracken <- factor(pH_lys$Bracken, levels = c("Present", "Absent"))
#reorder the habitats so they appear in the specified order
pH_lys$Habitat <- factor(pH_lys$Habitat, levels = c("Rainfall", "Grassland", "Heathland"))

#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
pH_lys_bxp <- ggboxplot(pH_lys, x = "Habitat", aes(y = `Lysate pH`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  expression("Lysate pH")
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

show(pH_lys_bxp)  
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_pH-lysate.svg"), width = 10, height= 5, pH_lys_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(pH_lys$`Lysate pH` ~ pH_lys$Bracken*pH_lys$Habitat)
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
anova <- aov(results_df$`CO2 flux (micromol CO2 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_29-09-2025.svg"), width = 10, height= 5, ch4_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CH4 flux (nmol CH4 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
anova <- aov(results_df$`CO2 flux (micromol CO2 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_01-10-2025.svg"), width = 10, height= 5, ch4_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CH4 flux (nmol CH4 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
anova <- aov(results_df$`CO2 flux (micromol CO2 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_03-10-2025.svg"), width = 10, height= 5, ch4_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CH4 flux (nmol CH4 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
anova <- aov(results_df$`CO2 flux (micromol CO2 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_06-10-2025.svg"), width = 10, height= 5, ch4_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CH4 flux (nmol CH4 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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












#### Analyse flux data based on time segments- 08-10-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/08-10-2025/gas-flux-08-10-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-08-10-2025.csv"))

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_08-10-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CO2 flux (micromol CO2 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_08-10-2025.svg"), width = 10, height= 5, ch4_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CH4 flux (nmol CH4 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
write_csv(results_df, "Processed Data/gas_flux_08-10-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 10-10-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/10-10-2025/gas-flux-10-10-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-10-10-2025.csv"))

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_10-10-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CO2 flux (micromol CO2 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_10-10-2025.svg"), width = 10, height= 5, ch4_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CH4 flux (nmol CH4 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
write_csv(results_df, "Processed Data/gas_flux_10-10-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 13-10-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/13-10-2025/gas-flux-13-10-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
#ensure the data is viewed as number, not character
flux_data$`[CO2]d_ppm` <- as.numeric(as.character(flux_data$`[CO2]d_ppm`))
flux_data$`[CH4]d_ppm` <- as.numeric(as.character(flux_data$`[CH4]d_ppm`))

# Combine Date and Time safely (for dd/mm/yyyy format with fractional seconds)
#correct the date from 08 to 13 - 10 -2025
flux_data <- flux_data %>%
  mutate(DateTime = as.POSIXct(
    sub("^08/10/2025", "13/10/2025", Time),
    format = "%d/%m/%Y %H:%M:%OS", tz = "UTC"
  ))

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
  here::here("Data", "Gas Flux Measurements", "measurement-times-13-10-2025.csv"))

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_13-10-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CO2 flux (micromol CO2 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_13-10-2025.svg"), width = 10, height= 5, ch4_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CH4 flux (nmol CH4 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
write_csv(results_df, "Processed Data/gas_flux_13-10-2025_slopes_output.csv")












#### Analyse flux data based on time segments- 15-10-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/15-10-2025/gas-flux-15-10-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-15-10-2025.csv"))

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
      if (p_value > 0.05) {
        p <- ggplot(data_chunk, aes(x = Time_sec, y = .data[[gas]])) +
          geom_point(color = "blue") +
          geom_smooth(method = "lm", se = FALSE, color = "red", na.rm = TRUE) +
          labs(
            title = paste("High p value:", round(p_value, 3), "| Sample:", sample_id, "| Gas:", gas_label),
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
#see spread of p values
hist(results_df$CO2_p_value)
hist(results_df$CH4_p_value)

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_15-10-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CO2 flux (micromol CO2 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_15-10-2025.svg"), width = 10, height= 5, ch4_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CH4 flux (nmol CH4 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
write_csv(results_df, "Processed Data/gas_flux_15-10-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 17-10-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/17-10-2025/gas-flux-17-10-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-17-10-2025.csv"))

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
      if (p_value > 0.05) {
        p <- ggplot(data_chunk, aes(x = Time_sec, y = .data[[gas]])) +
          geom_point(color = "blue") +
          geom_smooth(method = "lm", se = FALSE, color = "red", na.rm = TRUE) +
          labs(
            title = paste("High p value:", round(p_value, 3), "| Sample:", sample_id, "| Gas:", gas_label),
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
#see spread of p values
hist(results_df$CO2_p_value)
hist(results_df$CH4_p_value)

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_17-10-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CO2 flux (micromol CO2 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_17-10-2025.svg"), width = 10, height= 5, ch4_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CH4 flux (nmol CH4 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
write_csv(results_df, "Processed Data/gas_flux_17-10-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 20-10-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/20-10-2025/gas-flux-20-10-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-20-10-2025.csv"))

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
      if (p_value > 0.05) {
        p <- ggplot(data_chunk, aes(x = Time_sec, y = .data[[gas]])) +
          geom_point(color = "blue") +
          geom_smooth(method = "lm", se = FALSE, color = "red", na.rm = TRUE) +
          labs(
            title = paste("High p value:", round(p_value, 3), "| Sample:", sample_id, "| Gas:", gas_label),
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
#see spread of p values
hist(results_df$CO2_p_value)
hist(results_df$CH4_p_value)

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_20-10-2025.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CO2 flux (micromol CO2 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_20-10-2025.svg"), width = 10, height= 5, ch4_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$`CH4 flux (nmol CH4 per s per m2)` ~ results_df$Bracken*results_df$Habitat)
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
write_csv(results_df, "Processed Data/gas_flux_20-10-2025_slopes_output.csv")










#### read in processed gas flux data, combine to one dataframe ----

fluxes_29092025 <- read_csv("Processed Data/gas_flux_29-09-2025_slopes_output.csv")
fluxes_01102025 <- read_csv("Processed Data/gas_flux_01-10-2025_slopes_output.csv")
fluxes_03102025 <- read_csv("Processed Data/gas_flux_03-10-2025_slopes_output.csv")
fluxes_06102025 <- read_csv("Processed Data/gas_flux_06-10-2025_slopes_output.csv")
fluxes_08102025 <- read_csv("Processed Data/gas_flux_08-10-2025_slopes_output.csv")
fluxes_10102025 <- read_csv("Processed Data/gas_flux_10-10-2025_slopes_output.csv")
fluxes_13102025 <- read_csv("Processed Data/gas_flux_13-10-2025_slopes_output.csv")
fluxes_15102025 <- read_csv("Processed Data/gas_flux_15-10-2025_slopes_output.csv")
fluxes_17102025 <- read_csv("Processed Data/gas_flux_17-10-2025_slopes_output.csv")
fluxes_20102025 <- read_csv("Processed Data/gas_flux_20-10-2025_slopes_output.csv")

#combine the datasets
#combine the datasets
flux_data <- rbind(fluxes_29092025, fluxes_01102025,fluxes_03102025,fluxes_06102025, fluxes_08102025,fluxes_10102025,fluxes_13102025,fluxes_15102025,fluxes_17102025,fluxes_20102025)

#### plot average co2 and ch4 flux over time for each habitat ----
#remove unnecessary columns from dataframe

#remove data entries with p values < 0.05
#flux_data <- flux_data[flux_data$CO2_p_value <= 0.05,]
flux_timeseries <- flux_data[c(1,2, 12, 13, 28, 29)]

#get the date only in one column
flux_timeseries <- flux_timeseries %>%
  mutate(date = as.Date(start_time))

#group by date, Habitat, and Bracken, and calculate means
flux_avg <- flux_timeseries %>%
  group_by(date, Habitat, Bracken) %>%
  summarise(
    `Average CO2 flux` = mean(`CO2 flux (micromol CO2 per s per m2)`, na.rm = TRUE),
    `SD CO2 flux` = sd(`CO2 flux (micromol CO2 per s per m2)`, na.rm = TRUE),
    `Average CH4 flux` = mean(`CH4 flux (nmol CH4 per s per m2)`, na.rm = TRUE),
    `SD CH4 flux` = sd(`CH4 flux (nmol CH4 per s per m2)`, na.rm = TRUE),
    .groups = 'drop'
  )


# Create a combined grouping label
flux_avg <- flux_avg %>%
  mutate(Habitat_Bracken = paste(Habitat, Bracken, sep = " - "))

# plot the co2 flux over time
co2_avg_tms <- ggplot(flux_avg, aes(x = date, y = `Average CO2 flux`, color = Habitat_Bracken)) +
  geom_line() +
  geom_point()  +
  geom_errorbar(
    aes(ymin = `Average CO2 flux` - `SD CO2 flux`, 
        ymax = `Average CO2 flux` + `SD CO2 flux`),
    width = 0.2
  ) +
  labs(
    title = "Average CO2 Flux Over Time by Habitat and Bracken Presence",
    x = "Date",
    y = expression("Mean CO"[2] * " flux (" * mu * "mol" ~ s^{-1} ~ m^{-2} * ")"),
    color = "Group"
  ) +
  theme_minimal()
#show the plot
show(co2_avg_tms)
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux-timeseries.svg"), width = 10, height= 5, co2_avg_tms)

# plot the ch4 flux over time
ch4_avg_tms <- ggplot(flux_avg, aes(x = date, y = `Average CH4 flux`, color = Habitat_Bracken)) +
  geom_line() +
  geom_point() +
  geom_errorbar(
    aes(ymin = `Average CH4 flux` - `SD CH4 flux`, 
        ymax = `Average CH4 flux` + `SD CH4 flux`),
    width = 0.2
  ) +
  labs(
    title = "Average CH4 Flux Over Time by Habitat and Bracken Presence",
    x = "Date",
    y =  expression("Mean CH"[4] * " flux (" * "n" * "mol" ~ s^{-1} ~ m^{-2} * ")"),
    color = "Group"
  ) +
  theme_minimal()
#show the boxplot
show(ch4_avg_tms)

#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux-timeseries.svg"), width = 10, height= 5, ch4_avg_tms)

#linear regression to see if time affects gas fluxes

#first we need to convert the date to days since the start of the experiment
flux_avg <- flux_avg %>%
  mutate(
    days_since_start = as.numeric(date - min(date))
  )

#analyse effect of time (i.e. moisture loss) on average co2 flux
lm_co2_time <- lm(`Average CO2 flux` ~ days_since_start*Habitat*Bracken, data = flux_avg)
summary(lm_co2_time)

#analyse effect of time (i.e. moisture loss) on average ch4 flux
lm_ch4_time <- lm(`Average CH4 flux` ~ days_since_start*Habitat*Bracken, data = flux_avg)
summary(lm_ch4_time)



#### pool all time points and compare habitat types ----


#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
co2_bxp_tms <- ggboxplot(flux_avg, x = "Habitat", aes(y = `Average CO2 flux`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  expression("Mean CO"[2] * " flux (" * mu * "mol" ~ s^{-1} ~ m^{-2} * ")")
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

show(co2_bxp_tms)  
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-mean-flux_bxp.svg"), width = 10, height= 5, co2_bxp_tms)

#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
ch4_bxp_tms <- ggboxplot(flux_avg, x = "Habitat", aes(y = `Average CH4 flux`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  expression("Mean CH"[4] * " flux (" * "n" * "mol" ~ s^{-1} ~ m^{-2} * ")")
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

show(ch4_bxp_tms)  
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-mean-flux_bxp.svg"), width = 10, height= 5, ch4_bxp_tms)

#is there a different in fluxes between habitat types? - co2

#Type 1 two-way anova using data from all sites
anova <- aov(flux_avg$`Average CO2 flux` ~ flux_avg$Habitat*flux_avg$Bracken)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(flux_avg$`Average CO2 flux` ~ flux_avg$Habitat*flux_avg$Bracken)
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

#is there a different in fluxes between habitat types? - ch4

#Type 1 two-way anova using data from all sites
anova <- aov(flux_avg$`Average CH4 flux` ~ flux_avg$Habitat*flux_avg$Bracken)
#check homogeneity of variance
plot(anova, 1)
#levene test.  if p value < 0.05, there is evidence to suggest that the variance across groups is statistically significantly different.
leveneTest(flux_avg$`Average CO2 flux` ~ flux_avg$Habitat*flux_avg$Bracken)
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



#### filter for co2 r2 > 0.9 ----
results_df <- flux_data[flux_data$CO2_p_value < 0.05,]

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_compiled.svg"), width = 10, height= 5, co2_bxp)

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





#### power analysis ----
library(pwr)
pwr.f2.test(k=4,f=.25,sig.level=.05,power=.8)
#u = number of factors we have, f2 = R2/(1 - R2), 
pwr.f2.test(u = 2,  f2 = 0.3/(1-0.3), sig.level = 0.05, power = 0.8)
