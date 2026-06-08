#for generating maps of sample locations
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(ggmap)

# Load libraries
library(tidyverse)
library(here)
library(ggplot2)
#library(lubridate)
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
results_df$`CO2 flux (micromol CO2 per s per m2)` <- results_df$CO2_slope*((results_df$`p (atm)`*results_df$`V (L) no rhizons`)/(results_df$`A (m2)`*results_df$`R (L atm mol-1 K-1)`*results_df$`T (K)`))
#calculate gas flux - ch4
results_df$`CH4 flux (nmol CH4 per s per m2)` <- 1000*(results_df$CH4_slope*((results_df$`p (atm)`*results_df$`V (L) no rhizons`)/(results_df$`A (m2)`*results_df$`R (L atm mol-1 K-1)`*results_df$`T (K)`)))


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










#### Analyse flux data based on time segments- 22-10-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/22-10-2025/gas-flux-22-10-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-22-10-2025.csv"))

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_22-10-2025.svg"), width = 10, height= 5, co2_bxp)

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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_22-10-2025.svg"), width = 10, height= 5, ch4_bxp)

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
write_csv(results_df, "Processed Data/gas_flux_22-10-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 24-10-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/24-10-2025/gas-flux-24-10-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-24-10-2025.csv"))

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_24-10-2025.svg"), width = 10, height= 5, co2_bxp)

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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_24-10-2025.svg"), width = 10, height= 5, ch4_bxp)

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
write_csv(results_df, "Processed Data/gas_flux_24-10-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 27-10-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/27-10-2025/gas-flux-27-10-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-27-10-2025.csv"))

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_27-10-2025.svg"), width = 10, height= 5, co2_bxp)

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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_27-10-2025.svg"), width = 10, height= 5, ch4_bxp)

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
write_csv(results_df, "Processed Data/gas_flux_27-10-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 29-10-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/29-10-2025/gas-flux-29-10-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-29-10-2025.csv"))

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_29-10-2025.svg"), width = 10, height= 5, co2_bxp)

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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_29-10-2025.svg"), width = 10, height= 5, ch4_bxp)

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
write_csv(results_df, "Processed Data/gas_flux_29-10-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 31-10-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/31-10-2025/gas-flux-31-10-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-31-10-2025.csv"))

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_31-10-2025.svg"), width = 10, height= 5, co2_bxp)

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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_31-10-2025.svg"), width = 10, height= 5, ch4_bxp)

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
write_csv(results_df, "Processed Data/gas_flux_31-10-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 03-11-2025 ----


# --- Load flux data ---
# --- Load flux data ---
flux_data1 <- read_csv("Data/Gas Flux Measurements/03-11-2025/gas-flux-03-11-2025_pt-1.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
flux_data2 <- read_csv("Data/Gas Flux Measurements/03-11-2025/gas-flux-03-11-2025_pt-2.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)

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
  here::here("Data", "Gas Flux Measurements", "measurement-times-03-11-2025.csv"))

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_03-11-2025.svg"), width = 10, height= 5, co2_bxp)

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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_03-11-2025.svg"), width = 10, height= 5, ch4_bxp)

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
write_csv(results_df, "Processed Data/gas_flux_03-11-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 05-11-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/05-11-2025/gas-flux-05-11-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-05-11-2025.csv"))

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_05-11-2025.svg"), width = 10, height= 5, co2_bxp)

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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_05-11-2025.svg"), width = 10, height= 5, ch4_bxp)

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
write_csv(results_df, "Processed Data/gas_flux_05-11-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 07-11-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/07-11-2025/gas-flux-07-11-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-07-11-2025.csv"))

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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_07-11-2025.svg"), width = 10, height= 5, co2_bxp)

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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_07-11-2025.svg"), width = 10, height= 5, ch4_bxp)

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
write_csv(results_df, "Processed Data/gas_flux_07-11-2025_slopes_output.csv")










#### Analyse flux data based on time segments- 11-11-2025 ----


# --- Load flux data ---
flux_data <- read_csv("Data/Gas Flux Measurements/11-11-2025/gas-flux-11-11-2025.txt", col_types = cols(`[H2O]_ppm` = col_double(), `[CO2]_ppm` = col_double()), skip = 1)
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
  here::here("Data", "Gas Flux Measurements", "measurement-times-11-11-2025.csv"))

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


#reorder the treatments
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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-flux_11-11-2025.svg"), width = 10, height= 5, co2_bxp)

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
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux_11-11-2025.svg"), width = 10, height= 5, ch4_bxp)

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
write_csv(results_df, "Processed Data/gas_flux_11-11-2025_slopes_output.csv")










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
fluxes_22102025 <- read_csv("Processed Data/gas_flux_22-10-2025_slopes_output.csv")
fluxes_24102025 <- read_csv("Processed Data/gas_flux_24-10-2025_slopes_output.csv")
fluxes_27102025 <- read_csv("Processed Data/gas_flux_27-10-2025_slopes_output.csv")
fluxes_29102025 <- read_csv("Processed Data/gas_flux_29-10-2025_slopes_output.csv")
fluxes_31102025 <- read_csv("Processed Data/gas_flux_31-10-2025_slopes_output.csv")
fluxes_03112025 <- read_csv("Processed Data/gas_flux_03-11-2025_slopes_output.csv")
fluxes_05112025 <- read_csv("Processed Data/gas_flux_05-11-2025_slopes_output.csv")
fluxes_07112025 <- read_csv("Processed Data/gas_flux_07-11-2025_slopes_output.csv")
fluxes_11112025 <- read_csv("Processed Data/gas_flux_11-11-2025_slopes_output.csv")

#combine the datasets
#combine the datasets
flux_data <- rbind(fluxes_29092025, fluxes_01102025,fluxes_03102025,fluxes_06102025, fluxes_08102025,fluxes_10102025,fluxes_13102025,fluxes_15102025,fluxes_17102025,fluxes_20102025,fluxes_22102025,fluxes_24102025,fluxes_27102025,fluxes_29102025,fluxes_31102025,fluxes_03112025,fluxes_05112025,fluxes_07112025,fluxes_11112025)

write_csv(flux_data, "Processed Data/all_times_gas_fluxes.csv")


#### generate maps of sample locations (Figure 1) ----
#to get google satellite layer, Go to the Google Cloud Console.  Create a new project (or use an existing one). Enable the Maps Static API and Geocoding API.  Get your API key.

library(maptiles)
library(terra)
library(sf)
library(ggplot2)

d <- readr::read_csv(
  here::here("Data", "Field Data.csv")
) 


# Convert your coordinates to an sf object
coords_sf <- st_as_sf(d, 
                      coords = c("LongitudeE", "LatitudeN"), 
                      crs = 4326)

# Create a buffered bounding box to expand the map extent
bbox_buffered <- st_bbox(coords_sf)
buffer <- 0.002  # degrees — increase this for more surrounding land

bbox_buffered["xmin"] <- bbox_buffered["xmin"] - buffer
bbox_buffered["xmax"] <- bbox_buffered["xmax"] + buffer
bbox_buffered["ymin"] <- bbox_buffered["ymin"] - buffer
bbox_buffered["ymax"] <- bbox_buffered["ymax"] + buffer

# Convert buffered bbox to sf object for get_tiles
bbox_sf <- st_as_sfc(bbox_buffered)

# Get satellite tiles using the expanded extent
satellite_map <- get_tiles(bbox_sf, 
                           provider = "Esri.WorldImagery",
                           zoom = 17,
                           crop = TRUE)
# Convert tiles to data frame for ggplot
map_df <- as.data.frame(satellite_map, xy = TRUE)


# Robust approach — grab the first three non-xy columns regardless of their name
r_col <- names(map_df)[3]
g_col <- names(map_df)[4]
b_col <- names(map_df)[5]
# Plot
maps_figure <- ggplot() +
  geom_raster(data = map_df, 
              aes(x = x, y = y, 
                  fill = rgb(map_df[[r_col]], 
                             map_df[[g_col]], 
                             map_df[[b_col]], 
                             maxColorValue = 255))) +
  scale_fill_identity() +
  
  geom_sf(data = coords_sf,
          aes(shape = Bracken,
              colour  = Habitat),
          size = 3, stroke = 1.2) +
  
  scale_colour_manual(values = c("Grassland" = "green",
                                 "Heathland" = "purple")) +
  
  scale_shape_manual(values = c("Present" = 16,
                                "Absent" = 17)) +
  
  # Override the axis tick labels to plain numbers
  scale_x_continuous(labels = function(x) paste0(round(x, 3), "°E")) +
  scale_y_continuous(labels = function(y) paste0(round(y, 3), "°N")) +
  
  labs(colour = "Habitat",
       shape  = "Bracken",
       x = "Longitude",
       y = "Latitude") +
  
  theme_minimal() +
  theme(legend.position = "right",
        legend.box      = "vertical")


show(maps_figure)
#save the figure
ggsave("Figures/sample_locations_figure.svg", plot = maps_figure, width = 7.5, height = 6, dpi = 300)

# 
# 
# # Load the UK map from the rnaturalearth package
# uk_map <- ne_countries(scale = "medium", returnclass = "sf") %>%
#   filter(name == "United Kingdom")
# 
# #extract he sample coordinates from the dataframe
# sample_coordinates <- data.frame(
#   Longitude = d$LongitudeE,
#   Latitude = d$LatitudeN,
#   Habitat = d$Habitat,
#   Vegetation = d$Vegetation
# )
# # Calculate the bounding box surrounding the 30 sample coordinates
# longitude_range_sample <- range(sample_coordinates$Longitude)
# latitude_range_sample <- range(sample_coordinates$Latitude)
# 
# # Plot the UK map with only the bounding box
# UK_Hawes <- ggplot(data = uk_map) +
#   geom_sf(fill = "lightblue") +  # UK map with light blue fill
#   geom_rect(
#     aes(xmin = longitude_range_sample[1], xmax = longitude_range_sample[2], ymin = latitude_range_sample[1], ymax = latitude_range_sample[2]), 
#     color = "red", fill = NA, linewidth = 2) +  # Add bounding box around the points
#   theme_minimal()
# #show the map
# #show(UK_Hawes)
# #save the UK map with Haweswater highlighted
# #ggsave("Figures/UK_Haweswater_boundingbox.svg", width = 8, height = 6)
# 


#### plot average co2 and ch4 flux over time for each habitat ----
#remove unnecessary columns from dataframe
flux_data <- read_csv("Processed Data/all_times_gas_fluxes.csv")
#remove data entries with p values < 0.05
#flux_data <- flux_data[flux_data$CO2_p_value <= 0.05,]
flux_timeseries <- flux_data[c(1,2, 12, 13, 33, 34)]

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
  ) + geom_vline(xintercept = as.numeric(as.Date("2025-10-14")), 
               linetype = "dashed", color = "blue", linewidth = 0.5) +
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
  geom_vline(xintercept = as.numeric(as.Date("2025-10-14")), 
             linetype = "dashed", color = "blue", linewidth = 0.5) +
  theme_minimal()
#show the boxplot
show(ch4_avg_tms)

#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_ch4-flux-timeseries.svg"), width = 10, height= 5, ch4_avg_tms)

#### Question 1: is there a difference in CO2 and CH4 fluxes between vegetation and habitat?? ----
#remove unnecessary columns from dataframe
flux_data <- read_csv("Processed Data/all_times_gas_fluxes.csv")
#remove data entries with p values < 0.05
#flux_data <- flux_data[flux_data$CO2_p_value <= 0.05,]
flux_timeseries <- flux_data[c(1,2, 12, 13, 33, 34)]

#get the date only in one column
flux_timeseries <- flux_timeseries %>%
  mutate(date = as.Date(start_time))
#use all datapoints, not an everage
# Create a combined grouping label
flux_timeseries <- flux_timeseries %>%
  mutate(Habitat_Bracken = paste(Habitat, Bracken, sep = " - "))
#trim data to just pre-rainfall
#load in the covaraites
cov_data <- read_csv("Data/Covariates.csv")

#append the covariate data to the respective sample ID 

flux_timeseries <- flux_timeseries %>%
  left_join(cov_data, by = "sample_id")

#calculate actual soil moisture, not just relative change

#load the data, format data correctly
moisture_data <- read_csv("Data/Mesocosm Masses - Soil mass change as % of day 0.csv")
moisture_data$date <- as.Date(moisture_data$Date, format = "%d/%m/%Y")
#combine the two dataframes
flux_timeseries <- left_join(flux_timeseries, moisture_data, by = c("sample_id", "date"))
#new column containg actual soil moisture
flux_timeseries$`Estimated soil moisture (% wet mass)` <- flux_timeseries$`Soil moisture (% wet mass)`*(flux_timeseries$`Soil Moisture (% day 0)`/100)

flux_pre <- flux_timeseries %>%
  filter(date < as.Date("2025-10-15"))

flux_post <- flux_timeseries %>%
  filter(date > as.Date("2025-10-14"))



library(glmmTMB)
library(ggeffects)
#check if soil moisture varies between habitat/bracken
model_moist <- glmmTMB(
  `Soil moisture (% wet mass)` ~ Bracken*Habitat  + (1|`Soil volume (m3)`),
  data = flux_pre,
  family = gaussian()
)
summary(model_moist)


ggboxplot(flux_pre, x = "Habitat", aes(y = `Soil moisture (% wet mass)`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  "Soil Moisture (%)") + 
  theme(
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


#model whether methane and carbon dioxide flux change due to habitat, veg, with soil volume in the core as a random effect - before rainfall


model_co2 <- glmmTMB(
  `CO2 flux (micromol CO2 per s per m2)` ~ Bracken*Habitat + (1|`Soil volume (m3)`),
  data = flux_pre,
  family = gaussian()
)
summary(model_co2)


model_ch4 <- glmmTMB(
  `CH4 flux (nmol CH4 per s per m2)` ~ Bracken*Habitat + (1|`Soil volume (m3)`),
  data = flux_pre,
  family = gaussian()
)
summary(model_ch4)

#model whether methane and carbon dioxide flux change due to habitat, veg, with soil volume in the core as a random effect - after rainfall


model_co2 <- glmmTMB(
  `CO2 flux (micromol CO2 per s per m2)` ~ Bracken*Habitat + (1|`Soil volume (m3)`),
  data = flux_post,
  family = gaussian()
)
summary(model_co2)


model_ch4 <- glmmTMB(
  `CH4 flux (nmol CH4 per s per m2)` ~ Bracken*Habitat + (1|`Soil volume (m3)`),
  data = flux_post,
  family = gaussian()
)
summary(model_ch4)

#### Question 2: Is there a difference in fluxes between the simulated rainfall intensities ? ----
#flux data after rainfall was added

#assign a new variable, Pre/Post rainfall, to compare fluxes pre rainfall to fluxes post

#flux_post <-  flux_timeseries %>%
#  filter(date > as.Date("2025-10-15"))

#is there a differnce in flux between rainfall intensities?

model_co2 <- glmmTMB(
  `CO2 flux (micromol CO2 per s per m2)` ~ Bracken*Habitat*`Rainfall Intensity (ml)`,
  data = flux_post,
  family = gaussian()
)
summary(model_co2)

model_ch4 <- glmmTMB(
  `CH4 flux (nmol CH4 per s per m2)` ~ Bracken*Habitat*`Rainfall Intensity (ml)`,
  data = flux_post,
  family = gaussian()
)
summary(model_ch4)

#see if there is a difference in moisture between rainfall treatments
model_moist <- glmmTMB(
  `Estimated soil moisture (% wet mass)` ~ `Rainfall Intensity (ml)`*Bracken*Habitat  + (1|`Soil volume (m3)`),
  data = flux_post,
  family = gaussian()
)
summary(model_moist)


#see if there is a difference in moisture between treatments
model_moist <- glmmTMB(
  `Soil moisture (% wet mass)` ~ Bracken*Habitat  + (1|`Soil volume (m3)`),
  data = flux_post,
  family = gaussian()
)
summary(model_moist)


ggboxplot(flux_post, x = "Habitat", aes(y = `Soil moisture (% wet mass)`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  "Soil Moisture (%)") + 
  theme(
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

#### Question 3: Is there a difference in fluxes before vs after extreme rainfall simulation? ----

#assign a new variable, Pre/Post rainfall, to compare fluxes pre rainfall to fluxes post
flux_timeseries <- flux_timeseries %>%
  mutate(`pre/post` = if_else(date < as.Date("2025-10-15"),
                              "pre",
                              "post"))


model_co2 <- glmmTMB(
  `CO2 flux (micromol CO2 per s per m2)` ~ `pre/post`*Bracken*Habitat,
  data = flux_timeseries,
  family = gaussian()
)
summary(model_co2)

model_ch4 <- glmmTMB(
  `CH4 flux (nmol CH4 per s per m2)` ~ `pre/post`*Bracken*Habitat,
  data = flux_timeseries,
  family = gaussian()
)
summary(model_ch4)

#is there a difference in mositure pre/post treatment?
model_co2 <- glmmTMB(
 `Estimated soil moisture (% wet mass)` ~ `pre/post`,
  data = flux_timeseries,
  family = gaussian()
)
summary(model_co2)

model_ch4 <- glmmTMB(
  `CH4 flux (nmol CH4 per s per m2)` ~ `Soil Moisture (% day 0)`,
  data = flux_timeseries,
  family = gaussian()
)
summary(model_ch4)

# plot the co2 flux over time

#group by date, Habitat, and Bracken, and calculate means
flux_avg <- flux_timeseries %>%
  group_by(date, Habitat, Bracken) %>%
  summarise(
    `Average CO2 flux` = mean(`CO2 flux (micromol CO2 per s per m2)`, na.rm = TRUE),
    `SD CO2 flux` = sd(`CO2 flux (micromol CO2 per s per m2)`, na.rm = TRUE),
    `Average CH4 flux` = mean(`CH4 flux (nmol CH4 per s per m2)`, na.rm = TRUE),
    `SD CH4 flux` = sd(`CH4 flux (nmol CH4 per s per m2)`, na.rm = TRUE),
    `Average Moisture` = mean(`Estimated soil moisture (% wet mass)`, na.rm = TRUE),
    `SD Moisture` = sd(`Estimated soil moisture (% wet mass)`, na.rm = TRUE),
    .groups = 'drop'
  )

# Create a combined grouping label
flux_avg <- flux_avg %>%
  mutate(Habitat_Bracken = paste(Habitat, Bracken, sep = " - "))


# plot the co2 flux over time
moisture_avg_tms <- ggplot(flux_avg, aes(x = date, y = `Average Moisture`, color = Habitat_Bracken)) +
  geom_line() +
  geom_point()  +
  geom_errorbar(
    aes(ymin = `Average Moisture` - `SD Moisture`, 
        ymax = `Average Moisture` + `SD Moisture`),
    width = 0.2
  ) +
  labs(
    title = "Average Soil Moisture Over Time by Habitat and Bracken Presence",
    x = "Date",
    y = expression("Estimated average soil moisture (% wet mass)"),
    color = "Group"
  ) + geom_vline(xintercept = as.numeric(as.Date("2025-10-14")), 
                 linetype = "dashed", color = "blue", linewidth = 0.5) +
  theme_minimal()
#show the plot
show(moisture_avg_tms)
#save our plot
ggsave(path = "Figures", paste0(Sys.Date(), "_mesocosm-moisture-timeseries.svg"), width = 10, height= 5, moisture_avg_tms)



#### Old code ----

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


#### model whether methane and carbon dioxide flux change due to habitat, veg, Soil moisture, rack, rainfall mix, rainfall volume ----
library(glmmTMB)
library(ggeffects)

#reorder the sites so they show up on the plot from west (LHS) to east (RHS)
results_df$Bracken <- factor(results_df$Bracken, levels = c("Present", "Absent"))
#reorder the habitats so they appear in the specified order
results_df$Habitat <- factor(results_df$Habitat, levels = c("Grassland", "Heathland"))

#this uses data from the entire experiment - should we instead split the data into before and after the rainfall application

model_co2 <- glmmTMB(
  `CO2 flux (micromol CO2 per s per m2)` ~ Bracken*Habitat + (1|`Rainfall volume (ml)`),
  data = results_df,
  family = gaussian(link = "logit")
)
summary(model_co2)


model_ch4 <- glmmTMB(
  `CH4 flux (nmol CH4 per s per m2)` ~ Bracken*Habitat + (1|`Rainfall volume (ml)`),
  data = results_df,
  family = gaussian(link = "logit")
)
summary(model_ch4)


# Helper function to extract model summary
extract_model_summary <- function(model, response_name) {
  coefs <- summary(model)$coefficients$cond
  data.frame(
    Response = response_name,
    Effect = rownames(coefs),
    Estimate = coefs[, "Estimate"],
    Std_Error = coefs[, "Std. Error"],
    z_value = coefs[, "z value"],
    p_value = coefs[, "Pr(>|z|)"]
  ) %>%
    mutate(Significant = ifelse(p_value < 0.05, "Yes", "No"))
}

# Extract summaries for all models
summary_model <- extract_model_summary(model_co2, "CO2 flux (micromol CO2 per s per m2)")
summary_model
