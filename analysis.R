# Load libraries
library(tidyverse)
library(here)
library(ggplot2)
library(lubridate)
library(ggpubr)
library(multcompView) #for significant difference letters

flux_data <- readr::read_delim(
  here::here("Data", "Gas Flux Measurements", "gas-flux-29-09-2025.csv"))

# Convert Time to POSIXct with milliseconds
flux_data <- flux_data %>%
  mutate(Time = as.POSIXct(Time, format = "%H:%M:%OS", tz = "UTC"))

start_time <- as.POSIXct("16:34:16.000", format = "%H:%M:%OS", tz = "UTC")
end_time   <- as.POSIXct("16:36:25.000", format = "%H:%M:%OS", tz = "UTC")

filtered_data <- flux_data %>%
  filter(Time >= start_time & Time <= end_time) 


ggplot(filtered_data, aes(x = Time, y = `[CO2]d_ppm`, group = 1)) +
  geom_line(color = "blue") +
  geom_point(color = "red", size = 1.5) +
  labs(
    title = "[CO2]d_ppm over Time",
    x = "Time",
    y = "[CO2]d_ppm"
  ) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  theme_minimal()

#### Linear regression to extract gradient ----

# Create numeric Time in seconds since start of measurement
filtered_data <- filtered_data %>%
  mutate(Time_sec = as.numeric(Time - min(Time)))
#run linear regression
model <- lm(`[CO2]d_ppm` ~ Time_sec, data = filtered_data)
#the gradient of the slope (units ppm per second)
summary(model)$coefficients["Time_sec", "Estimate"]

#### Function to automatically trim flux data based on time segments

library(tidyverse)
library(lubridate)

# --- Load flux data ---
flux_data <- readr::read_delim(
  here::here("Data", "Gas Flux Measurements", "gas-flux-29-09-2025.csv"))

#remove commas from time row
flux_data$Time <- gsub(",", "", flux_data$Time)

# Combine Date and Time safely (for dd/mm/yyyy format with fractional seconds)
flux_data <- flux_data %>%
  mutate(
    DateTime = as.POSIXct(paste(Date, Time), format = "%d/%m/%Y %H:%M:%OS", tz = "UTC")
  )
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


# --- Initialize list to store regression results ---
results_list <- list()

# --- Loop through each time range ---
for (i in 1:nrow(time_ranges)) {
  sample_id <- time_ranges$Mesocosm[i]
  start_dt <- time_ranges$start_datetime[i]
  end_dt <- time_ranges$end_datetime[i]
  
  # Subset flux data for this time window
  data_chunk <- flux_data %>%
    filter(DateTime >= start_dt & DateTime <= end_dt) %>%
    mutate(Time_sec = as.numeric(DateTime - min(DateTime)))  # time relative to chunk start
  
  # Only run regression if enough data
  if (nrow(data_chunk) >= 2) {
    model <- lm(`[CO2]d_ppm` ~ Time_sec, data = data_chunk)
    slope <- coef(model)[["Time_sec"]]
  } else {
    slope <- NA
  }
  
  # Store the result
  results_list[[i]] <- tibble(
    sample_id = sample_id,
    start_time = start_dt,
    end_time = end_dt,
    slope_ppm_per_sec = slope
  )
}

# --- Combine results into one data frame ---
results_df <- bind_rows(results_list)

# --- View or export ---
print(results_df)


#append environmental factors 
factors <- readr::read_delim(
  here::here("Data", "Gas Flux Measurements", "mesocosm_id_factors.csv"))

# Left join to keep all rows from results_df and add matching Habitat and Vegetation
results_df <- results_df %>%
  left_join(factors, by = "sample_id")


#reorder the sites so they show up on the plot from west (LHS) to east (RHS)
results_df$Bracken <- factor(results_df$Bracken, levels = c("Present", "Absent"))
#reorder the habitats so they appear in the specified order
results_df$Habitat <- factor(results_df$Habitat, levels = c("Grassland", "Heathland"))

#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
co2_bxp <- ggboxplot(results_df, x = "Habitat", aes(y = `slope_ppm_per_sec`), color = "Bracken", lwd = 0.75)  + 
  labs(x = "Habitat",
       y = "CO2 rate (ppm per sec") + theme(
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
ggsave(path = "Figures", paste0(Sys.Date(), "_co2-ppm-per-sec.svg"), width = 10, height= 5, co2_bxp)

#Type 1 two-way anova using data from all sites
anova <- aov(results_df$slope_ppm_per_sec ~ results_df$Bracken*results_df$Habitat)
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
# write_csv(results_df, "co2_slopes_output.csv")


  