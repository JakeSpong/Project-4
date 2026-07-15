# Load libraries
library(ggplot2)
library(tidyr)
library(readxl) #for exporting tabs from xlsx to csv
library(tibble) #for renaming rows
library(readr) #for reading in many files
library(dplyr) #for making factor columns
library(stringr) #for making factor columns
library(glmmTMB) #for modelling results
library(sf) #for maniuplating sample coordinates
library(emmeans) # post-hoc comparisons
library(DHARMa)      # model diagnostics
library(multcomp) #for significance letters
library(purrr) # for mapping when running GLMMs
library(coin) #for checking effect sizes of Wilcox test
library(rstatix) #for wilcox test
library(svglite)

library(tidyverse)
library(here)

#library(lubridate)
library(ggpubr)
library(multcompView) #for significant difference letters
library(scales)
library(rlang)
#### calculate soil moisture content at end of mesocosm experiment NO LONGER NEEDED ----

#load in the moisture data, format data correctly
field_moisture <- read_csv("Data/Soil Moisture (field conditions).csv")
#field condition moisture
initial_moisture <- field_moisture[, c(1, 7)]
#moisture data
moisture_data <- read_csv("Data/mesocosm moisture - Soil mass change as % of day 0.csv")
#add %change in mass (i.e. moisture)
initial_moisture <- initial_moisture %>%
  left_join(moisture_data %>% select(1, 21), by = "Sample ID")
#calculate soil moisture as % wet mass for the day of sampling for enzymes
initial_moisture$`end Soil moisture (% wet mass)` <- 100*(1-((100/initial_moisture$`Soil Moisture (% day 0) destructive sampling (13 or 14/11/2025)`)*(1-(initial_moisture$`Soil moisture (% wet mass)`/100))))
#add the -M suffix
initial_moisture$`Sample ID` <- ifelse(!grepl("-M$", initial_moisture$`Sample ID`),
                                       paste0(initial_moisture$`Sample ID`, "-M"),
                                       initial_moisture$`Sample ID`)


#load covariates
covs <- read_csv("Data/Enzyme_variables.csv")

#add %change in mass (i.e. moisture)
covs <- covs %>%
  left_join(initial_moisture %>% select(1, 4), by = "Sample ID")

#put all soil moisture values in one column
covs <- covs %>%
  mutate(`Soil moisture (% wet mass)` = coalesce(`Soil moisture (% wet mass)`, `end Soil moisture (% wet mass)`))
#delete the unnecessary columns
covs <- covs[, 1:(ncol(covs)-2)]
#fill out dry soil mass row
covs$`Dry soil mass (g)` <- (covs$`Wet soil mass (g)`)*(1-(covs$`Soil moisture (% wet mass)`/100))

#save the updated covs file
write.csv(covs, "Data/Enzyme_variables.csv", row.names = FALSE)





 
#load in microplate, seprate tabs into csvs to analyse NO LONGER NEEDED ----

# Path to your Excel file
#excel_file <- "Data/Microplate Reader Data/26-06-2026-labelled_Microplate.xlsx"
excel_file <- "Data/Microplate Reader Data/01-07-2026-labelled_Microplate.xlsx"
# Get all sheet (tab) names
sheet_names <- excel_sheets(excel_file)

# Folder to save the CSVs into
output_dir <- "Data/Microplate Reader Data"

# Create the folder if it doesn't already exist
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
# Loop through each sheet and write it out as its own CSV
for (sheet in sheet_names) {
  df <- read_excel(excel_file, sheet = sheet)
  # clean the sheet name a bit so it's a safe filename (optional)
  safe_name <- gsub("[^A-Za-z0-9_-]", "_", sheet)
  write.csv(df, file = file.path(output_dir, paste0(safe_name, ".csv")), row.names = FALSE)
}
#tabs then sorted into folders for easy data analysis



## get the standard curve ---- NO LONGER NEEDED ----

# the standards dataframe
#standards <- read_csv("Data/Microplate Reader Data/26-06-2026 plates gain700/Standard_gain700.csv")
#data for column 4 and column 8 removed as column 4 is the same as column 3, and column 8 is lower than column 7...likely a pipetting mistake
standards <- read_csv("Data/Microplate Reader Data/01-07-2026 plates gain755/Standard_gain755.csv")
#just the data valyes
std_raw <- standards[13:20, 2:10]
# Ensure every column in a data frame is numeric (not character/factor)
std_raw[] <- lapply(std_raw, function(x) as.numeric(as.character(x)))

#the actual concentrations of the standards, in micromols
concentrations <- c(0, 0.0001728, 0.000432, 0.000864, 0.0012384,
                    0.0016416, 0.0020736, 0.0024768, 0.003312)
#set as matrix so we can have row names
std_raw <- as.matrix(std_raw)
#setr row and column names
rownames(std_raw) <- LETTERS[1:8]
colnames(std_raw) <- concentrations

#subtract column 1 (blank) from every column, row-wise -
blank <- std_raw[, 1]
std_corrected <- sweep(std_raw, 1, blank, FUN = "-")

#back to a df from a matrix
std_table <- as.data.frame(std_corrected)

#Per-row regression through the origin (LINEST with FALSE)
slopes <- apply(std_table, 1, function(y) {
  fit <- lm(y ~ 0 + concentrations)
  coef(fit)[["concentrations"]]
})
#name the columns so we know which slope belongs to which row of the microplate
names(slopes) <- rownames(std_raw)
#mean slope
mean_slope <- mean(slopes)

#plot the regression slopes

# reshape data into long format for plotting
plot_data <- std_table
plot_data$replicate <- rownames(std_table)
plot_long <- pivot_longer(plot_data, -replicate,
                          names_to = "concentration", values_to = "signal")
plot_long$concentration <- as.numeric(plot_long$concentration)

# build a data frame of fitted lines through the origin
line_data <- data.frame(
  replicate = names(slopes),
  slope = as.numeric(slopes)
)

ggplot(plot_long, aes(x = concentration, y = signal, color = replicate)) +
  geom_point(size = 2) +
  geom_abline(data = line_data,
              aes(intercept = 0, slope = slope, color = replicate)) +
  labs(x = "Concentration (umol)", y = "Absorbance",
       title = "Standard curves with fitted regression lines",
       color = "Replicate") +
  theme_minimal()


#export the slope
gain_value <- 755   # <- replace with gain value as per the standards being analysed

# Build a one-row data frame with the requested headers
new_row <- data.frame(
  Gain = gain_value,
  `Mean Slope` = mean_slope,
  check.names = FALSE
)

# Write to CSV
#write.csv(output_df, "Processed Data/mean_slope_output.csv", row.names = FALSE)
# Append to the existing CSV (no header, since it's already there) - ensure this file isn't open!
write.table(new_row, "Processed Data/mean_slope_output.csv",
            sep = ",", append = TRUE,
            row.names = FALSE, col.names = FALSE)

#### analyse the microplate data ---- NO LONGER NEEDED ----

#import slopes from standard curves
stnd_slopes <- read_csv("Processed Data/mean_slope_output.csv")

# Function to analyse a single microplate CSV
analyse_plate <- function(file_path, slope_value, gain_value) {
  
  raw <- read_csv(file_path, show_col_types = FALSE)
  
  # just the data we want
  df <- raw[12:20, 1:13]
  
  # name the columns from row 1, then drop that row
  colnames(df) <- as.character(df[1, ])
  df <- df[-1, ]
  
  # name the id column, set to row name
  names(df)[1] <- "row_id"
  df <- column_to_rownames(df, var = "row_id")
  
  # ensure numbers are coded as numbers (row names preserved via column_to_rownames)
  df <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))), check.names = FALSE)
  
  # ---- outlier detection ----
  # top 4 rows are sample (Extract), bottom 4 rows are control
  group_size <- 4
  n_groups <- ceiling(nrow(df) / group_size)
  group_id <- rep(1:n_groups, each = group_size, length.out = nrow(df))
  
  # quantile(..., type = 7) is R's default and matches Excel's QUARTILE.INC
  quartile_df <- data.frame(
    SampleID = character(), Group = integer(),
    Q1 = numeric(), Q3 = numeric(), IQR = numeric(),
    LowerBound = numeric(), UpperBound = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (col in colnames(df)) {
    for (g in unique(group_id)) {
      vals <- df[group_id == g, col]
      q1 <- quantile(vals, 0.25, type = 7, na.rm = TRUE)
      q3 <- quantile(vals, 0.75, type = 7, na.rm = TRUE)
      iqr <- q3 - q1
      quartile_df <- rbind(quartile_df, data.frame(
        SampleID = col, Group = g, Q1 = q1, Q3 = q3, IQR = iqr,
        LowerBound = q1 - 1.5 * iqr, UpperBound = q3 + 1.5 * iqr
      ))
    }
  }
  rownames(quartile_df) <- NULL
  # rename so we can easily see which data is the Sample (Extract) vs Control
  quartile_df$Group <- factor(quartile_df$Group, levels = c(1, 2), labels = c("Extract", "Control"))
  
  # remove any outliers, using the same Extract/Control grouping
  group_labels <- c(`1` = "Extract", `2` = "Control")
  df_no_outliers <- df
  
  outlier_log <- data.frame(
    RowLabel = character(), SampleID = character(), Group = character(),
    Value = numeric(), LowerBound = numeric(), UpperBound = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (col in colnames(df)) {
    for (g in unique(group_id)) {
      rows_in_group <- which(group_id == g)
      g_label <- group_labels[[as.character(g)]]
      bounds <- quartile_df[quartile_df$SampleID == col & quartile_df$Group == g_label, ]
      lower <- bounds$LowerBound
      upper <- bounds$UpperBound
      
      vals <- df_no_outliers[[col]][rows_in_group]
      is_outlier <- vals < lower | vals > upper
      is_outlier[is.na(is_outlier)] <- FALSE
      
      if (any(is_outlier)) {
        outlier_log <- rbind(outlier_log, data.frame(
          RowLabel = rownames(df)[rows_in_group][is_outlier],
          SampleID = col, Group = g_label,
          Value = vals[is_outlier],
          LowerBound = lower, UpperBound = upper
        ))
      }
      
      vals[is_outlier] <- NA
      df_no_outliers[[col]][rows_in_group] <- vals
    }
  }
  rownames(outlier_log) <- NULL
  
  # ---- average sample & control, subtract, convert to concentration ----
  plate_result <- data.frame(
    SampleID = character(), SampleMean = numeric(), ControlMean = numeric(),
    Difference = numeric(), Concentration = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (col in colnames(df_no_outliers)) {
    sample_mean  <- mean(df_no_outliers[[col]][group_id == 1], na.rm = TRUE)
    control_mean <- mean(df_no_outliers[[col]][group_id == 2], na.rm = TRUE)
    difference   <- sample_mean - control_mean
    concentration <- difference / slope_value
    
    plate_result <- rbind(plate_result, data.frame(
      SampleID = col, SampleMean = sample_mean, ControlMean = control_mean,
      Difference = difference, Concentration = concentration
    ))
  }
  
  rownames(plate_result) <- NULL
  plate_result$SourceFile <- basename(file_path)
  plate_result$Gain <- gain_value
  plate_result$SlopeUsed <- slope_value
  
  outlier_log$SourceFile <- basename(file_path)
  
  list(plate_result = plate_result, outlier_log = outlier_log)
 
}

# Find all plate CSVs to process (excluding standard curve files)
plate_files <- list.files(
  "Data/Microplate Reader Data",
  pattern = "\\.csv$",
  recursive = TRUE,
  full.names = TRUE
)

# drop any files that are standard curve files, not sample plates
plate_files <- plate_files[!grepl("Standard", plate_files, ignore.case = TRUE)]

#  Loop over every plate file, matching gain from its folder name, and combine all results into one table
#dataframe to store results
results_df <- data.frame()
#dataframe to store outliers
all_outliers_df <- data.frame()

#run the microplate analyser function, generating an output table of results, and a df containing the outliers that were removed 
for (file_path in plate_files) {
  
  # extract gain number from the folder name, e.g. ".../gain782/..." -> 782
  gain_match <- regmatches(file_path, regexpr("gain(\\d+)", file_path))
  gain_value <- as.numeric(gsub("gain", "", gain_match))
  
  # look up the matching slope for that gain
  slope_row <- stnd_slopes[stnd_slopes$Gain == gain_value, ]
  slope_value <- slope_row[[2]]  # "Mean Slope" column
  
  plate_output <- analyse_plate(file_path, slope_value, gain_value)
  results_df <- rbind(results_df, plate_output$plate_result)
  all_outliers_df <- rbind(all_outliers_df, plate_output$outlier_log)
}

#remove rows where the microplate wells were blank
results_df <- results_df[nchar(results_df$SampleID) > 2, ]

#reformat so we know which enzyme which data refers to
results_df$SampleType <- gsub("\\.csv$", "", results_df$SourceFile)  # remove ".csv"
results_df$SampleType <- gsub("[0-9]", "", results_df$SampleType)     # remove digits

#check all SampleIDs in the results_df are consistently formatted with their -M suffixes
results_df$SampleID <- gsub("(?<!-)\\s?M$", "-M", results_df$SampleID, perl = TRUE)
#ensure Sample ID is named correctly
colnames(results_df)[1] <- "Sample ID"


#add covariates to calculate EEA per hourr per g dry weight
covs <- read_csv("Data/Enzyme_variables.csv")
#join covariates to results_df in a new df - we want the date samples frozen, incubation time (hours), dry soil mass, along with volumes of liquids
results <- results_df %>%
  left_join(covs %>% select(1, 4, 6, 8, 9, 10, 11), by = "Sample ID")

#calculate dilution factor
results$`Dilution Factor` <- results$`Phosphate buffer volume added to dialysis tubes (microliters)`/results$`Volume sample added to microplate well (microliters)`
#calculate concentration factor
results$`Concentration Factor` <- results$`Phosphate buffer volume added to dialysis tubes (microliters)`/results$`Initial extract volume (microliters)`

#multiply by dilution factor
results$`moles enzyme in extract (micromol)` <- results$Concentration*results$`Dilution Factor`
#multiply by concentration factor
results$`correction for concentration step (sausages)` <- results$`moles enzyme in extract (micromol)`*results$`Concentration Factor`
#calculate extracellular enzyme activity per hour her g dry weight soil
results$`EEA per hour per g dry soil` <- ((results$`correction for concentration step (sausages)`/results$`Dry soil mass (g)`)*1000)/results$`Incubation time (hours)`

#add columns for the various factors
results <- results %>%
  mutate(
    #if M is present, it is from end of mesocosm study.  If no M, from field conditions
    Condition = if_else(str_detect(`Sample ID`, "-M$"), "Mesocosm", "Field"),
    
    # extract the 2-letter code (HB, HH, GG, GB) regardless of whether -M is present
    HabitatCode = str_extract(`Sample ID`, "[A-Z]{2}(?=(-M)?$)"),
    #if first letter of two-leter code H, it's heathland; if G, grassland
    Habitat = case_when(
      substr(HabitatCode, 1, 1) == "H" ~ "Heathland",
      substr(HabitatCode, 1, 1) == "G" ~ "Grassland",
      TRUE ~ NA_character_
    ),
    #if second letter of two-letter code is B, bracken is present, if G or H, bracken is absent 
    Bracken = case_when(
      substr(HabitatCode, 2, 2) == "B" ~ "Present",
      substr(HabitatCode, 2, 2) %in% c("G", "H") ~ "Absent",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-HabitatCode)  # drop the helper column, keep only the three requested


#save the data so we don't have to run the function again
write.csv(results, "Processed Data/EEA (nmol per h per g dry weight).csv", row.names = FALSE)


#### generate boxplot of EEAs under inital conditons ----

#load in the EEA data
results <- read_csv("Processed Data/EEA (nmol per h per g dry weight).csv")
#load in explanatory covaraite data
expcovs <- read_csv("Data/Explanatory covariates.csv")

#add explanatory covariate data to the results df, matching regardless of -M suffix
results <- results %>%
  mutate(join_id = str_remove(`Sample ID`, "-M$")) %>%
  left_join(expcovs, by = c("join_id" = "Sample ID")) %>%
  dplyr::select(-join_id)


# Split by Condition first
mesocosm_df <- results %>% filter(Condition == "Mesocosm")
field_df    <- results %>% filter(Condition == "Field")

# Strip the "-M" suffix so IDs are directly comparable to Field
mesocosm_df <- mesocosm_df %>%
  mutate(BaseID = gsub("-M$", "", `Sample ID`))

#generate a column to compare sample IDs from field and end of mesocosm
field_df <- field_df %>%
  mutate(BaseID = `Sample ID`)  # Field IDs already have no -M suffix

# Keep only IDs present in both sets
common_ids <- intersect(mesocosm_df$BaseID, field_df$BaseID)

#subset of mesocosm end samples to match those from initial conditions
mesocosm_subset <- mesocosm_df %>% filter(BaseID %in% common_ids)
#all the field samples
field_df    <- field_df %>% filter(BaseID %in% common_ids)
#all mesocosm end samples
mesocosm_all <- mesocosm_df

sampletype_labels <- c(
  "Gluc" = "β-glucosidase",
  "NAG" = "β-1,4-N-acetyl-glucosaminidase",
  "Phos" = "Acid phosphatase",
  "Xylo" = "β-xylosidase"
)

#plot the intial EEAs
initial <- ggplot(field_df, aes(x = Habitat, y = `EEA per hour per g dry soil`, fill = Bracken)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), 
             aes(color = Bracken), alpha = 0.6, size = 1.5) +
  facet_wrap(~ SampleType, labeller = labeller(SampleType = sampletype_labels)) +
  scale_fill_manual(values = c("Absent" = "sienna", "Present" = "limegreen")) +
  scale_color_manual(values = c("Absent" = "sienna", "Present" = "limegreen")) +
  coord_cartesian(ylim = c(0, 80)) +
  labs(x = "Habitat", 
       #tildes used to get spaces right
       y = expression("Extracellular enzyme activity (nmol" ~~ hr^-1 ~ "g"^-1 ~ "dry weight)")) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.background = element_rect(fill = "grey85", color = "black"),
        panel.grid = element_blank())
show(initial)
#save the figure
ggsave("Figures/initial_EEA_per_hour_by_habitat_bracken.svg", plot = initial, 
       width = 10, height = 7, dpi = 300)


#plot subset of final EEAs that match the initial (for a fair comparison of change)
final_subset <- ggplot(mesocosm_subset, aes(x = Habitat, y = `EEA per hour per g dry soil`, fill = Bracken)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), 
             aes(color = Bracken), alpha = 0.6, size = 1.5) +
  facet_wrap(~ SampleType, labeller = labeller(SampleType = sampletype_labels)) +
  scale_fill_manual(values = c("Absent" = "sienna", "Present" = "limegreen")) +
  scale_color_manual(values = c("Absent" = "sienna", "Present" = "limegreen")) +
  coord_cartesian(ylim = c(0, 80)) +
  labs(x = "Habitat", 
       #tildes used to get spaces right
       y = expression("Extracellular enzyme activity (nmol" ~~ hr^-1 ~ "g"^-1 ~ "dry weight)")) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.background = element_rect(fill = "grey85", color = "black"),
        panel.grid = element_blank())

show(final_subset)
#save the figure
ggsave("Figures/final_subset_EEA_per_hour_by_habitat_bracken.svg", plot = final_subset, 
       width = 10, height = 7, dpi = 300)



#plot final EEA from all reps (useful for telling differences between rainfall treatments)
final_all <- ggplot(mesocosm_all, aes(x = Habitat, y = `EEA per hour per g dry soil`, fill = Bracken)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), 
             aes(color = Bracken), alpha = 0.6, size = 1.5) +
  facet_wrap(~ SampleType, labeller = labeller(SampleType = sampletype_labels)) +
  scale_fill_manual(values = c("Absent" = "sienna", "Present" = "limegreen")) +
  scale_color_manual(values = c("Absent" = "sienna", "Present" = "limegreen")) +
  coord_cartesian(ylim = c(0, 80)) +
  labs(x = "Habitat", 
       #tildes used to get spaces right
       y = expression("Extracellular enzyme activity (nmol" ~~ hr^-1 ~ "g"^-1 ~ "dry weight)")) +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.background = element_rect(fill = "grey85", color = "black"),
        panel.grid = element_blank())

show(final_all)
#save the figure
ggsave("Figures/final_all_EEA_per_hour_by_habitat_bracken.svg", plot = final_all, 
       width = 10, height = 7, dpi = 300)

#run models to see if differences are significant NEEDS PREVIOUS TAB TO WORK----

#initial conditions

# convert longitude, latitude to sf object for initial field samples
sf_data <- st_as_sf(field_df, coords = c("LongitudeE", "LatitudeN"), crs = 4326)
# project to UTM (example: zone depends on your location)
sf_data <- st_transform(sf_data, crs = 32630)  # UK example
# extract projected coords
coords <- st_coordinates(sf_data)
#coordinates in m in the Universal Transverse Meractor cartesian coordinate system
field_df$x <- coords[,1]
field_df$y <- coords[,2]

#first split out data for each enzyme
sample_types <- unique(field_df$SampleType)
split_dfs <- field_df %>%
  group_by(SampleType) %>%
  group_split() %>%
  setNames(sample_types)

#check the data distripution for each enzyme - they all look right skewed!
hist(split_dfs[["Gluc"]]$`EEA per hour per g dry soil`)

#run glmms, using log function to account for right skewed distribution
models <- split_dfs %>%
  map(~ glmmTMB(`EEA per hour per g dry soil` ~ Habitat*Bracken, 
                family = Gamma(link = "log"), 
                data = .x))


#analyse each enzyme in turn
model <- models$Xylo
#look at the summary
summary(model)

# ── Confirm no spatial autocorrelation in residuals ───────────────────────────
sim_res <- simulateResiduals(model)
testSpatialAutocorrelation(sim_res, x = split_dfs$Xylo$x, y = split_dfs$Xylo$y)
sim_res <- simulateResiduals(model)
#check the plot for residuals, residuals vs predicted
plot(sim_res)
#check dispersion
testDispersion(sim_res)

# ── Type II significance of fixed effects ────────────────────────────────────
car::Anova(model, type = "II")

# ── Post-hoc comparisons ──────────────────────────────────────────────────────
emm <- emmeans(model, ~ Habitat*Bracken, type = "response")
# ── Pairwise contrasts to confirm which differences are significant ────────────
pairs(emm, adjust = "tukey")
#get signficance letters
cld_result <- cld(emm, Letters = letters, adjust = "tukey")
# 5. View results
print(cld_result)

# ── Or more targeted: effect of bracken within each habitat ───────────────────
emm_hab <- emmeans(model, ~ Bracken | Habitat, type = "response")
pairs(emm_hab, adjust = "tukey")

#subset of final conditions to match those samples used under initial conditions

# convert longitude, latitude to sf object for initial field samples
sf_data <- st_as_sf(mesocosm_subset, coords = c("LongitudeE", "LatitudeN"), crs = 4326)
# project to UTM (example: zone depends on your location)
sf_data <- st_transform(sf_data, crs = 32630)  # UK example
# extract projected coords
coords <- st_coordinates(sf_data)
#coordinates in m in the Universal Transverse Meractor cartesian coordinate system
mesocosm_subset$x <- coords[,1]
mesocosm_subset$y <- coords[,2]

#first split out data for each enzyme
sample_types <- unique(mesocosm_subset$SampleType)
split_dfs <- mesocosm_subset %>%
  group_by(SampleType) %>%
  group_split() %>%
  setNames(sample_types)

#check the data distripution for each enzyme - they all look right skewed!
hist(split_dfs[["Xylo"]]$`EEA per hour per g dry soil`)

#run glmms, using log function to account for right skewed distribution
models <- split_dfs %>%
  map(~ glmmTMB(`EEA per hour per g dry soil` ~ Habitat*Bracken, 
                family = Gamma(link = "log"), 
                data = .x))

#analyse each enzyme in turn
model <- models$Xylo
#look at the summary
summary(model)

# ── Confirm no spatial autocorrelation in residuals ───────────────────────────
sim_res <- simulateResiduals(model)
#change the enzymes when running the mdodels
testSpatialAutocorrelation(sim_res, x = split_dfs$Xylo$x, y = split_dfs$Xylo$y)
sim_res <- simulateResiduals(model)
#check the plot for residuals, residuals vs predicted
plot(sim_res)
#check dispersion
testDispersion(sim_res)

# ── Type II significance of fixed effects ────────────────────────────────────
car::Anova(model, type = "II")

# ── Post-hoc comparisons ──────────────────────────────────────────────────────
emm <- emmeans(model, ~ Habitat, type = "response")
# ── Pairwise contrasts to confirm which differences are significant ────────────
pairs(emm, adjust = "tukey")
#get signficance letters
cld_result <- cld(emm, Letters = letters, adjust = "tukey")
# 5. View results
print(cld_result)

# ── Or more targeted: effect of bracken within each habitat ───────────────────
emm_hab <- emmeans(model, ~ Bracken | Habitat, type = "response")
pairs(emm_hab, adjust = "tukey")


#allf final conditions to match those samples used under initial conditions

# convert longitude, latitude to sf object for initial field samples
sf_data <- st_as_sf(mesocosm_all, coords = c("LongitudeE", "LatitudeN"), crs = 4326)
# project to UTM (example: zone depends on your location)
sf_data <- st_transform(sf_data, crs = 32630)  # UK example
# extract projected coords
coords <- st_coordinates(sf_data)
#coordinates in m in the Universal Transverse Meractor cartesian coordinate system
mesocosm_all$x <- coords[,1]
mesocosm_all$y <- coords[,2]

#first split out data for each enzyme
sample_types <- unique(mesocosm_all$SampleType)
split_dfs <- mesocosm_all %>%
  group_by(SampleType) %>%
  group_split() %>%
  setNames(sample_types)

#check the data distripution for each enzyme - they all look right skewed!
hist(split_dfs[["Xylo"]]$`EEA per hour per g dry soil`)

#run glmms, using log function to account for right skewed distribution
models <- split_dfs %>%
  map(~ glmmTMB(`EEA per hour per g dry soil` ~ Habitat*Bracken, 
                family = Gamma(link = "log"), 
                data = .x))

#analyse each enzyme in turn
model <- models$NAG
#look at the summary
summary(model)

# ── Confirm no spatial autocorrelation in residuals ───────────────────────────
sim_res <- simulateResiduals(model)
testSpatialAutocorrelation(sim_res, x = split_dfs$NAG$x, y = split_dfs$NAG$y)
sim_res <- simulateResiduals(model)
#check the plot for residuals, residuals vs predicted
plot(sim_res)
#check dispersion
testDispersion(sim_res)

# ── Type II significance of fixed effects ────────────────────────────────────
car::Anova(model, type = "II")

# ── Post-hoc comparisons ──────────────────────────────────────────────────────
emm <- emmeans(model, ~ Habitat, type = "response")
# ── Pairwise contrasts to confirm which differences are significant ────────────
pairs(emm, adjust = "tukey")
#get signficance letters
cld_result <- cld(emm, Letters = letters, adjust = "tukey")
# 5. View results
print(cld_result)

# ── Or more targeted: effect of bracken within each habitat ───────────────────
emm_hab <- emmeans(model, ~ Bracken | Habitat, type = "response")
pairs(emm_hab, adjust = "tukey")

# run wilcox test to compare EEAs from field to end of mesocosm ----

#load in the EEA data
results <- read_csv("Processed Data/EEA (nmol per h per g dry weight).csv")
#load in explanatory covaraite data
expcovs <- read_csv("Data/Explanatory covariates.csv")

#add explanatory covariate data to the results df, matching regardless of -M suffix
results <- results %>%
  mutate(join_id = str_remove(`Sample ID`, "-M$")) %>%
  left_join(expcovs, by = c("join_id" = "Sample ID")) %>%
  dplyr::select(-join_id)

#define samples with -M suffix as being from end of experiment, those without as from the start of the experiment
results <- results %>%
  mutate(
    BaseID = str_remove(`Sample ID`, "-M$"),      # strip the -M suffix
    TimePoint = if_else(str_detect(`Sample ID`, "-M$"), "End", "Start")
  )

# Pivot wider, carrying Habitat and Bracken through as identifying columns
# (they should be constant within a BaseID, so no information is lost)
df_wide <- results %>%
  pivot_wider(
    id_cols = c(BaseID, SampleType, Habitat, Bracken),
    names_from = TimePoint,
    values_from = `EEA per hour per g dry soil`
  ) %>%
  filter(!is.na(Start) & !is.na(End))

#check for normality in the differences
normality_check <- df_wide %>%
  group_by(SampleType, Habitat, Bracken) %>%
  summarise(
    n_pairs = n(),
    shapiro_p = shapiro.test(End - Start)$p.value,
    .groups = "drop"
  )
#3 not normal, so use wilcox for consistency of tests
normality_check

#pivot back to long
df_long <- df_wide %>%
  pivot_longer(
    cols = c(Start, End),
    names_to = "TimePoint",
    values_to = "EEA per hour per g dry soil"
  ) %>%
  mutate(TimePoint = factor(TimePoint, levels = c("Start", "End")))
#run the wilcox tests
wilcox_results <- df_long %>%
  group_by(SampleType, Habitat, Bracken) %>%
  rstatix::wilcox_test(`EEA per hour per g dry soil` ~ TimePoint, paired = TRUE) %>%
  group_by(SampleType) %>%
  rstatix::adjust_pvalue(method = "BH") %>%
  rstatix::add_significance()
#view results
wilcox_results
#check effect sizes
effect_sizes <- df_long %>%
  group_by(SampleType, Habitat, Bracken) %>%
  rstatix::wilcox_effsize(`EEA per hour per g dry soil` ~ TimePoint, paired = TRUE)
#view effect sizes
effect_sizes




# Back to long format for rstatix (Start/End as levels of TimePoint again)
df_wide_long <- df_wide %>%
  pivot_longer(cols = c(Start, End), names_to = "TimePoint", values_to = "Activity") %>%
  mutate(TimePoint = factor(TimePoint, levels = c("Start", "End")))

#check normality of differences for each enzyme
results_wide %>%
  group_by(SampleType) %>%
  summarise(shapiro_p = shapiro.test(End - Start)$p.value)

#differences are not normal, so use a wilcox test, not t test
df_wide_long <- results_wide %>%
  pivot_longer(cols = c(Start, End), names_to = "TimePoint", values_to = "Activity")

#compute the wilcox test
df_wide_long %>%
  group_by(SampleType) %>%
  rstatix::wilcox_test(Activity ~ TimePoint, paired = TRUE) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()
#check effect size too
df_wide_long %>%
  group_by(SampleType) %>%
  wilcox_effsize(Activity ~ TimePoint, paired = TRUE)


ggplot(df_wide_long, aes(x = TimePoint, y = Activity)) +
  geom_boxplot(aes(fill = TimePoint), width = 0.4, alpha = 0.5, outlier.shape = NA) +
  geom_line(aes(group = BaseID), color = "grey50", alpha = 0.5) +
  geom_point(aes(color = TimePoint), size = 2, alpha = 0.7) +
  facet_wrap(~ SampleType, scales = "free_y") +
  scale_fill_manual(values = c("Start" = "#4C72B0", "End" = "#DD8452")) +
  scale_color_manual(values = c("Start" = "#4C72B0", "End" = "#DD8452")) +
  labs(
    title = "Paired enzyme activity: Start vs End",
    x = NULL,
    y = "Activity"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

#visualise
sample_comparisons <- ggplot(df_long, aes(x = TimePoint, y = `EEA per hour per g dry soil`)) +
  geom_line(aes(group = BaseID), color = "grey60", alpha = 0.5) +
  geom_point(aes(color = Bracken), size = 2, alpha = 0.8) +
  facet_grid(SampleType ~ Habitat + Bracken, scales = "free_y") +
  scale_color_manual(values = c("Absent" = "sienna", "Present" = "limegreen")) +
  labs(
    x = NULL,
    y = "EEA per hour per g dry soil",
    color = "Bracken"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "top",
    strip.background = element_rect(fill = "grey90"),
    strip.text = element_text(face = "bold", size = 9),
    panel.grid.minor = element_blank()
  )
show(sample_comparisons)
#save the figure
ggsave("Figures/EEA_sample-start-vs-end.svg", plot = sample_comparisons, 
       width = 10, height = 7, dpi = 300)



# do EEAs in final all differ with rainfall treatment? ----

#load in the EEA data
results <- read_csv("Processed Data/EEA (nmol per h per g dry weight).csv")
#load in explanatory covaraite data
expcovs <- read_csv("Data/Explanatory covariates.csv")

#add explanatory covariate data to the results df, matching regardless of -M suffix
results <- results %>%
  mutate(join_id = str_remove(`Sample ID`, "-M$")) %>%
  left_join(expcovs, by = c("join_id" = "Sample ID")) %>%
  dplyr::select(-join_id)

# Split by Condition first
mesocosm_df <- results %>% filter(Condition == "Mesocosm")

# convert longitude, latitude to sf object for initial field samples
sf_data <- st_as_sf(mesocosm_df, coords = c("LongitudeE", "LatitudeN"), crs = 4326)
# project to UTM (example: zone depends on your location)
sf_data <- st_transform(sf_data, crs = 32630)  # UK example
# extract projected coords
coords <- st_coordinates(sf_data)
#coordinates in m in the Universal Transverse Meractor cartesian coordinate system
mesocosm_df$x <- coords[,1]
mesocosm_df$y <- coords[,2]

#subset out for each enzyme
#first split out data for each enzyme
sample_types <- unique(mesocosm_df$SampleType)
split_dfs <- mesocosm_df %>%
  group_by(SampleType) %>%
  group_split() %>%
  setNames(sample_types)

#check the data distripution for each enzyme - they all look right skewed!
hist(split_dfs[["Xylo"]]$`EEA per hour per g dry soil`)

#model EEA as a result of rainfall treatment with habitat*bracken

#run glmms, using log function to account for right skewed distribution


models <- split_dfs %>%
  map(~ glmmTMB(`EEA per hour per g dry soil` ~ Habitat*Bracken*`Rainfall Intensity (ml)`, 
                family = Gamma(link = "log"), 
                data = .x))

model <- models$Gluc
#look at the summary
summary(model)

# Get slopes of Rainfall Intensity within each Habitat x Bracken combination,
# for every model in the list
trends_list <- models %>%
  map(~ emtrends(.x, ~ Habitat * Bracken, var = "Rainfall Intensity (ml)"))

# Get summaries with CIs and p-values for each
trend_summaries <- trends_list %>%
  map(~ summary(.x, infer = c(TRUE, TRUE)))

#check converenge is not a problem - all 0s so we're fine
map(models, ~ .x$fit$convergence)

#check dispersion
# Simulate residuals for every model in the list
sim_outputs <- models %>%
  map(~ simulateResiduals(fittedModel = .x, n = 1000))

# Run the dispersion test on each simulation output
dispersion_tests <- sim_outputs %>%
  map(~ testDispersion(.x, plot = FALSE))

# Pull out the key stats (statistic and p-value) into one tidy dataframe
dispersion_summary <- imap_dfr(dispersion_tests, ~ data.frame(
  Model      = .y,
  Dispersion = .x$statistic,
  p_value    = .x$p.value
))
#check dispersion of models - all looks good
dispersion_summary

#view results - suggest differences under bracken...
trend_summaries

#compare with models...reveal differences not significant
models %>% map(~ Anova(.x, type = 2))
