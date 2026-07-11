# Load libraries
library(ggplot2)
library(tidyr)
library(readxl) #for exporting tabs from xlsx to csv
library(tibble) #for renaming rows
library(readr) #for reading in many files
library(dplyr) #for making factor columns
library(stringr) #for making factor columns
library(glmmTMB) #for modelling results
library(svglite)

library(tidyverse)
library(here)

#library(lubridate)
library(ggpubr)
library(multcompView) #for significant difference letters
library(scales)
library(dplyr)
library(rlang)

#load in microplate, seprate tabs into csvs to analyse ----

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



## get the standard curve ---- ----

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

## plot the regression slopes ---- ----


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



# export the slope ----
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

#### analyse the microplate data ----

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

#reformce so we know which enzyme which data refers to
results_df$SampleType <- gsub("\\.csv$", "", results_df$SourceFile)  # remove ".csv"
results_df$SampleType <- gsub("[0-9]", "", results_df$SampleType)     # remove digits
#remove rows where the microplate wells were blank
results_df <- results_df[nchar(results_df$SampleID) > 2, ]

#load and append covariates needed e.g. moisture, soil mass, dilution factor ----
#covaraites
covs <- read_csv("Data/Enzyme_variables.csv")
#check all SampleIDs in the results_df are consistently formatted with their -M suffixes
results_df$SampleID <- gsub("(?<!-)\\s?M$", "-M", results_df$SampleID, perl = TRUE)
#ensure Sample ID is named correctly
colnames(results_df)[1] <- "Sample ID"
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
#calculate extracellular enzyme activity per hour
results$`EEA per hour` <- ((results$`correction for concentration step (sausages)`/results$`Dry soil mass (g)`)*1000)/results$`Incubation time (hours)`

#add columns for the various factors
results <- results %>%
  mutate(
    #if M is present, it is from end of mesocosm study.  If no M, from field conditions
    Condition = if_else(str_detect(`Sample ID`, "-M$"), "Mesocosm", "Field"),
    
    # extract the 2-letter code (HB, HH, GG, GB) regardless of whether -M is present
    HabitatCode = str_extract(`Sample ID`, "[A-Z]{2}(?=(-M)?$)"),
    
    Habitat = case_when(
      substr(HabitatCode, 1, 1) == "H" ~ "Heathland",
      substr(HabitatCode, 1, 1) == "G" ~ "Grassland",
      TRUE ~ NA_character_
    ),
    
    Bracken = case_when(
      substr(HabitatCode, 2, 2) == "B" ~ "Present",
      substr(HabitatCode, 2, 2) %in% c("G", "H") ~ "Absent",
      TRUE ~ NA_character_
    )
  ) %>%
  select(-HabitatCode)  # drop the helper column, keep only the three requested

results
#### generate boxplot of EEAs under inital conditons ----
#EEA from initial field conditions
field_df <- results %>%
  filter(Condition == "Field")

#plot the intial EEAs
initial <- ggplot(field_df, aes(x = Habitat, y = `EEA per hour`, fill = Bracken)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), 
             aes(color = Bracken), alpha = 0.6, size = 1.5) +
  facet_wrap(~ SampleType) +
  scale_fill_manual(values = c("Absent" = "sienna", "Present" = "limegreen")) +
  scale_color_manual(values = c("Absent" = "sienna", "Present" = "limegreen")) +
  labs(x = "Habitat", y = "EEA per hour") +
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        strip.background = element_rect(fill = "grey85", color = "black"))

show(initial)
#save the figure
ggsave("Figures/initial_EEA_per_hour_by_habitat_bracken.svg", plot = initial, 
       width = 10, height = 7, dpi = 300)

#### EEAs under end of mesocosm conditions ----

#load in the moisture data, format data correctly
field_moisture <- read_csv("Data/Soil Moisture (field conditions).csv")
#field condition moisture
initial_moisture <- field_moisture[, c(1, 7)]

moisture_data <- read_csv("Data/mesocosm moisture - Soil mass change as % of day 0.csv")
                    


#add %change in mass (i.e. moisture)
initial_moisture <- initial_moisture %>%
  left_join(moisture_data %>% select(1, 21), by = "Sample ID")

initial_moisture$`Gravimetric moisture content when sampling for enzymes` <- initial_moisture$`Soil moisture (% wet mass)`*(initial_moisture$`Soil Moisture (% day 0) destructive sampling (13 or 14/11/2025)`/100)

initial_moisture$`Sample ID` <- ifelse(!grepl("-M$", initial_moisture$`Sample ID`),
                             paste0(initial_moisture$`Sample ID`, "-M"),
                             initial_moisture$`Sample ID`)
      

#load covariates
covs <- read_csv("Data/Enzyme_variables.csv")

#filter to final day of mesocosm study







#run models to see if differences are significant ----

#first split out data for each enzyme
sample_types <- unique(field_df$SampleType)
split_dfs <- field_df %>%
  group_by(SampleType) %>%
  group_split() %>%
  setNames(sample_types)

gluc <- split_dfs[["Gluc"]] 

#check the data distripution for each enzyme - they all look right skewed!
hist(split_dfs[["Xylo"]]$`EEA per hour`)
print(split_dfs[["Xylo"]], n = 31)

library(purrr)

models <- split_dfs %>%
  map(~ glmmTMB(`EEA per hour` ~ Habitat*Bracken, 
                family = Gamma(link = "log"), 
                data = .x))

# View summary for each model
models %>% map(summary)

#INITIAL CODE WORKING ON ONE MICROPLATE
#### now analyse each microplate gain 782----
Gluc1 <- read_csv("Data/Microplate Reader Data/26-06-2026 plates gain782/Gluc1.csv")
#just the data we want
df <- Gluc1[12:20, 1:13]
#name the columns and rows
colnames(df) <- as.character(df[1, ])
#delete the row of names
df <- df[-1, ]
#name the id column, set to row name
names(df)[1] <- "row_id"
df <- column_to_rownames(df, var = "row_id")
#convert to df
df <- as.data.frame(df)
#ensure numbers are coded as numbers
df <- as.data.frame(lapply(df, function(x) as.numeric(as.character(x))), check.names = FALSE)

#compute outliers
#top 4 rows are sample, bottom 4 rows are control
group_size <- 4
n_groups <- ceiling(nrow(df) / group_size)
group_id <- rep(1:n_groups, each = group_size, length.out = nrow(df))

# quantile(..., type = 7) is R's default and matches Excel's QUARTILE.INC
quartile_df <- data.frame(
  SampleID  = character(),
  Group     = integer(),
  Q1        = numeric(),
  Q3        = numeric(),
  IQR       = numeric(),
  LowerBound = numeric(),
  UpperBound = numeric(),
  stringsAsFactors = FALSE
)

for (col in colnames(df)) {
  for (g in unique(group_id)) {
    vals <- df[group_id == g, col]
    q1 <- quantile(vals, 0.25, type = 7, na.rm = TRUE)
    q3 <- quantile(vals, 0.75, type = 7, na.rm = TRUE)
    iqr <- q3 - q1
    lower <- q1 - 1.5 * iqr
    upper <- q3 + 1.5 * iqr
    quartile_df <- rbind(quartile_df,
                         data.frame(SampleID = col, Group = g,
                                    Q1 = q1, Q3 = q3, IQR = iqr,
                                    LowerBound = lower, UpperBound = upper))
  }
}

rownames(quartile_df) <- NULL
#rename so we can easily see which data is the Sample (top 4 microplate rows), which is the control (bottom 4 microplate rows)
quartile_df$Group <- factor(quartile_df$Group, levels = c(1, 2), labels = c("Extract", "Control"))
quartile_df

# Remove outliers, using the same Sample/Control grouping
group_labels <- c(`1` = "Extract", `2` = "Control")

df_no_outliers <- df  # start fresh from the original data
#log to store the outliers removed
outlier_log <- data.frame(
  RowLabel   = character(),
  SampleID   = character(),
  Group      = character(),
  Value      = numeric(),
  LowerBound = numeric(),
  UpperBound = numeric(),
  stringsAsFactors = FALSE
)
#remove any outliers
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
        RowLabel   = rownames(df)[rows_in_group][is_outlier],
        SampleID   = col,
        Group      = g_label,
        Value      = vals[is_outlier],
        LowerBound = lower,
        UpperBound = upper
      ))
    }
    
    vals[is_outlier] <- NA
    df_no_outliers[[col]][rows_in_group] <- vals
  }
}

rownames(outlier_log) <- NULL
df_no_outliers
outlier_log

#average the sample and control for each

#load in the standard curve slope

#### now analyse each microplate gain 782----
stnd_slopes <- read_csv("Processed Data/mean_slope_output.csv")

#create empty dataframe
result_df <- data.frame(
  SampleID    = character(),
  SampleMean  = numeric(),
  ControlMean = numeric(),
  Difference  = numeric(),
  Concentration = numeric(),
  stringsAsFactors = FALSE
)

for (col in colnames(df_no_outliers)) {
  ##average the extract and cnotrol abordbances
  sample_mean  <- mean(df_no_outliers[[col]][group_id == 1], na.rm = TRUE)
  control_mean <- mean(df_no_outliers[[col]][group_id == 2], na.rm = TRUE)
  #subtract control from extract
  difference   <- sample_mean - control_mean
  #divide by standarc curve slope to get moles per well
  concentration <- difference / stnd_slopes[1, 2]
  #bind together into table
  result_df <- rbind(result_df, data.frame(
    SampleID    = col,
    SampleMean  = sample_mean,
    ControlMean = control_mean,
    Difference  = difference,
    Concentration = concentration
  ))
}

rownames(result_df) <- NULL
result_df
