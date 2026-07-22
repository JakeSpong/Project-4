#### create combined dataframe for RDA ----
#load in the EEA data
results <- read_csv("Processed Data/EEA (nmol per h per g dry weight).csv")
#load in explanatory covaraite data
expcovs <- read_csv("Data/Explanatory covariates.csv")
#add explanatory covariate data to the results df, matching regardless of -M suffix
results <- results %>%
  mutate(join_id = str_remove(`Sample ID`, "-M$")) %>%
  left_join(expcovs, by = c("join_id" = "Sample ID")) %>%
  dplyr::select(-join_id)
# enzymes from inital field conditions
df <- results %>% filter(Condition == "Field")
#just the data we acutally need
df <- df[c(1,9,20,22,23, 24, 31, 32, 33,34, 35)]

#load in gas flux data
flux_data <- read_csv("Processed Data/all_times_gas_fluxes.csv")
#remove data entries with p values < 0.05
#flux_data <- flux_data[flux_data$CO2_p_value <= 0.05,]
flux_timeseries <- flux_data[c(1,2, 12, 13, 32, 33, 34)]

#get the date only in one column
flux_timeseries <- flux_timeseries %>%
  mutate(date = as.Date(start_time))

#just inital gas fluxes
flux <- flux_timeseries %>%
  filter(date == as.Date("2025-09-29"))
#just the sample IDs and flux values
df_flux <- flux[c(1, 6, 7)]
#rename the sample ID column so we can merge dfs
colnames(df_flux)[colnames(df_flux) == "sample_id"] <- "Sample ID"
#left join
df <- left_join(df, df_flux, by = c("Sample ID"))
#get in the field soil moisture data!
moist <- read_csv("Data/Soil Moisture (field conditions).csv")
#just the data we want
moist_df <- moist[c(1, 7)]
#left join
df <- left_join(df, moist_df, by = c("Sample ID"))
#convert from long to wide for EEA data
df_wide <- df %>%
  pivot_wider(
    names_from = SampleType,
    values_from = `EEA per hour per g dry soil`
  )
#save combined df
write.csv(df_wide, "Processed Data/RDA data.csv", row.names = FALSE)

#### load in RDA data ----
#load in the EEA data
d <- read_csv("Processed Data/RDa data.csv")
#make sure our variables are coded as factors
d$Bracken <- factor(d$Bracken, levels = c("Present", "Absent"), labels = c("Present", "Absent"))
d$Habitat <- factor(d$Habitat, levels = c("Heathland", "Grassland"), labels = c("Heathland", "Grassland"))
df_PCA <- d[, -c(1, 2, 3, 8, 9)]
#ensure all columns are numeric
df_PCA[] <- lapply(df_PCA, as.numeric)

library(stats) #for prcomp, for the PCA
PCA <- prcomp(df_PCA, retx = TRUE, center = TRUE, scale. = TRUE, tol = NULL, rank. = NULL)
library(ggbiplot)

ggbiplot(PCA) +
  geom_point(aes(colour = d$Bracken, shape = d$Habitat, size = 2.0)) +
  labs(colour = "Bracken") +
  theme_bw()




#dbRDA working but pointless...we wnt PCA

#d <- as.matrix(d)
#replace row index with sample names
#rownames(d) <- d[,1]
#d <- as.tibble(d)
#explanatory data frame: paramters that differed significantly.  THese have different base units so we may want to standardize them e.g. by z scoring
#gas fluxes and enzyume data
cdf <- d[, c(10, 11, 13, 14, 15, 16)]
#the explanaotary variables
edf <- d[, c(4, 5, 6, 7, 12)]
colnames(edf) <- c("elevation", "latitude", "longitude", "pH", "soil moisture")
#ensure all columns are numeric
edf <- edf %>%
  mutate(across(c("elevation", "latitude", "longitude", "pH", "soil moisture"), as.numeric))
#PCA of edf to see if any explanatory factors are highly correlated
correlation_check <- cor(edf)
cor(edf)
# Step 2: Remove columns with all NA
env_data_clean <- edf[, colSums(!is.na(edf)) > 0]
# Step 4: Drop rows with any NA (PCA can't handle them)
env_data_clean <- na.omit(env_data_clean)

# Step 5: Scale 
scaled_env_data <- scale(env_data_clean)

#check factors for correlation e.g. are elevation and LongitudeE correlated? seems likely!
scaled_env_data <- as.data.frame(scaled_env_data)
library(vegan) #for dbRDA

#dbRDA using normaliesd environmental factors, removed julian
dbrda_summary <- dbrda(formula = cdf ~ `elevation` + `latitude` + `longitude` + `pH` + `soil moisture`, scaled_env_data, distance = "euclidean", sqrt.dist = FALSE, add = FALSE, dfun = vegdist, metaMDSdist = FALSE, na.action = na.exclude, subset = NULL)


#just the environmental factors that are significant
#dbrda_summary <- dbrda(formula = cdf ~ `Elevation (m)` + `LongitudeE` + `CN ratio` + `Groundfrost days` + `Relative humidity (%)`, edf, distance = "euclidean", sqrt.dist = FALSE, add = FALSE, dfun = vegdist, metaMDSdist = FALSE, na.action = na.exclude, subset = NULL)

summary(dbrda_summary)


# Define treatment variable and convert to factor
#treatment <- as.factor(edf$treatment)
d$Treatment <- paste(d$Habitat, d$Bracken, sep = "-")
treatment <- factor(d$Treatment, levels = unique(d$Treatment))
# Define custom colors and point shapes (pch) for the 6 treatments
treatment_levels <- levels(treatment)
#colours for each treatment
#colors =c("#999999","#999999", "#E69F00", "#E69F00", "#56B4E9","#56B4E9", "#009E73","#009E73", "#CC79A7", "#CC79A7", "#0072B2","#0072B2")[seq_along(treatment_levels)]
colors = c("#AA4499", "#AA4499", "limegreen", "limegreen")
#shapes for point codes
#pchs<- c(15, 0, 16, 1, 17,2, 18, 3, 19, 4, 20, 5)[seq_along(treatment_levels)]
pchs<- c(16, 15, 16, 15)[seq_along(treatment_levels)]
# Define cex values (make squares smaller, circles normal)
cexs <- c(1.2, 1.2, 1.2, 1.2)[seq_along(treatment_levels)]  # adjust sizes as needed

# Get variance explained by each axis
eig_vals <- eigenvals(dbrda_summary)
var_explained <- eig_vals / sum(eig_vals) * 100
axis_labels <- paste0("dbRDA", 1:2, " (", round(var_explained[1:2], 1), "%)")



dev.new()
#pdf("figures/2026-02-17_dbRDA_plot.pdf", width = 9, height = 7)  
# Base plot: empty dbRDA ordination
plot(dbrda_summary, type = "n", scaling = 2, , 
     xlab = axis_labels[1], ylab = axis_labels[2])  # Use scaling = 2 for species-environment biplot

# Add site points, colored by treatment
for (i in seq_along(treatment_levels)) {
  sel <- treatment == treatment_levels[i]
  points(scores(dbrda_summary, display = "sites", scaling = 2)[sel, ], 
         col = colors[i], pch = pchs[i], cex = cexs[i])
} 
# --- Add ellipses around treatment groups ---
ordiellipse(dbrda_summary, groups = treatment, display = "sites", kind = "se", conf = 0.95, draw = "polygon", col = colors, border = colors,lwd = 1.5,lty = 1, alpha = 60)  # transparency (0–255); needs vegan >= 2.6-4


# Add legend
legend(x = 5, y = 10, legend = treatment_levels, 
       col = colors, pch = pchs)

# Extract biplot scores of environmental variables
env_vectors <- scores(dbrda_summary, display = "bp", scaling = 2)
# Scale factor to adjust vector length visually
vec_multiplier <- ordiArrowMul(env_vectors)  # automatic scaling
# Add arrows
apply(env_vectors, 1, function(row) {
  arrows(0, 0, row[1] * vec_multiplier, row[2] * vec_multiplier, 
         length = 0.1, col = "black")
})
# Add vector labels
text(env_vectors * vec_multiplier, labels = rownames(env_vectors), col = "black", pos = 4, cex = 0.8)


dev.off()

summary(dbrda_summary)
#is the model significant?
anova(dbrda_summary)
#test axes for significance
anova(dbrda_summary, by = "axis", perm.max = 500)
#test environmental variables for significance
anova(dbrda_summary, by = "terms", perm.max = 500)


