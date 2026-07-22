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
df_RDA <- read_csv("Processed Data/RDa data.csv")


d <- as.data.frame(readr::read_csv(
  here::here("data", "2026-01-29_all-variables_masterfile.csv")
))
d <- d[-1]
#make sure our variables are coded as factors
d$Vegetation <- factor(d$Vegetation, levels = c("Bracken", "Heather"), labels = c("Bracken", "Heather"))

#convert the sample date to an appropriate format
d$`Sample Date` <- as.Date(d$`Date Sampled`, format = "%d/%m/%Y")
#convert to Julian date for linear temporal trends
d$julian <- scale(as.numeric(d$`Sample Date`))


#replace row index with sample names
rownames(d) <- d[,1]
#the community data frame for standardized abundances - with thsi we can explain approx 35% of the variance
cdf <- d[,(37:44)]
#community dataframe for raw abundances - with this we can explain only 15% of the variance
#cdf <- d[, (37:44)]

#replace null (empty excell cell) with "0"
cdf[is.na(cdf)] <- 0
cdf <- as.matrix(cdf)

#explanatory data frame: paramters that differed significantly.  THese have different base units so we may want to standardize them e.g. by z scoring
#elevation, %bracken, %heather, moisture, pH, latitude, longitude, TNb, total C, total N, C:N ratio, alpha, SUVA.  Add WEOC, WEN to dataframe too, and date
edf <- d[, c(5, 6, 7, 22, 28, 29, 30, 31, 32, 33, 34, 35, 36, 78, 79, 80, 81, 82, 84)]
colnames(edf) <- c("elevation", "% bracken", "% heather", "longitude", "soil moisture", "pH", "WEOC", "WEN", "total carbon", "total nitrogen", "C:N", "alpha", "SUVA", "groundfrost days", "relative humidity", "total rainfall", "partial pressure of water vapour", "snow lying days", "sampling date")
#PCA of edf to see if any explanatory factors are highly correlated
correlation_check <- cor(edf)
cor(edf)
# Step 2: Remove columns with all NA
env_data_clean <- edf[, colSums(!is.na(edf)) > 0]

# Step 3: Remove columns with zero or near-zero variance
# Use a small threshold to catch nearly constant columns
# env_data_clean <- env_data_clean[, apply(env_data_clean, 2, function(x) {
#   x <- x[is.finite(x)]  # remove Infs
#   if (length(x) == 0) return(FALSE)
#   sd_val <- sd(x, na.rm = TRUE)
#   return(!is.na(sd_val) && sd_val > 1e-8)
# })]

# Step 4: Drop rows with any NA (PCA can't handle them)
env_data_clean <- na.omit(env_data_clean)

# Step 5: Scale and run PCA
scaled_env_data <- scale(env_data_clean)
pca_result <- prcomp(scaled_env_data, center = TRUE, scale. = TRUE)

# Step 6: View results
summary(pca_result)

#sample_ids <- d$`Sample ID`[as.numeric(rownames(env_data_clean))]

pca_scores <- as.data.frame(pca_result$x)
pca_scores$SampleID <- d$`Sample ID`
pca_scores$Vegetation <- d$Vegetation

#plot the loadings
# Extract loadings for PC1 and PC2
loadings <- as.data.frame(pca_result$rotation[, 1:2])
loadings$Variable <- rownames(loadings)

# Plot loadings as arrows from origin
ggplot(loadings, aes(x = PC1, y = PC2, label = Variable)) +
  geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")), color = "blue") +
  geom_text(hjust = 0.5, vjust = -0.7, size = 4) +
  xlim(c(min(loadings$PC1, 0) * 1.1, max(loadings$PC1, 0) * 1.1)) +
  ylim(c(min(loadings$PC2, 0) * 1.1, max(loadings$PC2, 0) * 1.1)) +
  labs(title = "PCA Loadings Plot (PC1 & PC2)", x = "PC1", y = "PC2") +
  theme_minimal()


#visualize the biplot
# Extract PCA scores and add SampleID if you have it
#pca_scores <- as.data.frame(pca_result$x[, 1:2])
#if (exists("sample_ids")) pca_scores$SampleID <- sample_ids

# Scale loadings for plotting (to match score scale roughly)
loadings_scaled <- loadings
scale_factor <- min(
  (max(pca_scores$PC2) - min(pca_scores$PC2)) / (max(loadings$PC2) - min(loadings$PC2)),
  (max(pca_scores$PC1) - min(pca_scores$PC1)) / (max(loadings$PC1) - min(loadings$PC1))
)
loadings_scaled$PC1 <- loadings$PC1 * scale_factor * 0.7
loadings_scaled$PC2 <- loadings$PC2 * scale_factor * 0.7

pca_scores$Vegetation <- factor(pca_scores$Vegetation)
#add variance explained by each PC
# Calculate proportion of variance explained
pca_var <- pca_result$sdev^2
pca_var_explained <- round(100 * pca_var / sum(pca_var), 1)  # in %
#number sites so we can easily tell which are from which site
id_map <- c("Bridestones" = 6, "Scarth Wood Moor" = 5, "Brimham Moor" = 4, "Widdybanks" = 3, "Haweswater" = 2, "Whiteside" = 1)
pca_scores$Site <- d$Site
pca_scores$Site <- id_map[pca_scores$Site]


library(ggrepel)
# Plot
ggplot() +
  geom_point(data = pca_scores, aes(x = PC1, y = PC2, color = Vegetation), size = 3) +
  geom_text(data = pca_scores, aes(x = PC1, y = PC2, label = Site), alpha = 0.3, vjust = -1) +
  geom_segment(data = loadings_scaled, aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.3, "cm")), color = "red") +
  geom_text_repel(data = loadings_scaled, aes(x = PC1, y = PC2, label = Variable),
                  color = "black", hjust = 0.5, vjust = -0.7) +
  labs(x = paste0("PC1 (", pca_var_explained[1], "%)"),
       y = paste0("PC2 (", pca_var_explained[2], "%)"),
       color = "Vegetation"
  ) + scale_color_manual(values = c("Bracken" = "limegreen", "Heather" = "#AA4499")) + 
  theme_minimal() 

#do linear regressions of longitude again hydrological variables to see the strength that we can see indicated by the PCA
# Regression model
model <- lm(`partial pressure of water vapour` ~ longitude, data = scaled_env_data)
summary(model)


#assign the treatments to relevant rows of the dataframe
#edf$treatment <- c(rep("Brimham Bracken",10),rep("Brimham Heather",10), rep("Bridestones Bracken",10),rep("Bridestones Heather",10), rep("Haweswater Bracken", 10), rep("Haweswater Heather", 10), rep("Scarth Wood Bracken",10),rep("Scarth Wood Heather",10), rep("Widdybanks Bracken",10),rep("Widdybanks Heather",10), rep("Whiteside Bracken", 10), rep("Whiteside Heather", 10))
#edf$treatment <- c(rep("Bracken",10),rep("Heather",10), rep("Bracken",10),rep("Heather",10), rep("Bracken", 10), rep("Heather", 10), rep("Bracken",10),rep("Heather",10), rep("Bracken",10),rep("Heather",10), rep("Bracken", 10), rep("Heather", 10))
#all explanatory factors modelled

#check factors for correlation e.g. are elevation and LongitudeE correlated? seems likely!
scaled_env_data <- as.data.frame(scaled_env_data)
#dbRDA using normaliesd environmental factors, removed julian
dbrda_summary <- dbrda(formula = cdf ~ `pH` + `longitude` + `WEOC` + `WEN` + `total carbon` + `C:N` + `alpha` + `SUVA` + `relative humidity` + `total rainfall` + `snow lying days` +`sampling date` + `partial pressure of water vapour` + `soil moisture`, scaled_env_data, distance = "euclidean", sqrt.dist = FALSE, add = FALSE, dfun = vegdist, metaMDSdist = FALSE, na.action = na.exclude, subset = NULL)


#just the environmental factors that are significant
#dbrda_summary <- dbrda(formula = cdf ~ `Elevation (m)` + `LongitudeE` + `CN ratio` + `Groundfrost days` + `Relative humidity (%)`, edf, distance = "euclidean", sqrt.dist = FALSE, add = FALSE, dfun = vegdist, metaMDSdist = FALSE, na.action = na.exclude, subset = NULL)

summary(dbrda_summary)


# Define treatment variable and convert to factor
#treatment <- as.factor(edf$treatment)

treatment <- factor(d$Vegetation, levels = unique(d$Vegetation))
# Define custom colors and point shapes (pch) for the 6 treatments
treatment_levels <- levels(treatment)
#colours for each treatment
#colors =c("#999999","#999999", "#E69F00", "#E69F00", "#56B4E9","#56B4E9", "#009E73","#009E73", "#CC79A7", "#CC79A7", "#0072B2","#0072B2")[seq_along(treatment_levels)]
colors = c("limegreen", "#AA4499")
#shapes for point codes
#pchs<- c(15, 0, 16, 1, 17,2, 18, 3, 19, 4, 20, 5)[seq_along(treatment_levels)]
pchs<- c(15, 16)[seq_along(treatment_levels)]
# Define cex values (make squares smaller, circles normal)
cexs <- c(1.2, 1.2)[seq_along(treatment_levels)]  # adjust sizes as needed

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
legend(x = 12, y = 23, legend = treatment_levels, 
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


#ellipses for the sites only


# Define treatment variable and convert to factor
#treatment <- as.factor(edf$treatment)

#rename the sites
d <- d %>%
  mutate(Site = case_when(
    Site == "Brimham Moor" ~ "Brimham Moor SSSI",
    Site == "Bridestones" ~ "Bridestone Moor SSSI",
    Site == "Haweswater" ~ "Naddle Forest SSSI",
    Site == "Scarth Wood Moor" ~ "Scarth Wood Moor SSSI",
    Site == "Widdybanks"        ~ "Moor House - Upper Teesdale NNR",
    Site == "Whiteside"         ~ "Whiteside SSSI",
    TRUE ~ Site  # keeps anything else unchanged
  ))


treatment <- factor(d$Site, levels = unique(d$Site))
treatment <- factor(treatment, levels = c("Bridestone Moor SSSI", "Scarth Wood Moor SSSI", "Brimham Moor SSSI", "Moor House - Upper Teesdale NNR", "Naddle Forest SSSI", "Whiteside SSSI"))
# Define custom colors and point shapes (pch) for the 6 treatments
treatment_levels <- levels(treatment)
#colours for each treatment
colors =c("red","orange", "yellow", "green", "blue","violet")[seq_along(treatment_levels)]
#colors = c("limegreen", "#AA4499")
#shapes for point codes
#pchs<- c(15, 0, 16, 1, 17,2, 18, 3, 19, 4, 20, 5)[seq_along(treatment_levels)]
pchs<- c(17, 17, 17, 17, 17, 17)[seq_along(treatment_levels)]
# Define cex values (make squares smaller, circles normal)
#cexs <- c(1.2, 1.2, 1.2, 1.2,1.2,1.2)[seq_along(treatment_levels)]  # adjust sizes as needed

# Get variance explained by each axis
eig_vals <- eigenvals(dbrda_summary)
var_explained <- eig_vals / sum(eig_vals) * 100
axis_labels <- paste0("dbRDA", 1:2, " (", round(var_explained[1:2], 1), "%)")



dev.new()
#pdf("figures/_2026-02-17_dbRDA_plot_sites.pdf", width = 9, height = 7)  
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
ordiellipse(dbrda_summary, groups = treatment, display = "sites", kind = "se", conf = 0.95, draw = "polygon", col = colors, border = colors,lwd = 1.5,lty = 1, alpha = 80)  # transparency (0–255); needs vegan >= 2.6-4


# Add legend
legend(x = -37, y = -15, legend = treatment_levels, 
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




#ellipses for the sites x bracken - a complete mess
#create new column of Vegeation:Site
# Create a new column combining Vegetation and Site
d$Veg.Site <- paste(d$Site, d$Vegetation, sep = "-")


# Define treatment variable and convert to factor
#treatment <- as.factor(edf$treatment)

treatment <- factor(d$Veg.Site, levels = unique(d$Veg.Site))
treatment <- factor(treatment, levels = c("Bridestones-Bracken","Bridestones-Heather", "Scarth Wood Moor-Bracken", "Scarth Wood Moor-Heather", "Brimham Moor-Bracken", "Brimham Moor-Heather", "Widdybanks-Bracken", "Widdybanks-Heather", "Haweswater-Bracken", "Haweswater-Heather", "Whiteside-Bracken", "Whiteside-Heather"))
# Define custom colors and point shapes (pch) for the 6 treatments
treatment_levels <- levels(treatment)
#colours for each treatment
colors =c("red", "red", "orange", "orange", "yellow", "yellow", "green", "green", "blue", "blue", "violet", "violet")[seq_along(treatment_levels)]
#colors = c("limegreen", "#AA4499")
#shapes for point codes
#pchs<- c(15, 0, 16, 1, 17,2, 18, 3, 19, 4, 20, 5)[seq_along(treatment_levels)]
pchs<- c(15, 16, 15, 16, 15, 16, 15, 16, 15, 16, 15, 16)[seq_along(treatment_levels)]
# Define cex values (make squares smaller, circles normal)
cexs <- c(1.2, 1.2, 1.2, 1.2,1.2,1.2, 1.2, 1.2, 1.2, 1.2,1.2,1.2)[seq_along(treatment_levels)]  # adjust sizes as needed

# Get variance explained by each axis
eig_vals <- eigenvals(dbrda_summary)
var_explained <- eig_vals / sum(eig_vals) * 100
axis_labels <- paste0("dbRDA", 1:2, " (", round(var_explained[1:2], 1), "%)")



dev.new()
pdf("figures/_2026-02-17_dbRDA_plot_sitesxbracken.pdf", width = 9, height = 7)  
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
legend(x = 10, y = 18, legend = treatment_levels, 
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
