library(maptiles)
library(terra)
library(sf)
library(ggplot2)

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

#for models
library(glmmTMB)
library(splines) #to allow for curved models
library(ggeffects)
library(emmeans)

library(marginaleffects) #for converting model output to change in flux per day
library(broom.mixed) #for emtrends
library(DHARMa)

#### generate maps of sample locations (Figure 1) ----
#to get google satellite layer, Go to the Google Cloud Console.  Create a new project (or use an existing one). Enable the Maps Static API and Geocoding API.  Get your API key.


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


#### pH ----
pH <- read_csv("Data/pH.csv")


#reorder the sites so they show up on the plot from west (LHS) to east (RHS)
pH$Bracken <- factor(pH$Bracken, levels = c("Present", "Absent"))
#reorder the habitats so they appear in the specified order
pH$Habitat <- factor(pH$Habitat, levels = c("Rainfall", "Grassland", "Heathland"))

#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
pH_lys_bxp <- ggboxplot(pH, x = "Habitat", aes(y = `Lysate pH`), color = "Bracken", lwd = 0.75)  +
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


pH <- pH[-(1:3),]
#boxplot the data. Use aes() with backticks (``) so avoid an error with our column name
pH_bxp <- ggboxplot(pH, x = "Habitat", aes(y = `Field soil pH`), color = "Bracken", lwd = 0.75)  +
  labs(
    x = "Habitat",
    y =  expression("Field soil pH")
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

show(pH_bxp)  

hist(pH$`Field soil pH`)


#trim off the rainfall data
model_data <- pH_lys[-(1:3),]
hist(model_data$`Lysate pH`)

model <- glmmTMB(
  `Lysate pH` ~ Bracken*Habitat,
  data = model_data,
  family = gaussian()
)
summary(model)

#check for model overdispersion
sim_res <- simulateResiduals(model)
plot(sim_res)   # KS test + dispersion plot
testDispersion(sim_res)
testUniformity(sim_res) # checks overall residual uniformity
testOutliers(sim_res)   # checks for outliers


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


#### model CO2 flux before and after rainfall FIGURE 2 TOP----
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

#model only prior to the rainfall treatment

# Create a combined grouping label
flux_avg <- flux_avg %>%
  mutate(Habitat_Bracken = paste(Habitat, Bracken, sep = " - "))
#rename the variables to avoid modelling issue with glmmTMB
colnames(flux_avg)[colnames(flux_avg) == "Average CO2 flux"] <- "average_CO2_flux"
#rename the variables to avoid modelling issue with glmmTMB
colnames(flux_avg)[colnames(flux_avg) == "Average CH4 flux"] <- "average_CH4_flux"
#convert dates to numbers
flux_avg$date <- as.Date(flux_avg$date)
flux_avg$date_num <- as.numeric(flux_avg$date)
#ensure factors are coded properly
flux_avg$Habitat <- factor(flux_avg$Habitat)
flux_avg$Bracken <- factor(flux_avg$Bracken)

#assign a new variable, Pre/Post rainfall, to compare fluxes pre rainfall to fluxes post
flux_before <- flux_avg %>%
  filter(date < as.Date("2025-10-14"))



hist(flux_before$average_CO2_flux)

#splines allows for curved model.  df = 1 is linear, increasing number increases curviness - check AIC comparisons below, and check visualisations, to determine best model
co2_model <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 1) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)
# check residuals 
simulateResiduals(co2_model, plot = TRUE)
#check model outputs
summary(co2_model)
#get actual change in flux per day

#convert model slope to change in flux per day
co2_slopes_before <- slopes(
  co2_model,
  variables = "date_num",
  by = c("Habitat", "Bracken")
)


#create prediction data
newdatbefore <- expand.grid(
  date_num = seq(
    min(flux_before$date_num),
    max(flux_before$date_num),
    length.out = 200
  ),
  Habitat = levels(flux_before$Habitat),
  Bracken = levels(flux_before$Bracken)
)

newdatbefore$Habitat <- factor(
  newdatbefore$Habitat,
  levels = levels(flux_before$Habitat)
)

newdatbefore$Bracken <- factor(
  newdatbefore$Bracken,
  levels = levels(flux_before$Bracken)
)

#use the model to predict
predbefore <- predict(
  co2_model,
  newdata = newdatbefore,
  type = "response",
  se.fit = TRUE
)

newdatbefore$fit <- predbefore$fit
newdatbefore$lower <- predbefore$fit - 1.96 * predbefore$se.fit
newdatbefore$upper <- predbefore$fit + 1.96 * predbefore$se.fit

newdatbefore$date <- as.Date(
  newdatbefore$date_num,
  origin = "1970-01-01"
)

#create grouping variable
newdatbefore$Group <- interaction(
  newdatbefore$Habitat,
  newdatbefore$Bracken,
  sep = "-"
)

flux_before$Group <- interaction(
  flux_before$Habitat,
  flux_before$Bracken,
  sep = "-"
)

#set shape and line formatting for the different treatemnts
shape_vals <- c("Absent" = 17,  # triangle
                "Present" = 16) # circle

line_vals <- c("Absent" = "dashed",
               "Present" = "solid")



co2_before <- ggplot() +
  
  # observed data
  geom_point(
    data = flux_before,
    aes(
      x = date,
      y = average_CO2_flux,
      colour = Habitat,
      shape = Bracken
    ),
    alpha = 0.8,
    size = 2
  ) +
  
  # model ribbon
  geom_ribbon(
    data = newdatbefore,
    aes(
      x = date,
      ymin = lower,
      ymax = upper,
      fill = Habitat,
      group = interaction(Habitat, Bracken)
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  # model lines
  geom_line(
    data = newdatbefore,
    aes(
      x = date,
      y = fit,
      colour = Habitat,
      linetype = Bracken,
      group = interaction(Habitat, Bracken)
    ),
    linewidth = 1
  ) + 
  coord_cartesian(ylim = c(-0.1, 3)) +
  
  # scales
  scale_colour_manual(values = c("Grassland" = "#1b9e77",
                                 "Heathland" = "#AA4499")) +
  
  scale_fill_manual(values = c("Grassland" = "#1b9e77",
                               "Heathland" = "#AA4499")) +
  
  scale_shape_manual(values = shape_vals) +
  
  scale_linetype_manual(values = line_vals) +
  
  scale_x_date(date_breaks = "3 days", date_labels = "%d %b") +
  
  theme_bw() +
  
  labs(
    x = "Date",
    y = expression("Mean CO"[2] * " flux (" * mu * "mol" ~ s^{-1} ~ m^{-2} * ")"),
    colour = "Habitat",
    fill = "Habitat",
    shape = "Bracken",
    linetype = "Bracken"
  )
#show the figure
show(co2_before)
#save the figure
#ggsave("Figures/co2_before.svg", plot = co2_before, width = 7.5, height = 6, dpi = 300)


#compare AICs
m1 <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 1) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)
m2 <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 2) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)

m3 <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 3) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)

m4 <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 4) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)

m5 <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 5) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)
#linear model bascially the same as the others, apart from overfitted m5, so just use linear models for this bit
AIC(m1, m2, m3, m4, m5)

# model CO2 flux after rainfall 
#filter for fluxes after rainfall added
flux_after <- flux_avg %>%
  filter(date > as.Date("2025-10-14"))


hist(flux_after$average_CO2_flux)

#library(splines)
#splines allows for curved model
co2_model <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 1) * Habitat * Bracken,
  data = flux_after,
  family = gaussian()
)
# check residuals 
simulateResiduals(co2_model, plot = TRUE)

#check model outputs
summary(co2_model)
#get actual change in flux per day

#convert model slope to change in flux per day
co2_slopes_after <- slopes(
  co2_model,
  variables = "date_num",
  by = c("Habitat", "Bracken")
)
#show the change in slopes per day
summary(slopes)



#create prediction data
newdatafter <- expand.grid(
  date_num = seq(
    min(flux_after$date_num),
    max(flux_after$date_num),
    length.out = 200
  ),
  Habitat = levels(flux_after$Habitat),
  Bracken = levels(flux_after$Bracken)
)

newdatafter$Habitat <- factor(
  newdatafter$Habitat,
  levels = levels(flux_after$Habitat)
)

newdatafter$Bracken <- factor(
  newdatafter$Bracken,
  levels = levels(flux_after$Bracken)
)

#use the model to predict
predafter <- predict(
  co2_model,
  newdata = newdatafter,
  type = "response",
  se.fit = TRUE
)

newdatafter$fit <- predafter$fit
newdatafter$lower <- predafter$fit - 1.96 * predafter$se.fit
newdatafter$upper <- predafter$fit + 1.96 * predafter$se.fit

newdatafter$date <- as.Date(
  newdatafter$date_num,
  origin = "1970-01-01"
)

#create grouping variable
newdatafter$Group <- interaction(
  newdatafter$Habitat,
  newdat$Bracken,
  sep = "-"
)

flux_after$Group <- interaction(
  flux_after$Habitat,
  flux_after$Bracken,
  sep = "-"
)

# #set shape and line formatting for the different treatemnts
# shape_vals <- c("Absent" = 17,  # triangle
#                 "Present" = 16) # circle
# 
# line_vals <- c("Absent" = "dashed",
#                "Present" = "solid")
# 


co2_after <- ggplot() +
  
  # observed data
  geom_point(
    data = flux_after,
    aes(
      x = date,
      y = average_CO2_flux,
      colour = Habitat,
      shape = Bracken
    ),
    alpha = 0.8,
    size = 2
  ) +
  
  # model ribbon
  geom_ribbon(
    data = newdatafter,
    aes(
      x = date,
      ymin = lower,
      ymax = upper,
      fill = Habitat,
      group = interaction(Habitat, Bracken)
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  # model lines
  geom_line(
    data = newdatafter,
    aes(
      x = date,
      y = fit,
      colour = Habitat,
      linetype = Bracken,
      group = interaction(Habitat, Bracken)
    ),
    linewidth = 1
  ) +
  coord_cartesian(ylim = c(-0.1, 3)) +
  
  # scales
  scale_colour_manual(values = c("Grassland" = "#1b9e77",
                                 "Heathland" = "#AA4499")) +
  
  scale_fill_manual(values = c("Grassland" = "#1b9e77",
                               "Heathland" = "#AA4499")) +
  
  scale_shape_manual(values = shape_vals) +
  
  scale_linetype_manual(values = line_vals) +
  
  scale_x_date(date_breaks = "3 days", date_labels = "%d %b") +
  
  theme_bw() +
  
  labs(
    x = "Date",
    y = NULL,
    colour = "Habitat",
    fill = "Habitat",
    shape = "Bracken",
    linetype = "Bracken"
  )
#show the figure
show(co2_after)
#panel figure

#combine plots into a single figure
co2_models <- ggarrange(co2_before, co2_after,
                        labels = c("A", "B"),
                        ncol = 2, nrow = 1,
                        #the width of each panel of the multifigure plot
                        widths = c(1,1.6),
                        legend = FALSE)

show(co2_models)
#save the figure
#ggsave("Figures/co2_models.svg", plot = co2_models, width = 10, height = 6, dpi = 300)


#compare AICs
m1 <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 1) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)
m2 <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 2) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)

m3 <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 3) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)

m4 <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 4) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)

m5 <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 5) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)
#linear model bascially the same as the others, apart from overfitted m5, so just use linear models for this bit
AIC(m1, m2, m3, m4, m5)

#interpret the model

summary(co2_model)



#### model CH4 flux before and after rainfall - needs prior tab - FIGURE 2 bottom ----


hist(flux_before$average_CH4_flux)

#library(splines)
#splines allows for curved model
ch4_model <- glmmTMB(
  average_CH4_flux ~ ns(date_num, df = 1) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)

# check residuals 
simulateResiduals(ch4_model, plot = TRUE)
#check model outputs
summary(ch4_model)
#get actual change in flux per day

#convert model slope to change in flux per day
ch4_slopes_before <- slopes(
  ch4_model,
  variables = "date_num",
  by = c("Habitat", "Bracken")
)



#create prediction data
newdatbefore <- expand.grid(
  date_num = seq(
    min(flux_before$date_num),
    max(flux_before$date_num),
    length.out = 200
  ),
  Habitat = levels(flux_before$Habitat),
  Bracken = levels(flux_before$Bracken)
)

newdatbefore$Habitat <- factor(
  newdatbefore$Habitat,
  levels = levels(flux_before$Habitat)
)

newdatbefore$Bracken <- factor(
  newdatbefore$Bracken,
  levels = levels(flux_before$Bracken)
)

#use the model to predict
predbefore <- predict(
  ch4_model,
  newdata = newdatbefore,
  type = "response",
  se.fit = TRUE
)

newdatbefore$fit <- predbefore$fit
newdatbefore$lower <- predbefore$fit - 1.96 * predbefore$se.fit
newdatbefore$upper <- predbefore$fit + 1.96 * predbefore$se.fit

newdatbefore$date <- as.Date(
  newdatbefore$date_num,
  origin = "1970-01-01"
)

#create grouping variable
newdatbefore$Group <- interaction(
  newdatbefore$Habitat,
  newdatbefore$Bracken,
  sep = "-"
)

flux_before$Group <- interaction(
  flux_before$Habitat,
  flux_before$Bracken,
  sep = "-"
)

#set shape and line formatting for the different treatemnts
shape_vals <- c("Absent" = 17,  # triangle
                "Present" = 16) # circle

line_vals <- c("Absent" = "dashed",
               "Present" = "solid")



ch4_before <- ggplot() +
  
  # observed data
  geom_point(
    data = flux_before,
    aes(
      x = date,
      y = average_CH4_flux,
      colour = Habitat,
      shape = Bracken
    ),
    alpha = 0.8,
    size = 2
  ) +
  
  # model ribbon
  geom_ribbon(
    data = newdatbefore,
    aes(
      x = date,
      ymin = lower,
      ymax = upper,
      fill = Habitat,
      group = interaction(Habitat, Bracken)
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  # model lines
  geom_line(
    data = newdatbefore,
    aes(
      x = date,
      y = fit,
      colour = Habitat,
      linetype = Bracken,
      group = interaction(Habitat, Bracken)
    ),
    linewidth = 1
  ) + 
  coord_cartesian(ylim = c(-2.0, 0.5)) +
  
  # scales
  scale_colour_manual(values = c("Grassland" = "#1b9e77",
                                 "Heathland" = "#AA4499")) +
  
  scale_fill_manual(values = c("Grassland" = "#1b9e77",
                               "Heathland" = "#AA4499")) +
  
  scale_shape_manual(values = shape_vals) +
  
  scale_linetype_manual(values = line_vals) +
  
  scale_x_date(date_breaks = "3 days", date_labels = "%d %b") +
  
  theme_bw() +
  
  labs(
    x = "Date",
    y = expression("Mean CH"[4] * " flux (nmol" ~ s^{-1} ~ m^{-2} * ")"),
    colour = "Habitat",
    fill = "Habitat",
    shape = "Bracken",
    linetype = "Bracken"
  )
#show the figure
show(ch4_before)
#save the figure
#ggsave("Figures/co2_before.svg", plot = co2_before, width = 7.5, height = 6, dpi = 300)


#compare AICs
m1 <- glmmTMB(
  average_CH4_flux ~ ns(date_num, df = 1) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)
m2 <- glmmTMB(
  average_CH4_flux ~ ns(date_num, df = 2) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)

m3 <- glmmTMB(
  average_CH4_flux ~ ns(date_num, df = 3) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)

m4 <- glmmTMB(
  average_CH4_flux ~ ns(date_num, df = 4) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)

m5 <- glmmTMB(
  average_CH4_flux ~ ns(date_num, df = 5) * Habitat * Bracken,
  data = flux_before,
  family = gaussian()
)
#linear model bascially the same as the others, apart from overfitted m5, so just use linear models for this bit
AIC(m1, m2, m3, m4, m5)



hist(flux_after$average_CH4_flux)

#library(splines)
#splines allows for curved model
ch4_model <- glmmTMB(
  average_CH4_flux ~ ns(date_num, df = 1) * Habitat * Bracken,
  data = flux_after,
  family = gaussian()
)
# check residuals 
simulateResiduals(ch4_model, plot = TRUE)

#convert model slope to change in flux per day
ch4_slopes_after <- slopes(
  ch4_model,
  variables = "date_num",
  by = c("Habitat", "Bracken")
)


#create prediction data
newdatafter <- expand.grid(
  date_num = seq(
    min(flux_after$date_num),
    max(flux_after$date_num),
    length.out = 200
  ),
  Habitat = levels(flux_after$Habitat),
  Bracken = levels(flux_after$Bracken)
)

newdatafter$Habitat <- factor(
  newdatafter$Habitat,
  levels = levels(flux_after$Habitat)
)

newdatafter$Bracken <- factor(
  newdatafter$Bracken,
  levels = levels(flux_after$Bracken)
)

#use the model to predict
predafter <- predict(
  ch4_model,
  newdata = newdatafter,
  type = "response",
  se.fit = TRUE
)

newdatafter$fit <- predafter$fit
newdatafter$lower <- predafter$fit - 1.96 * predafter$se.fit
newdatafter$upper <- predafter$fit + 1.96 * predafter$se.fit

newdatafter$date <- as.Date(
  newdatafter$date_num,
  origin = "1970-01-01"
)

#create grouping variable
newdatafter$Group <- interaction(
  newdatafter$Habitat,
  newdatafter$Bracken,
  sep = "-"
)

flux_after$Group <- interaction(
  flux_after$Habitat,
  flux_after$Bracken,
  sep = "-"
)

#set shape and line formatting for the different treatemnts
shape_vals <- c("Absent" = 17,  # triangle
                "Present" = 16) # circle

line_vals <- c("Absent" = "dashed",
               "Present" = "solid")



ch4_after <- ggplot() +
  
  # observed data
  geom_point(
    data = flux_after,
    aes(
      x = date,
      y = average_CH4_flux,
      colour = Habitat,
      shape = Bracken
    ),
    alpha = 0.8,
    size = 2
  ) +
  
  # model ribbon
  geom_ribbon(
    data = newdatafter,
    aes(
      x = date,
      ymin = lower,
      ymax = upper,
      fill = Habitat,
      group = interaction(Habitat, Bracken)
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  # model lines
  geom_line(
    data = newdatafter,
    aes(
      x = date,
      y = fit,
      colour = Habitat,
      linetype = Bracken,
      group = interaction(Habitat, Bracken)
    ),
    linewidth = 1
  ) +
  coord_cartesian(ylim = c(-2.0, 0.5)) +
  
  # scales
  scale_colour_manual(values = c("Grassland" = "#1b9e77",
                                 "Heathland" = "#AA4499")) +
  
  scale_fill_manual(values = c("Grassland" = "#1b9e77",
                               "Heathland" = "#AA4499")) +
  
  scale_shape_manual(values = shape_vals) +
  
  scale_linetype_manual(values = line_vals) +
  
  scale_x_date(date_breaks = "3 days", date_labels = "%d %b") +
  
  theme_bw() +
  
  labs(
    x = "Date",
    y = NULL,
    colour = "Habitat",
    fill = "Habitat",
    shape = "Bracken",
    linetype = "Bracken"
  )
#show the figure
show(ch4_after)
#panel figure

#combine plots into a single figure
ch4_models <- ggarrange(ch4_before, ch4_after,
                        labels = c("A", "B"),
                        ncol = 2, nrow = 1,
                        #the width of each panel of the multifigure plot
                        widths = c(1,1.6),
                        legend = FALSE)

show(ch4_models)
#save the figure
#ggsave("Figures/ch4_models.svg", plot = ch4_models, width = 10, height = 6, dpi = 300)


#combine plots into a single figure
gasflux_models <- ggarrange(co2_before, co2_after, ch4_before, ch4_after,
                            labels = c("A", "B", "C", "D"),
                            ncol = 2, nrow = 2,
                            #the width of each panel of the multifigure plot
                            widths = c(1,1.6),
                            legend = FALSE)

show(gasflux_models)
#save the figure
ggsave("Figures/gasflux_models.svg", plot = gasflux_models, width = 10, height = 6, dpi = 300)



#### model differences after rainfall applied between rainfall treatments for CO2 (FIGURE 3) and CH4 (Figure 4) ----
flux_data <- read_csv("Processed Data/all_times_gas_fluxes.csv")
#remove data entries with p values < 0.05
#flux_data <- flux_data[flux_data$CO2_p_value <= 0.05,]
flux_timeseries <- flux_data[c(1,2, 12, 13, 32, 33, 34)]

#get the date only in one column
flux_timeseries <- flux_timeseries %>%
  mutate(date = as.Date(start_time))

#group by date, Habitat, and Bracken, and calculate means
flux_avg <- flux_timeseries %>%
  group_by(date, Habitat, Bracken, `Rainfall volume (ml)`) %>%
  summarise(
    `Average CO2 flux` = mean(`CO2 flux (micromol CO2 per s per m2)`, na.rm = TRUE),
    `SD CO2 flux` = sd(`CO2 flux (micromol CO2 per s per m2)`, na.rm = TRUE),
    `Average CH4 flux` = mean(`CH4 flux (nmol CH4 per s per m2)`, na.rm = TRUE),
    `SD CH4 flux` = sd(`CH4 flux (nmol CH4 per s per m2)`, na.rm = TRUE),
    .groups = 'drop'
  )

#model only prior to the rainfall treatment


# Create a combined grouping label
flux_avg <- flux_avg %>%
  mutate(Habitat_Bracken = paste(Habitat, Bracken, sep = " - "))

#rename the variables to avoid modelling issue with glmmTMB
colnames(flux_avg)[colnames(flux_avg) == "Average CO2 flux"] <- "average_CO2_flux"
#rename the variables to avoid modelling issue with glmmTMB
colnames(flux_avg)[colnames(flux_avg) == "Average CH4 flux"] <- "average_CH4_flux"
#rename the variables to avoid modelling issue with glmmTMB
colnames(flux_avg)[colnames(flux_avg) == "Rainfall volume (ml)"] <- "rainfall"


#convert dates to numbers
flux_avg$date <- as.Date(flux_avg$date)
flux_avg$date_num <- as.numeric(flux_avg$date)

flux_avg$Habitat <- factor(flux_avg$Habitat)
flux_avg$Bracken <- factor(flux_avg$Bracken)

#gas fluxes after rainfall treatment applied
flux_after <- flux_avg %>%
  filter(date > as.Date("2025-10-14"))

hist(flux_after$average_CO2_flux)

# ensure rainfall is a factor
flux_after$rainfall <- factor(flux_after$rainfall)

co2_after_rainfall_model <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 1)*Habitat*Bracken*rainfall,
  data = flux_after,
  family = gaussian()
)


#create prediction data
newdatafter <- expand.grid(
  date_num = seq(
    min(flux_after$date_num),
    max(flux_after$date_num),
    length.out = 200
  ),
  Habitat = levels(flux_after$Habitat),
  Bracken = levels(flux_after$Bracken),
  rainfall = unique(flux_after$rainfall)
)

newdatafter$Habitat <- factor(
  newdatafter$Habitat,
  levels = levels(flux_after$Habitat)
)
newdatafter$Bracken <- factor(
  newdatafter$Bracken,
  levels = levels(flux_after$Bracken)
)
newdatafter$rainfall <- factor(
  newdatafter$rainfall,
  levels = unique(flux_after$rainfall)
)

# use the model to predict
predafter <- predict(
  co2_after_rainfall_model,
  newdata = newdatafter,
  type = "response",
  se.fit = TRUE
)

newdatafter$fit <- predafter$fit
newdatafter$lower <- predafter$fit - 1.96 * predafter$se.fit
newdatafter$upper <- predafter$fit + 1.96 * predafter$se.fit
newdatafter$date <- as.Date(
  newdatafter$date_num,
  origin = "1970-01-01"
)

# create grouping variable (now includes rainfall)
newdatafter$Group <- interaction(
  newdatafter$Habitat,
  newdatafter$Bracken,
  newdatafter$rainfall,
  sep = "-"
)
flux_after$Group <- interaction(
  flux_after$Habitat,
  flux_after$Bracken,
  flux_after$rainfall,
  sep = "-"
)
#set shape and line formatting for the different treatemnts
shape_vals <- c("Absent" = 17,  # triangle
                "Present" = 16) # circle

line_vals <- c("Absent" = "dashed",
               "Present" = "solid")

#plot showing flux for each rainfall treatment
co2_after <- ggplot() + 
  # observed data
  geom_point(data = flux_after, aes(x = date, y = average_CO2_flux, colour = rainfall), alpha = 0.8, size = 2) +
    # model ribbon
  geom_ribbon(data = newdatafter, aes(x = date, ymin = lower, ymax = upper, fill = rainfall, group = rainfall), alpha = 0.15, colour = NA) +
    # model lines
  geom_line(data = newdatafter, aes(x = date, y = fit, colour = rainfall, group = rainfall), linewidth = 1) +
    # four panels: one per Habitat x Bracken combination
  facet_grid(Bracken ~ Habitat) +
    # scales
  scale_colour_viridis_d(name = "Rainfall (ml)", option = "plasma", direction = -1) +
  
  scale_fill_viridis_d(name = "Rainfall (ml)", option = "plasma", direction = -1) +
  
  scale_x_date(date_breaks = "3 days", date_labels = "%d %b") +
  
  theme_bw() +
  
  labs(x = "Date",
    y = expression("Mean CO"[2] * " flux (" * mu * "mol" ~ s^{-1} ~ m^{-2} * ")"),
    colour = "Rainfall (ml)",
    fill = "Rainfall (ml)"
  )
#show the figure
show(co2_after)
#save the figure
#ggsave("Figures/co2_rainfall_models.svg", plot = co2_after, width = 10, height = 6, dpi = 300)


#add for CH4
# Fit CH4 model
ch4_after_rainfall_model <- glmmTMB(
  average_CH4_flux ~ ns(date_num, df = 1)*Habitat*Bracken*rainfall,
  data = flux_after,
  family = gaussian()
)

# Predict from CH4 model
predafter_ch4 <- predict(
  ch4_after_rainfall_model,
  newdata = newdatafter,
  type = "response",
  se.fit = TRUE
)

# Store CH4 predictions in a separate dataframe
newdatafter_ch4 <- newdatafter
newdatafter_ch4$fit   <- predafter_ch4$fit
newdatafter_ch4$lower <- predafter_ch4$fit - 1.96 * predafter_ch4$se.fit
newdatafter_ch4$upper <- predafter_ch4$fit + 1.96 * predafter_ch4$se.fit

# Now plot using newdatafter_ch4
ch4_after <- ggplot() +
  #data points
  geom_point(
    data = flux_after,
    aes(x = date, y = average_CH4_flux, colour = rainfall),
    alpha = 0.8, size = 2
  ) +
  #the model ribbon
  geom_ribbon(
    data = newdatafter_ch4,
    aes(x = date, ymin = lower, ymax = upper,
        fill = rainfall, group = rainfall),
    alpha = 0.15, colour = NA
  ) +
  #model line
  geom_line(
    data = newdatafter_ch4,
    aes(x = date, y = fit, colour = rainfall, group = rainfall),
    linewidth = 1
  ) +
  
  facet_grid(Bracken ~ Habitat) +
  coord_cartesian(ylim = c(-3.0, 0.7)) +
  scale_colour_viridis_d(name = "Rainfall (ml)", option = "plasma", direction = -1) +
  scale_fill_viridis_d(name = "Rainfall (ml)", option = "plasma", direction = -1) +
  scale_x_date(date_breaks = "3 days", date_labels = "%d %b") +
  theme_bw() +
  labs(x = "Date", y = expression("Mean CH"[4] * " flux (nmol" ~ s^{-1} ~ m^{-2} * ")"), colour = "Rainfall (ml)", fill = "Rainfall (ml)")

show(ch4_after)
#save the figure
#ggsave("Figures/ch4_rainfall_models.svg", plot = ch4_after, width = 10, height = 6, dpi = 300)


# check to see if rainfaill intensity changes the model slope (i.e. co2 flux over time)
emtrends(co2_after_rainfall_model, 
         pairwise ~ rainfall | Habitat * Bracken, 
         var = "date_num")

# check to see if rainfaill intensity changes the model slope (i.e. ch4 flux over time)
emtrends(ch4_after_rainfall_model, 
         pairwise ~ rainfall | Habitat * Bracken, 
         var = "date_num")
#### model differences prior to rainfall addition for CO2 (supp FIGURE 1) and CH4 (supp Figure 2) NEEDS PRIOR TAB TO WORK ----
#gas fluxes after rainfall treatment applied
flux_before <- flux_avg %>%
  filter(date < as.Date("2025-10-14"))

hist(flux_before$average_CO2_flux)

# ensure rainfall is a factor
flux_before$rainfall <- factor(flux_before$rainfall)

co2_before_rainfall_model <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 1)*Habitat*Bracken*rainfall,
  data = flux_before,
  family = gaussian()
)


#create prediction data
newdatbefore <- expand.grid(
  date_num = seq(
    min(flux_before$date_num),
    max(flux_before$date_num),
    length.out = 200
  ),
  Habitat = levels(flux_before$Habitat),
  Bracken = levels(flux_before$Bracken),
  rainfall = unique(flux_before$rainfall)
)

newdatbefore$Habitat <- factor(
  newdatbefore$Habitat,
  levels = levels(flux_before$Habitat)
)
newdatbefore$Bracken <- factor(
  newdatbefore$Bracken,
  levels = levels(flux_before$Bracken)
)
newdatbefore$rainfall <- factor(
  newdatbefore$rainfall,
  levels = unique(flux_before$rainfall)
)

# use the model to predict
predbefore <- predict(
  co2_before_rainfall_model,
  newdata = newdatbefore,
  type = "response",
  se.fit = TRUE
)

newdatbefore$fit <- predbefore$fit
newdatbefore$lower <- predbefore$fit - 1.96 * predbefore$se.fit
newdatbefore$upper <- predbefore$fit + 1.96 * predbefore$se.fit
newdatbefore$date <- as.Date(
  newdatbefore$date_num,
  origin = "1970-01-01"
)

# create grouping variable (now includes rainfall)
newdatbefore$Group <- interaction(
  newdatbefore$Habitat,
  newdatbefore$Bracken,
  newdatbefore$rainfall,
  sep = "-"
)
flux_before$Group <- interaction(
  flux_before$Habitat,
  flux_before$Bracken,
  flux_before$rainfall,
  sep = "-"
)
#set shape and line formatting for the different treatemnts
shape_vals <- c("Absent" = 17,  # triangle
                "Present" = 16) # circle

line_vals <- c("Absent" = "dashed",
               "Present" = "solid")




co2_before <- ggplot() +
  
  # observed data
  geom_point(
    data = flux_before,
    aes(
      x = date,
      y = average_CO2_flux,
      colour = rainfall
    ),
    alpha = 0.8,
    size = 2
  ) +
  
  # model ribbon
  geom_ribbon(
    data = newdatbefore,
    aes(
      x = date,
      ymin = lower,
      ymax = upper,
      fill = rainfall,
      group = rainfall
    ),
    alpha = 0.15,
    colour = NA
  ) +
  
  # model lines
  geom_line(
    data = newdatbefore,
    aes(
      x = date,
      y = fit,
      colour = rainfall,
      group = rainfall
    ),
    linewidth = 1
  ) +
  
  # four panels: one per Habitat x Bracken combination
  facet_grid(Bracken ~ Habitat) +
  
  # scales
  scale_colour_viridis_d(
    name = "Rainfall (ml)",
    option = "plasma",
    direction = -1
  ) +
  
  scale_fill_viridis_d(
    name = "Rainfall (ml)",
    option = "plasma", direction = -1
  ) +
  
  scale_x_date(date_breaks = "3 days", date_labels = "%d %b") +
  
  theme_bw() +
  
  labs(
    x = "Date",
    y = expression("Mean CO"[2] * " flux (" * mu * "mol" ~ s^{-1} ~ m^{-2} * ")"),
    colour = "Rainfall (ml)",
    fill = "Rainfall (ml)"
  )
#show the figure
show(co2_before)
#save the figure
ggsave("Figures/before-rain-applied_co2_rainfall_models.svg", plot = co2_before, width = 10, height = 6, dpi = 300)


#add for CH4
# Fit CH4 model
ch4_before_rainfall_model <- glmmTMB(
  average_CH4_flux ~ ns(date_num, df = 1)*Habitat*Bracken*rainfall,
  data = flux_before,
  family = gaussian()
)

# Predict from CH4 model
predbefore_ch4 <- predict(
  ch4_before_rainfall_model,
  newdata = newdatbefore,
  type = "response",
  se.fit = TRUE
)

# Store CH4 predictions in a separate dataframe
newdatbefore_ch4 <- newdatbefore
newdatbefore_ch4$fit   <- predbefore_ch4$fit
newdatbefore_ch4$lower <- predbefore_ch4$fit - 1.96 * predbefore_ch4$se.fit
newdatbefore_ch4$upper <- predbefore_ch4$fit + 1.96 * predbefore_ch4$se.fit

# Now plot using newdatafter_ch4
ch4_before <- ggplot() +
  
  geom_point(
    data = flux_before,
    aes(x = date, y = average_CH4_flux, colour = rainfall),
    alpha = 0.8, size = 2
  ) +
  
  geom_ribbon(
    data = newdatbefore_ch4,
    aes(x = date, ymin = lower, ymax = upper,
        fill = rainfall, group = rainfall),
    alpha = 0.15, colour = NA
  ) +
  
  geom_line(
    data = newdatbefore_ch4,
    aes(x = date, y = fit, colour = rainfall, group = rainfall),
    linewidth = 1
  ) +
  
  facet_grid(Bracken ~ Habitat) +
  coord_cartesian(ylim = c(-3.0, 0.7)) +
  scale_colour_viridis_d(name = "Rainfall (ml)", option = "plasma", direction = -1) +
  scale_fill_viridis_d(name = "Rainfall (ml)", option = "plasma", direction = -1) +
  scale_x_date(date_breaks = "3 days", date_labels = "%d %b") +
  theme_bw() +
  labs(x = "Date", y = expression("Mean CH"[4] * " flux (nmol" ~ s^{-1} ~ m^{-2} * ")"), colour = "Rainfall (ml)", fill = "Rainfall (ml)")

show(ch4_before)
#save the figure
ggsave("Figures/before-rainfall-applied_ch4_rainfall_models.svg", plot = ch4_before, width = 10, height = 6, dpi = 300)

# check to see if rainfaill intensity changes the model slope (i.e. co2 flux over time)
emtrends(co2_before_rainfall_model, 
         pairwise ~ rainfall | Habitat * Bracken, 
         var = "date_num")

# check to see if rainfaill intensity changes the model slope (i.e. ch4 flux over time)
emtrends(ch4_before_rainfall_model, 
         pairwise ~ rainfall | Habitat * Bracken, 
         var = "date_num")


#### compare before and after models - NEEDS PRIOR TWO TABS TO WORK -   don't think we need this now, not included in papr results ----

# Create combined before/after dataframe
flux_combined <- bind_rows(
  flux_before %>% mutate(Period = "Before"),
  flux_after  %>% mutate(Period = "After")
) %>%
  mutate(Period = factor(Period, levels = c("Before", "After")))

#check data distribution
hist(flux_combined$average_CO2_flux)
# Models on combined data - Period replaces rainfall as the key fixed effect
co2_period_model <- glmmTMB(
  average_CO2_flux ~ ns(date_num, df = 1) * Habitat * Bracken * Period,
  data = flux_combined,
  family = gaussian()
)

# check residuals - looks great.  Well done me.
simulateResiduals(co2_period_model, plot = TRUE)

# Compare average flux before vs after, per Habitat x Bracken
emmeans(co2_period_model,
        pairwise ~ Period | Habitat * Bracken,
        at = list(date_num = mean(flux_combined$date_num)))

# Compare slope before vs after, per Habitat x Bracken
emtrends(co2_period_model,
         pairwise ~ Period | Habitat * Bracken,
         var = "date_num")




#check data distribution - skewed, so cannot use gaussian
hist(flux_combined$average_CH4_flux)

#try t_family...
ch4_period_model <- glmmTMB(
  average_CH4_flux ~ ns(date_num, df = 1) * Habitat * Bracken * Period,
  data = flux_combined,
  family = t_family()
)

# residuals perfect now we are using t_family()
simulateResiduals(ch4_period_model, plot = TRUE)

emmeans(ch4_period_model,
        pairwise ~ Period | Habitat * Bracken,
        at = list(date_num = mean(flux_combined$date_num)))



emtrends(ch4_period_model,
         pairwise ~ Period | Habitat * Bracken,
         var = "date_num")

