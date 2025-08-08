
# === 1. Packages ####
rm(list=ls()) 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("phyloseq") 

#remotes::install_github("david-barnett/microViz")

suppressPackageStartupMessages({
  library(devtools)
  library(phyloseq)
  library(multcomp)
  library(microViz)
  library(vegan)
  library(ggpubr)
  library(dplyr)
  library(sf)
  library(geosphere)
  library(leaflet)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ggplot2)
  library(ggthemes)
  library(tinytex)
  library(dplyr)
  library(emmeans)
  library(performance)
  library(lmerTest)
  library(lme4)
  library(sf)
  library(raster)  
  library(terra)
  library(car)
  library(ggeffects)
  library(dplyr)
  library(broom.mixed)
  library(purrr)
  library(tibble)
})

#--- 2. Data transformation --- ####

# Load your phyloseq object here 
#pw <- readRDS("Data/MiCoDaV2_PrimerV4_Water_OnlyOcean.rds")
ps <- readRDS("data/physeq_subset_soil_Filtered2025_updated.rds")

# Metadata

names(ps@sam_data)

# 4. Plot the map
# Base world map
metadata <- ps@sam_data


df.20 <- metadata[1:729,]

# Aggregate data to count samples per locality (longitude and latitude)
df.aggregated <- df.20 %>%
  group_by(Longitude, Latitude, environment_.biome.) %>%
  summarise(sample_count = n(), .groups = "drop")

# Convert to an sf object for spatial plotting
df.sf <- st_as_sf(df.aggregated, coords = c("Longitude", "Latitude"), crs = 4326)



#--- 3. Köppen Climate Zones --- ####

library(raster)
library(sf)
library(dplyr)

# Load Köppen climate shapefile
# Replace the path with the actual folder where you saved it
koppen_raster <- raster("Data/1991_2020/koppen_geiger_0p5.tif")  # adjust filename if needed

# Make sure CRS matches
df.sf <- st_transform(df.sf, crs(koppen_raster))

# Extract Köppen zone values
koppen_values <- raster::extract(koppen_raster, df.sf)

# Add to your spatial data
df.sf$Koppen_Zone <- koppen_values

metadata$Koppen_Zone <- df.sf$Koppen_Zone
sample_data(ps)$Koppen_Zone <- metadata$Koppen_Zone


#geom_sf(data = df.sf, aes(color = Koppen_Zone))

#match to actual zone names
koppen_lookup <- c(
  "1" = "Af", "2" = "Am", "3" = "Aw",  # tropical
  "4" = "BWh", "5" = "BWk", "6" = "BSh", "7" = "BSk",  # arid
  "8" = "Csa", "9" = "Csb", "10" = "Csc",  # temperate
  "11" = "Cwa", "12" = "Cwb", "13" = "Cwc",
  "14" = "Cfa", "15" = "Cfb", "16" = "Cfc",
  "17" = "Dsa", "18" = "Dsb", "19" = "Dsc", "20" = "Dsd",
  "21" = "Dwa", "22" = "Dwb", "23" = "Dwc", "24" = "Dwd",
  "25" = "Dfa", "26" = "Dfb", "27" = "Dfc", "28" = "Dfd",
  "29" = "ET", "30" = "EF"
)

# Assign readable labels
df.sf$Koppen_Label <- koppen_lookup[as.character(df.sf$Koppen_Zone)]

metadata$Koppen_Label <- df.sf$Koppen_Label
sample_data(ps)$Koppen_Label <- metadata$Koppen_Label


#Broader Köppen Zones 

# 1. Define broad zone mapping with full names
koppen_broad_lookup <- c(
  # Tropical
  "Af" = "Tropical", "Am" = "Tropical", "Aw" = "Tropical",
  
  # Arid
  "BWh" = "Arid", "BWk" = "Arid", "BSh" = "Arid", "BSk" = "Arid",
  
  # Temperate
  "Csa" = "Temperate", "Csb" = "Temperate", "Csc" = "Temperate",
  "Cwa" = "Temperate", "Cwb" = "Temperate", "Cwc" = "Temperate",
  "Cfa" = "Temperate", "Cfb" = "Temperate", "Cfc" = "Temperate",
  
  # Cold (formerly "D")
  "Dsa" = "Cold", "Dsb" = "Cold", "Dsc" = "Cold", "Dsd" = "Cold",
  "Dwa" = "Cold", "Dwb" = "Cold", "Dwc" = "Cold", "Dwd" = "Cold",
  "Dfa" = "Cold", "Dfb" = "Cold", "Dfc" = "Cold", "Dfd" = "Cold",
  
  # Polar
  "ET" = "Polar", "EF" = "Polar"
)

# 2. Assign broad climate zone
df.sf$Koppen_Broad <- koppen_broad_lookup[df.sf$Koppen_Label]

# 3. Add to metadata and phyloseq object
metadata$Koppen_Broad <- df.sf$Koppen_Broad
sample_data(ps)$Koppen_Broad <- metadata$Koppen_Broad

#Prepare df for plotting Koppen Map
koppen_df <- as.data.frame(koppen_raster, xy = TRUE)
koppen_df <- na.omit(koppen_df)
koppen_df$zone_code <- as.character(koppen_df$koppen_geiger_0p5)  # adjust column name accordingly
koppen_df$zone_label <- koppen_lookup[koppen_df$zone_code]
koppen_df$zone_broad <- koppen_broad_lookup[koppen_df$zone_label]


# --- 4. Shannon Index --- ####

# Alpha Diversity and Richness
# Calculate diversity (Observed richness and Shannon/Jaccard diversity)
alpha_div <- estimate_richness(ps, measures = c("Observed", "Shannon")) 
div_measures <- alpha_div
div_measures$SampleID <- rownames(div_measures) 

# Extract metadatas
metadata <- data.frame(sample_data(ps))
metadata$SampleID <- rownames(metadata)

# Merge diversity and metadata
data_plot <- merge(div_measures, metadata, by = "SampleID")

#Create separate lat/long columns
#data_lat_lon <- data_plot_w %>%
#  separate(lat_long, into = c("Latitude", "Longitude"), sep = ",", convert = TRUE)

# Create Absolute Latitude column
data_plot$Abs_Latitude <- abs(data_plot$Latitude)




# --- 5. World Map --- ####

# Load world map with gray background
world <- ne_countries(scale = "medium", returnclass = "sf")


# Define color mapping with a colorblind-friendly palette
color_map <- c(
  "Water" = "#377eb8",    # DodgerBlue
  "Organism" = "#984ea3", # DarkOrange
  "Other" = "#4daf4a",    # LimeGreen
  "Mineral" = "#ff7f00"   # Crimson
) 

# Create the first map --> data distribution map

ggplot() +
  geom_sf(data = world, fill = "gray90", color = "gray70") +
  geom_sf(data = df.sf, color = "#4DAF4A", alpha = 0.7, size = 3) +
  coord_sf(
    xlim = c(-170, -30),   # longitude range for Americas
    ylim = c(-60, 80),     # latitude range to include N and S America
    expand = FALSE
  ) +
  theme_minimal() +
  labs(
    title = "MiCoDa Version 2.0",
    subtitle = "Primer V4, Soil"
  ) +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12)
  )

  

#koppen zones map
ggplot() +
  geom_tile(data = koppen_df, aes(x = x, y = y, fill = zone_label)) +
  scale_fill_viridis_d(
    name = "Köppen Climate Zone",
    option = "viridis",
    direction = -1  # reverses the color order
  ) +
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Global Köppen-Geiger Climate Zones",
    subtitle = "Based on 1991–2020 data"
  ) +
  theme(
    legend.position = "right",         # place legend next to the map
    legend.key.size = unit(0.8, "lines"), # smaller legend keys
    legend.text = element_text(size = 8), # smaller legend text
    legend.title = element_text(size = 10), # smaller legend title
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )


#Add Shannon Index data points
ggplot() +
  # Köppen climate background
  geom_tile(data = koppen_df, aes(x = x, y = y, fill = zone_broad)) +
  scale_fill_manual(
    name = "Köppen Climate Zone",
    values = c(
      "Polar" = "#cccccc",  
      "Cold" = "#999999",      
      "Temperate" = "#777777", 
      "Arid" = "#555555",      
      "Tropical" = "#444444"      
    )
  ) +
  
  # Shannon diversity points
  # Fill layer
  geom_point(
    data = data_plot,
    aes(x = Longitude, y = Latitude, size = Shannon),
    color = "#4DAF4A",
    fill = "#4DAF4A", # required for shape 21
    alpha = 0.7,
    shape = 21
  )+
  scale_size_continuous(
    name = "Shannon Diversity",
    range = c(0.8, 5)  # Adjust size range as needed
  ) +
  
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Köppen Climate Zones with Shannon Diversity",
    subtitle = "MiCoDa V2 - Primer V4 (Soil)",
    x = "",
    y = ""
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )


# Shannon Diversity Analysis by Continent
ggplot(data_plot, aes(x = Continent, y = Shannon, fill = Continent)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # no outliers, more transparent boxes
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +  # points for all samples
  theme_minimal() +
  labs(
    title = "Shannon Diversity by Continent",
    y = "Shannon Index",
    x = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_brewer(palette = "Set1")  # optional: nice color palette


# Shannon Diversity Analysis by broad Climate Zone
data_plot$Koppen_Broad <- metadata$Koppen_Broad
# Remove rows where Koppen_Broad is NA
data_plot <- data_plot[!is.na(data_plot$Koppen_Broad), ]

ggplot(data_plot, aes(x = Koppen_Broad, y = Shannon, fill = Koppen_Broad)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +  # no outliers, more transparent boxes
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +  # points for all samples
  theme_minimal() +
  labs(
    title = "Shannon Diversity by Broad Climate Zones",
    y = "Shannon Index",
    x = "Köppen Climate Zone"
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_brewer(palette = "Set1")  # Color palette for discrete variables



# --- 6. Abiotic Predictors --- ####

# * _ PH Global ####

library(terra)
library(sf)
library(ggplot2)
library(viridis)

# --- Load pH raster ---
ph_raster <- rast("Data/ph_soil_30cm.tif")

# --- Prepare spatial sample points ---
# Assuming your sample data with lat/lon is in 'data_plot'

# Convert data_plot to sf object (if not already)
df.sf <- st_as_sf(data_plot, coords = c("Longitude", "Latitude"), crs = 4326)

# Reproject sample points to raster CRS
df.sf <- st_transform(df.sf, crs(ph_raster))

# Convert sf to SpatVector for terra
df_vect <- vect(df.sf)

# --- Extract pH values from raster at sample points ---
ph_values <- terra::extract(ph_raster, df_vect, ID=FALSE)[, 1]  # extract values, drop ID

# Add extracted pH values to sf object
df.sf$pH <- ph_values

# --- Handle scaling ---
# If raster values are stored as integer * 10 (check this for your raster!)
df.sf$pH <- df.sf$pH / 10

# Also add pH to original dataframe for modeling
data_plot$pH_from_raster <- ph_values / 10
sample_data(ps)$pH_from_raster <- metadata$pH_from_raster

# --- Prepare raster for plotting ---
# Aggregate raster to lower resolution for faster plotting (factor 40 is example)
ph_raster_lowres <- terra::aggregate(ph_raster, fact = 40, fun = mean)

# Convert to dataframe for ggplot
ph_df <- as.data.frame(ph_raster_lowres, xy = TRUE, na.rm = TRUE)
colnames(ph_df) <- c("x", "y", "pH")

# Scale raster pH values as well
ph_df$pH <- ph_df$pH / 10

# --- Plotting ---
x11()
ggplot() +
  geom_tile(data = ph_df, aes(x = x, y = y, fill = pH)) +  # raster layer
  geom_sf(data = world, fill = NA, color = "gray70", size = 0.3) + # borders
  # geom_sf(data = df.sf, aes(color = "#4DAF4A"), size = 2) +  # points colored by extracted pH
  scale_fill_viridis_c(name = "Soil pH (Raster)", option = "plasma", direction = -1) +
  #scale_color_viridis_c(name = "Soil pH (Points)", option = "plasma", direction = -1) +
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Global Soil pH (0–30 cm)",
    subtitle = "Source: SoilGrids",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )


ggplot() +
  geom_tile(data = ph_df, aes(x = x, y = y, fill = pH)) +  # raster layer
  geom_sf(data = world, fill = NA, color = "gray70", size = 0.3) + # borders
  geom_sf(data = cluster_sf, aes(), color = "darkgreen", size = 2, alpha = 0.7) +
  geom_hline(yintercept = 9, color = "red", linetype = "dashed", size = 1) +  # line at 9° latitude
  scale_fill_viridis_c(name = "Soil pH (Raster)", option = "plasma", direction = -1) +
  coord_sf(xlim = c(-170, -30), ylim = c(-60, 75), expand = FALSE) +  # zoom to Americas
  theme_minimal() +
  labs(
    title = "Soil pH in the Americas (0–30 cm)",
    subtitle = "Source: SoilGrids",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )


# --- Check for missing extracted values ---
summary(data_plot$pH_from_raster)
cat("Number of samples outside raster extent or missing values:", sum(is.na(data_plot$pH_from_raster)), "\n")

# Optional: list samples with missing pH
missing_samples <- data_plot[is.na(data_plot$pH_from_raster), ]
if(nrow(missing_samples) > 0) {
  print(missing_samples[, c("SampleID", "Latitude", "Longitude")])
}



# * _ Soil Temperature 0-5 cm ####

# Load your temp raster (e.g., 0–5 cm depth layer)
temp_raster <- rast("Data/SBIO1_Annual_Mean_Temperature_0_5cm.tif")


# Reproject sample data to match pH raster
df.sf <- st_transform(df.sf, crs(temp_raster))


# Convert to 'SpatVector' for terra extraction
df_vect <- vect(df.sf)

# Extract temp values
temp_values <- terra::extract(temp_raster, df_vect)[,2]  # [,2] to drop ID column

# Add to sf object
df.sf$mean_temp <- temp_values

# Add to metadata and phyloseq object
metadata$mean_temp <- temp_values
sample_data(ps)$mean_temp <- metadata$mean_temp

#plot global mean temp
# Convert the raster to a dataframe for plotting
# Reduce resolution by a factor of 40
temp_raster_lowres <- terra::aggregate(temp_raster, fact = 40, fun = mean)

# Now convert
temp_df <- as.data.frame(temp_raster_lowres, xy = TRUE, na.rm = TRUE)
colnames(temp_df) <- c("x", "y", "mean_temp")

# Remove NA values
temp_df <- na.omit(temp_df) 


# Plot global soil temp
  ggplot() +
  geom_tile(data = temp_df, aes(x = x, y = y, fill = mean_temp)) +  # raster layer
  geom_sf(data = world, fill = NA, color = "gray70", size = 0.3) +  # borders
  geom_sf(data = df.sf, color = "#4DAF4A", size = 2) +              # fixed green color_
  
  scale_fill_viridis_c(
    name = "Soil Temperature (Raster)",
    option = "plasma",
    direction = 1
  ) +
  scale_color_viridis_c(
    name = "Soil Temperature (Points)",
    option = "plasma",
    direction = 1
  ) +
  coord_sf() +   
 # geom_point(data = data_plot, aes(x = Longitude, y = Latitude, size = Shannon),color = "#4DAF4A", alpha = 0.7 ) +
  theme_minimal() +
  labs(
    title = "Global Soil Temperature (0–5 cm)",
    subtitle = "Source: ",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank()
  )

#prepare temp data for modeling
#read csv
  dat_Soil_2 <- read.csv2("Data/dat_Soil_2.csv")
  # Rename columns
  names(dat_Soil_2)[names(dat_Soil_2) == "Run"] <- "SampleID"
  names(dat_Soil_2)[names(dat_Soil_2) == "SBIO1_Annual_Mean_Temperature_0_5cm"] <- "mean_temp"
  
  # Join datasets
  data_plot <- left_join(data_plot, dat_Soil_2, by = "SampleID")  
  
  data_plot$mean_temp <- as.numeric(as.character(data_plot$mean_temp))
  dat_Soil_2$mean_temp <- as.numeric(as.character(dat_Soil_2$mean_temp))
  
  
  

# --- 7. Balancing Data --- ####
  
# * _ Continents ####

#data preparation
str(data_plot$Koppen_Broad)

data_plot %>%
  count(Continent.x)

#Exclude "Antarctica", "Oceania", "Europe" due to sparse data
data_continent_filtered <- data_plot %>%
  filter(!Continent.x %in% c("Antarctica", "Oceania", "Europe", "Africa", "Central America"),
         !is.na(Continent.x))

data_continent_filtered %>%
  count(Continent.x)

# Find the minimum group size
min_n <- data_continent_filtered %>%
  count(Continent.x) %>%
  summarise(min_n = min(n)) %>%
  pull(min_n)

# Before downsampling continents
set.seed(123)  # use any fixed number you like
continent_balanced <- data_continent_filtered %>%
  group_by(Continent.x) %>%
  slice_sample(n = min_n) %>%
  ungroup()

# Downsample each group to the same size
continent_balanced <- data_continent_filtered %>%
  group_by(Continent.x) %>%
  slice_sample(n = min_n) %>%
  ungroup()

continent_balanced %>%
  count(Continent.x)


# * _ Koppen Zones ####

data_plot %>%
  count(data_plot$Koppen_Broad)

continent_balanced %>%
  count(continent_balanced$Koppen_Broad)

#Exclude low sample zones
data_koppen_filtered <- continent_balanced %>%
  filter(!Koppen_Broad %in% c("Arid", "Polar"),   # exclude Arid and Polar
         !is.na(Koppen_Broad))                    # exclude NA values in Koppen_Broad

#balance climate zones
data_koppen_filtered %>%
  count(Koppen_Broad)

# Find the minimum group size
min_n <- data_koppen_filtered %>%
  count(Koppen_Broad) %>%
  summarise(min_n = min(n)) %>%
  pull(min_n)

# Before downsampling Koppen zones
set.seed(123)  # same seed ensures consistent selection
global_balanced <- data_koppen_filtered %>%
  group_by(Koppen_Broad) %>%
  slice_sample(n = min_n) %>%
  ungroup()

# Downsample each group to the same size
global_balanced <- data_koppen_filtered %>%
  group_by(Koppen_Broad) %>%
  slice_sample(n = min_n) %>%
  ungroup()

global_balanced %>%
  count(Koppen_Broad)



# * _ America ####

#get America data only
america_filtered <- data_plot %>%
  filter(
    Continent.x %in% c("North America", "Central America", "South America"),
    Latitude >= -60, Latitude <= 85,          # roughly Americas latitude range
    Longitude >= -170, Longitude <= -30       # roughly Americas longitude range
  )

#Set border between north and south america
america_panama <- america_filtered %>%
  mutate(Continent.x = if_else(
    Latitude >= 9 | (Latitude < 9 & Longitude < -79.9),
    "North America",
    "South America"
  ))


america_panama %>%
  count(Continent.x)

# Find the minimum group size
min_n <- america_panama %>%
  count(Continent.x) %>%
  summarise(min_n = min(n)) %>%
  pull(min_n)

# Downsample each group to the same size
america_balanced <- america_panama %>%
  group_by(Continent.x) %>%
  slice_sample(n = min_n) %>%
  ungroup()

america_balanced %>%
  count(Continent.x)



# --- 8. Modeling --- ####

# * _ Global Multiple Regression _ ####

# --- LMM global ---
global_model <- lmer(Shannon ~ pH_from_raster + mean_temp + Koppen_Broad + Continent.x + (1|City_Region.x) + (1|Authors.x), data = global_balanced)

#check for multicollinearity
vif(global_model)

#diagnostic plots
check_model(global_model, check = c("homogeneity", "linearity", "outliers", "normality"))

# a) Calculate Cook's Distance
cooksd <- cooks.distance(global_model)

# Use number of rows from the *data*, not the model object
n <- nrow(global_balanced)

# Define threshold
threshold <- 4 / n

# Identify influential points
influential <- which(cooksd > threshold)
length(influential)  # How many influential points

# Remove influential points from the original data
global_cleaned <- global_balanced[-influential, ]

                  
# --- refit cleaned ---
global_model_cleaned <- lmer(Shannon ~ pH_from_raster + mean_temp + Koppen_Broad + Continent.x + (1|City_Region.x) + (1|Authors.x), data = global_cleaned)

#check for multicollinearity
vif(global_model_cleaned)

#diagnostic plots
check_model(global_model_cleaned, check = c("homogeneity", "linearity", "outliers", "normality"))

summary(global_model_cleaned)

#refit with mean_temp and Koppen_Broad
global_model_re <- lmer(Shannon ~ mean_temp + Koppen_Broad + (1|City_Region.x) + (1|Authors.x), data = global_cleaned)

#check for multicollinearity
vif(global_model_re)

#diagnostic plots
check_model(global_model_re, check = c("homogeneity", "linearity", "outliers", "normality"))

summary(global_model_re)
r2(global_model_re)



# * __ Plot Model Results __ ####

#get predictions
pred_temp_gl <- ggpredict(global_model_re, terms = "mean_temp")
pred_koppen_gl <- ggpredict(global_model_re, terms = "Koppen_Broad")


#plot global temperature
ggplot() +
  geom_ribbon(data = pred_temp_gl, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "#CBB3E6", alpha = 0.25) +
  geom_point(data = global_balanced, aes(x = mean_temp, y = Shannon, color = Koppen_Broad),
             size = 2, alpha = 0.7, shape = 16) +
  geom_line(data = pred_temp_gl, aes(x = x, y = predicted), color = "black", size = 1.3) +
  labs(
    title = "Effect of Mean Temperature on Shannon Index",
    subtitle = "p = 1.14e-06 *** \nConditional R2: 0.952
     Marginal R2: 0.036",
    x = "Mean Temperature",
    y = "Predicted Shannon Index",
    color = "Climate Zone"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.key.size = unit(1, "lines")
  )

#plot global pH
ggplot() +
  geom_ribbon(data = pred_pH_gl, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "#CBB3E6", alpha = 0.25) +
  geom_point(data = global_balanced, aes(x = pH_from_raster, y = Shannon, color = Continent.x),
             size = 2, alpha = 0.7, shape = 16) +
  geom_line(data = pred_pH_gl, aes(x = x, y = predicted), color = "black", size = 1.3) +
  labs(
    title = "Effect of Mean Temperature on Shannon Index",
    subtitle = "p = 0.09182 . \nConditional R2: 0.879
     Marginal R2: 0.018",
    x = "pH",
    y = "Predicted Shannon Index (sqrt)",
    color = "Continent"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.key.size = unit(1, "lines")
  )



### --- Look at Interaction of mean_temp and Koppen_Broad ---

global_koppen <- lmer(Shannon ~ Koppen_Broad*mean_temp + (1|City_Region.x) + (1|Authors.x), data = global_balanced)
#diagnostic plots
check_model(global_koppen, check = c("homogeneity", "linearity", "outliers", "normality"))
### Remove influential outliers
cooksd <- cooks.distance(global_koppen)
# Threshold for influence
n <- nrow(global_balanced)
threshold <- 4 / n

# Find influential observations
influential_points <- which(cooksd > threshold)

# Remove and refit
cleaned_data <- global_balanced[-influential_points, ]
global_koppen_clean <- update(global_koppen, data = cleaned_data)

summary(global_koppen_clean)
anova(global_koppen_clean)
r2(global_koppen_clean)

emtrends(global_koppen, pairwise ~ Koppen_Broad, var = "mean_temp")


#Plot interaction effect
# Predict interaction between Koppen_Broad and mean_temp
interaction_preds <- ggpredict(global_koppen_clean, terms = c("mean_temp", "Koppen_Broad"))

ggplot(interaction_preds, aes(x = x, y = predicted, color = group)) +
  geom_line(size = 1) +
 # geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group), alpha = 0.2, color = NA) +
  labs(
    title = "Interaction: Climate Zone × Temperature",
    subtitle = "
Climate Zone             p= 0.003713 ** \n
Tempeature               p= 6.824e-10 ***\n
Climate Zone:Tempeature p= 0.003719**",
    x = "Mean Temperature",
    y = "Predicted Shannon Index",
    color = "Climate Zone",
    fill = "Climate Zone"
  ) +
  theme_minimal(base_size = 13)


# --- ANOVA of Climate Zones ---

global_koppen <- lmer(Shannon ~ Koppen_Broad + mean_temp + (1|City_Region.x) + (1|Authors.x), data = global_balanced)
#diagnostic plots
check_model(global_koppen, check = c("homogeneity", "linearity", "outliers", "normality"))
### Remove influential outliers
cooksd <- cooks.distance(global_koppen)
# Threshold for influence
n <- nrow(global_balanced)
threshold <- 4 / n

# Find influential observations
influential_points <- which(cooksd > threshold)

# Remove and refit
cleaned_data <- global_balanced[-influential_points, ]
global_koppen_clean <- update(global_koppen, data = cleaned_data)

summary(global_koppen_clean)
anova(global_koppen_clean)
r2(global_koppen_clean)
emm <- emmeans(global_koppen_clean, ~ Koppen_Broad)

# Pairwise comparisons with Tukey adjustment
posthoc_results <- pairs(emm, adjust = "tukey")

print(posthoc_results)

# Add letters to each mean
(model_means <- emmeans(object = global_koppen_clean ,
                        specs = ~ Koppen_Broad)) 

## Main effect of Type 
(model_means_cld <- cld(object = model_means,
                        adjust = "sidak", #applies the Sid?k correction, which adjusts p-values to reduce the likelihood of false positives when performing multiple pairwise comparisons
                        reversed=T, # makes "a" the letter for the highest mean
                        Letter = letters, #Uses lowercase letters (a, b, c, ...) instead of default uppercase
                        alpha = 0.05)) #Sets the significance level at 5% for determining group differences.

model_means_cld_df <- as.data.frame(model_means_cld)

# Boxplot
ggplot(cleaned_data, aes(x = Koppen_Broad, y = Shannon, fill = Koppen_Broad)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Koppen_Broad),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.4, size = 1.5) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  geom_text(
    data = model_means_cld_df,
    aes(x = Koppen_Broad, y = 8, label = .group),
    size = 6,
    inherit.aes = FALSE
  ) +
  labs(title = "Shannon-Index predicted by Climate Zone",
       subtitle = "p = 2.2e-16 ***  F = 89.24\nConditional R²: 0.989   Marginal R²: 0.034",
       x = "Climate Zone",
       y = "Shannon Index") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none" 
  )



# *  _ America Multiple Regression _ ####

#use LMM

#Random: City_Region.x , Authors.x , 

america_balanced$Shannon_sqrt <- sqrt(america_balanced$Shannon)
table(america_balanced$Shannon)

length(america_balanced)

ph_america <- lmer(Shannon_sqrt ~ pH_from_raster + mean_temp + Koppen_Broad + Continent.x + (1|City_Region.x) + (1|Authors.x), data = america_balanced)

#check for multicollinearity
vif(ph_america)

#diagnostic plots
check_model(ph_america, check = c("homogeneity", "linearity", "outliers", "normality"))

#a)
#Calculate Cook's Distance
cooksd <- cooks.distance(ph_america)
n <- nrow(america_balanced)
threshold <- 4/n
#threshold <- 1
influential <- which(cooksd > threshold)
length(influential)  # How many influential points
america_no_outliers <- america_balanced[-influential, ]

#b)
#filtering outliers
america_no_outliers <- america_no_outliers %>%
  filter(Shannon_sqrt >= 1)

#refit model
ph_america_clean <- lmer(Shannon_sqrt ~ pH_from_raster + mean_temp + Koppen_Broad + Continent.x + (1|City_Region.x) + (1|Authors.x), data = america_no_outliers)

#diagnostic plots 
check_model(ph_america_clean, check = c("homogeneity", "linearity", "outliers", "normality"))

summary(ph_america_clean)

plot(ph_america)

#remove mean_temp 
#refit model
regression_model_a <- lmer(Shannon_sqrt ~ pH_from_raster + Continent.x + (1|City_Region.x) + (1|Authors.x), data = america_no_outliers)

check_model(regression_model_a, check = c("homogeneity", "linearity", "outliers", "normality"))


summary(regression_model_a)
r2(regression_model_a)

#plot regression

# Add predicted values to your dataframe
# Get predicted values
library(ggeffects)
pred_ph <- ggpredict(regression_model_a, terms = "pH_from_raster")
pred_continent <- ggpredict(regression_model_a, terms = "Continent.x")



# * __ Plot Model Results __ ####

#plot ph america
ggplot() +
  geom_point(data = america_no_outliers, aes(x = pH_from_raster, y = Shannon_sqrt, color = Continent.x),
             size = 2, alpha = 0.7, shape = 16) +
  geom_ribbon(data = pred_ph, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "#CBB3E6", alpha = 0.25) +
  geom_line(data = pred_ph, aes(x = x, y = predicted), color = "black", size = 1.3) +
  scale_color_manual(values = c(
    "South America" = "#89CFF0",
    "North America" = "#FFB6B9"
  )) +
  labs(
    title = "Effect of Soil pH on Shannon Index (sqrt-transformed)",
    subtitle = "p = 0.01931 * \nConditional R2: 0.913
     Marginal R2: 0.029",
    x = "Soil pH",
    y = "Predicted Shannon Index (sqrt)",
    color = "Continent"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.key.size = unit(1, "lines")
  )


# Boxplot Koppen_Broad America
ggplot(america_no_outliers, aes(x = Continent.x, y = Shannon_sqrt, fill = Continent.x)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Continent.x),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.4, size = 1.5) +
  scale_fill_brewer(palette = "Set2") +
  scale_color_brewer(palette = "Set2") +
  labs(title = "Shannon-Index predicted by Continent",
       subtitle = "p = 0.00386 ** ***  Estimate = -0.10830\nConditional R2: 0.913
     Marginal R2: 0.029",
       x = "",
       y = "Shannon Index (sqrt-transformed)") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.position = "none" 
  )



# * __ Cluster Argentinia __ ####
cluster_points <- america_balanced %>%
  filter(pH_from_raster >= 6, pH_from_raster <= 7.5,
         Shannon_sqrt >= 1.5, Shannon_sqrt <= 2)

table(cluster_points$Authors.x)
cluster_points %>%
  count(Authors.x) %>%
  arrange(desc(n))

cluster_sf <- st_as_sf(cluster_points, coords = c("Longitude", "Latitude"), crs = 4326)

table(cluster_sf$VegetationType_LandUse.x) 

ggplot(data = world) +
  geom_sf(fill = "gray90", color = "gray50") +
  geom_sf(data = cluster_sf, aes(), color = "red", size = 2, alpha = 0.7) +
  theme_minimal() +
  labs(title = "Data points with pH 6-7.5 and Shannon_sqrt 1.5-2",
       subtitle = "Filtered cluster points",
       x = "Longitude",
       y = "Latitude")

#plot with pH
#x11()
ggplot() +
  geom_tile(data = ph_df, aes(x = x, y = y, fill = pH)) +  # raster layer
  geom_sf(data = world, fill = NA, color = "gray70", size = 0.3) + # borders
  geom_sf(data = cluster_sf, aes(), color = "green", size = 3, alpha = 1) +
  scale_fill_viridis_c(name = "Soil pH (Raster)", option = "plasma", direction = -1) +
  #scale_color_viridis_c(name = "Soil pH (Points)", option = "plasma", direction = -1) +
  coord_sf(xlim = c(-170, -30), ylim = c(-60, 75), expand = FALSE) +  # zoom to Americas
  theme_minimal() +
  labs(
    title = "Global Soil pH (0–30 cm)",
    subtitle = "Source: SoilGrids",
    x = "",
    y = ""
  ) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none",
  )



# --- 9. Distance Decay America --- ####

# Extract coordinates
coords <- dplyr::select(america_balanced, SampleID, Latitude, Longitude)

# Compute geographic distance matrix (Haversine distance)
geo_dist <- distm(coords[, c("Longitude", "Latitude")], fun = distHaversine)
rownames(geo_dist) <- coords$SampleID
colnames(geo_dist) <- coords$SampleID

# Bray-curtis distances
bray_dist <- readRDS("Data/Soil_bray_dist.rds")

# Convert Bray-Curtis distance to matrix
comm_dist <- as.matrix(bray_dist)

# Vectorize both matrices
geo_vec <- as.vector(geo_dist[lower.tri(geo_dist)])
comm_vec <- as.vector(comm_dist[lower.tri(comm_dist)])

# Check sample names in both distance matrices
samples_geo <- rownames(geo_dist)
samples_comm <- rownames(comm_dist)

# Find common samples
common_samples <- intersect(samples_geo, samples_comm)

# Subset both matrices to common samples and in the same order
geo_dist_sub <- geo_dist[common_samples, common_samples]
comm_dist_sub <- comm_dist[common_samples, common_samples]

# Now vectorize lower triangles
geo_vec <- as.vector(geo_dist_sub[lower.tri(geo_dist_sub)])
comm_vec <- as.vector(comm_dist_sub[lower.tri(comm_dist_sub)])

# Create dataframe
decay_df <- data.frame(
  GeoDistance = geo_vec / 1000,  # km
  CommunityDistance = comm_vec
)


# Plot distance decay relationship
ggplot(decay_df, aes(x = GeoDistance, y = CommunityDistance)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  theme_minimal() +
  labs(title = "Distance Decay",
       x = "Geographic Distance (km)",
       y = "Community Dissimilarity (Bray-Curtis)")



# Create all pair combinations of sample IDs
pairs <- t(combn(coords$SampleID, 2))

# Build pairs_df with sample pairs and their distances
pairs_df <- data.frame(
  Sample1 = pairs[,1],
  Sample2 = pairs[,2],
  GeoDistance = geo_vec,
  CommunityDistance = comm_vec
)

continent_map <- setNames(america_balanced$Continent.x, america_balanced$SampleID)
pairs_df$Continent1 <- continent_map[pairs_df$Sample1]
pairs_df$Continent2 <- continent_map[pairs_df$Sample2]

pairs_df$ContinentComparison <- apply(
  pairs_df[, c("Continent1", "Continent2")],
  1,
  function(x) paste(sort(x), collapse = "–")  # ensures "North–South", not "South–North"
)



ggplot(pairs_df, aes(x = GeoDistance / 1000, y = CommunityDistance, color = ContinentComparison)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm", se = FALSE, color = "red") +
  scale_color_manual(
    values = c(
      "North America–North America" = "#0072B2",  # blue
      "North America–South America" = "black",  # orange
      "South America–South America" = "#E69F00"   # green
    )
  ) +
  facet_wrap(~ ContinentComparison) +
  theme_minimal() +
  theme( panel.grid = element_blank())+
  labs(
    title = "Distance Decay by Continental Comparison",
    x = "Geographic Distance (km)",
    y = "Community Dissimilarity (Bray–Curtis)",
    color = "Comparison Type"
  )

#final plot
# 1. Duplicate the dataset for the "All Comparisons" panel
pairs_all <- pairs_df
pairs_all$ContinentComparison <- "All Comparisons"

# 2. Combine with original data
pairs_faceted <- rbind(pairs_df, pairs_all)

# 3. Plot
ggplot(pairs_faceted, aes(x = GeoDistance / 1000, y = CommunityDistance, color = ContinentComparison)) +
  geom_point(alpha = 0.4, size = 1.2) +
  geom_smooth(method = "lm", se = FALSE, size = 1, alpha = 0.8, color = "red") +
  scale_color_manual(
    values = c(
      "North America–North America" = "#0072B2",  # blue
      "North America–South America" = "#E69F00",  # orange
      "South America–South America" = "#009E73",  # green
      "All Comparisons" = "gray30"                # neutral dark gray
    )
  ) +
  facet_wrap(~ ContinentComparison, ncol = 2) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none",
    strip.text = element_text(face = "bold")
  ) +
  labs(
    title = "Distance Decay by Continental Comparison",
    x = "Geographic Distance (km)",
    y = "Community Dissimilarity (Bray–Curtis)"
  )



# Join continent info for both samples
pairs_df <- pairs_df %>%
  left_join(america_balanced %>% dplyr::select(SampleID, Continent.x), by = c("Sample1" = "SampleID")) %>%
  rename(Continent1 = Continent.x) %>%
  left_join(america_balanced %>% dplyr::select(SampleID, Continent.x), by = c("Sample2" = "SampleID")) %>%
  rename(Continent2 = Continent.x)

pairs_df_same_continent <- pairs_df %>%
  filter(Continent1 == Continent2)

ggplot(pairs_df, aes(x = GeoDistance, y = CommunityDistance, color = Continent1)) +
  geom_point(alpha = 0.3) +
  geom_smooth(method = "lm") +
  labs(
    title = "Distance Decay Relationship by Continent",
    x = "Geographic Distance (km)",
    y = "Community Dissimilarity (Bray-Curtis)",
    color = "Continent"
  ) +
  theme_minimal()


# --- Final Plots --- ####
library(ggplot2)
library(patchwork)


# * Graphs ####

# Use your custom theme everywhere
my_custom_theme <- theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black"),
    legend.key.size = unit(1, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# Define consistent climate zone colors (example; adjust levels & colors as needed)
climate_colors <- c(
  "Tropical" = "#1B9E77",    # teal-green
  "Arid" = "",        # bright orange
  "Temperate" = "blue",   # purple
  "Cold" = "#E7298A",        # pink/magenta
  "Polar" = "#66A61E"        # olive green
)


# Plot 1: Effect of Mean Temperature on Shannon Index
p1 <- ggplot() +
  geom_ribbon(data = pred_temp_gl, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "#CBB3E6", alpha = 0.25) +
  geom_point(data = global_balanced, aes(x = mean_temp, y = Shannon, color = Koppen_Broad),
             size = 2, alpha = 0.7, shape = 16) +
  geom_line(data = pred_temp_gl, aes(x = x, y = predicted), color = "black", size = 1.3) +
  labs(
    title = "Effect of Mean Temperature on Shannon Index",
    subtitle = "p = 1.14e-06 *** \nConditional R2: 0.952\nMarginal R2: 0.036",
    x = "Mean Temperature",
    y = "Predicted Shannon Index",
    color = "Climate Zone"
  ) +
  my_custom_theme +
  scale_color_manual(values = climate_colors)

# Plot 2: Interaction Climate Zone × Temperature
p2 <- ggplot(interaction_preds, aes(x = x, y = predicted, color = group)) +
  
  geom_point(data = global_balanced, aes(x = mean_temp, y = Shannon, color = Koppen_Broad),
             size = 2, alpha = 0.3, shape = 16) +
  geom_line(size = 1) +
  labs(
    title = "Interaction: Climate Zone × Temperature",
    subtitle = "Climate Zone p= 0.003713 ** \nTemperature p= 6.824e-10 ***\nClimate Zone:Temperature p= 0.003719**",
    x = "Mean Temperature",
    y = "Predicted Shannon Index",
    color = "Climate Zone"
  ) +
  my_custom_theme +
  scale_color_manual(values = climate_colors)

# Plot 3: Effect of Soil pH on Shannon Index )
p3 <- ggplot() +
  geom_point(data = america_no_outliers, aes(x = pH_from_raster, y = Shannon_sqrt, color = Continent.x),
             size = 2, alpha = 0.7, shape = 16) +
  geom_ribbon(data = pred_ph, aes(x = x, ymin = conf.low, ymax = conf.high), fill = "#CBB3E6", alpha = 0.25) +
  geom_line(data = pred_ph, aes(x = x, y = predicted), color = "black", size = 1.3) +
  scale_color_manual(values = continent_colors) +
  labs(
    title = "Effect of Soil pH on Shannon Index (sqrt-transformed)",
    subtitle = "p = 0.01931 * \nConditional R2: 0.913\nMarginal R2: 0.029",
    x = "Soil pH",
    y = "Predicted Shannon Index (sqrt)",
    color = "Continent"
  ) +
  my_custom_theme

# Plot 4: Shannon Index predicted by Continent 
p4 <- ggplot(america_no_outliers, aes(x = Continent.x, y = Shannon_sqrt, fill = Continent.x)) +
  geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.7,
               position = position_dodge(width = 0.8)) +
  geom_jitter(aes(color = Continent.x),
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.4, size = 1.5) +
  scale_fill_manual(values = continent_colors) +
  scale_color_manual(values = continent_colors) +
  labs(
    title = "Shannon-Index predicted by Continent",
    subtitle = "p = 0.00386 **\nEstimate = -0.10830\nConditional R2: 0.913\nMarginal R2: 0.029",
    x = "",
    y = "Predicted Shannon Index (sqrt-transformed)"
  ) +
  my_custom_theme +
  theme(legend.position = "none")  # no legend for boxplot


#shannon by continent boxplot
p5 <- ggplot(data_plot, aes(x = Continent.x, y = Shannon, fill = Continent.x)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Shannon Diversity by Continent",
    y = "Shannon Index",
    x = "",
    fill = "Continent"  # ✅ correct way to set legend title
  ) +
  my_custom_theme +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_brewer(palette = "Set1")


#shannon by continent boxplot
p6 <- ggplot(data_plot, aes(x = Koppen_Broad, y = Shannon, fill = Koppen_Broad)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  theme_minimal() +
  labs(
    title = "Shannon Diversity by Climate Zone",
    y = "Shannon Index",
    x = "",
    fill = "Koppen_Broad"  
  ) +
  my_custom_theme +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  scale_fill_brewer(palette = "Set1")


# Combine all plots into a 2x2 grid
combined_plot <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Summary of Ecological Model Results",
    theme = theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))
  )
# Print combined plot
print(combined_plot)

# Save combined plot
ggsave("Graphs/combined_plot.png", combined_plot, width = 16, height = 12, dpi = 300)

# Define a list of plots and filenames
plot_list <- list(
  p1 = p1,
  p2 = p2,
  p3 = p3,
  p4 = p4,
  p5 = p5,
  p6 = p6
)

# Save each plot
for (plot_name in names(plot_list)) {
  ggsave(
    filename = file.path("Graphs", paste0(plot_name, ".png")),
    plot = plot_list[[plot_name]],
    width = 8, height = 6, dpi = 300
  )
}

# * Maps ####

#theme
map_theme <- theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )


m1 <- ggplot() +
  geom_sf(data = world, fill = "gray90", color = "gray70") +
  geom_sf(data = df.sf, aes(color = Authors), alpha = 0.7, size = 2) +
  theme_minimal() +
  labs(
    title = "Distribution of Soil Microbes",
    subtitle = "Primer V4, Soil",
    color = "Author"  # optional legend label
  ) +
  map_theme

#koppen zones map
m2 <- ggplot() +
  geom_tile(data = koppen_df, aes(x = x, y = y, fill = zone_broad)) +
  scale_fill_viridis_d(
    name = "Köppen Climate Zone",
    option = "viridis",
    direction = -1  # reverses the color order
  ) +
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Global Köppen-Geiger Climate Zones",
    subtitle = "Based on 1991–2020 data"
  ) +
  map_theme+
  theme(
    legend.position = "right",         # place legend next to the map
    legend.key.size = unit(0.8, "lines"), # smaller legend keys
    legend.text = element_text(size = 8), # smaller legend text
    legend.title = element_text(size = 10), # smaller legend title
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )


#Add Shannon Index data points
m3 <- ggplot() +
  # Köppen climate background
  geom_tile(data = koppen_df, aes(x = x, y = y, fill = zone_broad)) +
  scale_fill_manual(
    name = "Köppen Climate Zone",
    values = c(
      "Polar" = "#cccccc",  
      "Cold" = "#999999",      
      "Temperate" = "#777777", 
      "Arid" = "#555555",      
      "Tropical" = "#444444"      
    )
  ) +
  
  # Shannon diversity points
  # Fill layer
  geom_point(
    data = data_plot,
    aes(x = Longitude, y = Latitude, size = 2),
    color = "#4DAF4A",
    fill = "#4DAF4A", # required for shape 21
    alpha = 0.7,
    shape = 21
  )+
  coord_equal() +
  theme_minimal() +
  labs(
    title = "Köppen Climate Zones with Shannon Diversity",
    subtitle = "",
    x = "",
    y = ""
  ) +
  map_theme


#pH Global Map
m4 <-ggplot() +
  geom_tile(data = ph_df, aes(x = x, y = y, fill = pH)) +  # raster layer
  geom_sf(data = world, fill = NA, color = "gray70", size = 0.3) + # borders
  # geom_sf(data = df.sf, aes(color = "#4DAF4A"), size = 2) +  # points colored by extracted pH
  scale_fill_viridis_c(name = "Soil pH (Raster)", option = "plasma", direction = -1) +
  #scale_color_viridis_c(name = "Soil pH (Points)", option = "plasma", direction = -1) +
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Global Soil pH (0–30 cm)",
    subtitle = "Source: SoilGrids",
    x = "Longitude",
    y = "Latitude"
  ) +
  map_theme

#pH Global Map with datapoints
m6 <- ggplot() +
  geom_tile(data = ph_df, aes(x = x, y = y, fill = pH)) +  # raster layer
  geom_sf(data = world, fill = NA, color = "gray70", size = 0.3) + # borders
 geom_sf(data = df.sf, aes(), size = 2, color = "#4DAF4A") +  # points colored by extracted pH
  scale_fill_viridis_c(name = "Soil pH (Raster)", option = "plasma", direction = -1) +
  #scale_color_viridis_c(name = "Soil pH (Points)", option = "plasma", direction = -1) +
  coord_sf() +
  theme_minimal() +
  labs(
    title = "Global Soil pH (0–30 cm)",
    subtitle = "Source: SoilGrids",
    x = "Longitude",
    y = "Latitude"
  ) +
  map_theme


#pH america all data
m9 <- ggplot() +
  geom_tile(data = ph_df, aes(x = x, y = y, fill = pH)) +  # raster layer
  geom_sf(data = world, fill = NA, color = "gray70", size = 0.3) + # borders
  geom_sf(data = , aes(), color = "#4DAF4A", size = 2, alpha = 0.7) +
  #geom_hline(yintercept = 9, color = "red", linetype = "dashed", size = 1) +  # line at 9° latitude
  scale_fill_viridis_c(name = "Soil pH (Raster)", option = "plasma", direction = -1) +
  coord_sf(xlim = c(-170, -30), ylim = c(-60, 75), expand = FALSE) +  # zoom to Americas
  theme_minimal() +
  labs(
    title = "America Soil pH (0–30 cm)",
    subtitle = "Source: SoilGrids",
    x = "Longitude",
    y = "Latitude"
  ) +
  map_theme



#pH america cluster data
m7 <- ggplot() +
  geom_tile(data = ph_df, aes(x = x, y = y, fill = pH)) +  # raster layer
  geom_sf(data = world, fill = NA, color = "gray70", size = 0.3) + # borders
  geom_sf(data = cluster_sf, aes(), color = "green", size = 4, alpha = 0.7) +
  #geom_hline(yintercept = 9, color = "red", linetype = "dashed", size = 1) +  # line at 9° latitude
  scale_fill_viridis_c(name = "Soil pH (Raster)", option = "plasma", direction = -1) +
  coord_sf(xlim = c(-170, -30), ylim = c(-60, 75), expand = FALSE) +  # zoom to Americas
  theme_minimal() +
  labs(
    title = "America Soil pH (0–30 cm)",
    subtitle = "Source: SoilGrids",
    x = "Longitude",
    y = "Latitude"
  ) +
  map_theme


# Plot global soil temp
m8 <- ggplot() +
  geom_tile(data = temp_df, aes(x = x, y = y, fill = mean_temp)) +  # raster layer
  geom_sf(data = world, fill = NA, color = "gray70", size = 0.3) +  # borders
  geom_sf(data = df.sf, color = "#4DAF4A", size = 2) +              # fixed green color_
  
  scale_fill_viridis_c(
    name = "Soil Temperature (Raster)",
    option = "plasma",
    direction = 1
  ) +
  scale_color_viridis_c(
    name = "Soil Temperature (Points)",
    option = "plasma",
    direction = 1
  ) +
  coord_sf() +   
  # geom_point(data = data_plot, aes(x = Longitude, y = Latitude, size = Shannon),color = "#4DAF4A", alpha = 0.7 ) +
  theme_minimal() +
  labs(
    title = "Global Soil Temperature (0–5 cm)",
    subtitle = "Source: ",
    x = "Longitude",
    y = "Latitude"
  )+
  map_theme

# Define a named list of map figures
figures <- list(
  m1 = m1,
  m2 = m2,
  m3 = m3,
  m4 = m4,
  m6 = m6,
  m7 = m7,
  m8 = m8,
  m9 = m9
)

# Create output directory if it doesn't exist
output_dir <- "figures"
dir.create(output_dir, showWarnings = FALSE)

# Save each figure using consistent formatting
for (name in names(figures)) {
  ggsave(
    filename = file.path(output_dir, paste0(name, ".png")),
    plot = figures[[name]],
    width = 8,
    height = 6,
    dpi = 300,
    bg = "white"
  )
}



# 10. Cross Validation ####

# Global ####

run_model_once <- function(data_plot) {
  
  # Step 1: Filter out sparse continents
  data_continent_filtered <- data_plot %>%
    filter(!Continent.x %in% c("Antarctica", "Oceania", "Europe", "Africa", "Central America"),
           !is.na(Continent.x))
  
  # Step 2: Downsample continents
  min_n_continent <- data_continent_filtered %>%
    count(Continent.x) %>%
    summarise(min_n = min(n)) %>%
    pull(min_n)
  
  continent_balanced <- data_continent_filtered %>%
    group_by(Continent.x) %>%
    slice_sample(n = min_n_continent) %>%
    ungroup()
  
  # Step 3: Filter out sparse climate zones
  data_koppen_filtered <- continent_balanced %>%
    filter(!Koppen_Broad %in% c("Arid", "Polar"),
           !is.na(Koppen_Broad))
  
  # Step 4: Downsample climate zones
  min_n_koppen <- data_koppen_filtered %>%
    count(Koppen_Broad) %>%
    summarise(min_n = min(n)) %>%
    pull(min_n)
  
  global_balanced <- data_koppen_filtered %>%
    group_by(Koppen_Broad) %>%
    slice_sample(n = min_n_koppen) %>%
    ungroup()
  
  # Step 5: Fit the model
  model <- lmer(Shannon ~ Koppen_Broad * mean_temp + 
                  (1 | City_Region.x) + (1 | Authors.x),
                data = global_balanced)
  
  # Step 6: Extract fixed effects and p-values
  tidy(model, effects = "fixed") %>%
    dplyr::select(term, estimate, p.value)
}


# --- Run simulation over multiple random samplings ---
set.seed(42)
n_iter <- 100

model_results <- map_dfr(1:n_iter, ~run_model_once(data_plot), .id = "replicate")

# --- Summarize robustness of model results ---
summary_stats <- model_results %>%
  group_by(term) %>%
  summarise(
    mean_estimate = mean(estimate),
    sd_estimate = sd(estimate),
    prop_significant = mean(p.value < 0.05),
    .groups = "drop"
  )

# --- View results ---
print(summary_stats)

library(knitr)
library(kableExtra)


summary_stats %>%
  kable(
    caption = "Summary of Model Estimates and Significance Proportions",
    digits = 4,
    col.names = c("Term", "Mean Estimate", "SD of Estimate", "Proportion Significant")
  ) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed"))




# Americas ####

run_model_once <- function(data_plot) {
  
  #get America data only
  america_filtered <- data_plot %>%
    filter(
      Continent.x %in% c("North America", "Central America", "South America"),
      Latitude >= -60, Latitude <= 85,          # roughly Americas latitude range
      Longitude >= -170, Longitude <= -30       # roughly Americas longitude range
    )
  
  #Set border between north and south america
  america_panama <- america_filtered %>%
    mutate(Continent.x = if_else(
      Latitude >= 9 | (Latitude < 9 & Longitude < -79.9),
      "North America",
      "South America"
    ))
  
  
  # Find the minimum group size
  min_n <- america_panama %>%
    count(Continent.x) %>%
    summarise(min_n = min(n)) %>%
    pull(min_n)
  
  # Downsample each group to the same size
  america_balanced <- america_panama %>%
    group_by(Continent.x) %>%
    slice_sample(n = min_n) %>%
    ungroup()
  
  america_balanced %>%
    count(Continent.x)
  
  america_balanced$Shannon_sqrt <- sqrt(america_balanced$Shannon)
 
  ph_america <- lmer(Shannon_sqrt ~ pH_from_raster + mean_temp + Koppen_Broad + Continent.x + (1|City_Region.x) + (1|Authors.x), data = america_balanced)
  
 
  # Fit the model
  model <- lmer(Shannon_sqrt ~ pH_from_raster + mean_temp + Koppen_Broad + Continent.x + (1|City_Region.x) + (1|Authors.x), data = america_balanced)
  
  # Step 6: Extract fixed effects and p-values
  tidy(model, effects = "fixed") %>%
    dplyr::select(term, estimate, p.value)
}


# --- Run simulation over multiple random samplings ---
set.seed(42)
n_iter <- 100

model_results <- map_dfr(1:n_iter, ~run_model_once(data_plot), .id = "replicate")

# --- Summarize robustness of model results ---
summary_stats <- model_results %>%
  group_by(term) %>%
  summarise(
    mean_estimate = mean(estimate),
    sd_estimate = sd(estimate),
    prop_significant = mean(p.value < 0.05),
    .groups = "drop"
  )

# --- View results ---
print(summary_stats)

library(knitr)
library(kableExtra)


summary_stats %>%
  kable(
    caption = "Summary of Model Estimates and Significance Proportions",
    digits = 4,
    col.names = c("Term", "Mean Estimate", "SD of Estimate", "Proportion Significant")
  ) %>%
  kable_styling(full_width = FALSE, bootstrap_options = c("striped", "hover", "condensed"))







