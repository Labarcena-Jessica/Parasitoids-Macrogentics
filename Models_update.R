#Retrieveing climatic data from WorldClim using the cell centre coordinates for the modelling of the Nucleotide Diversity values. In this case I download Bioclimatic variables. 

library(raster)
library(sp)

clim_data <- getData("worldclim",var="bio",res=5)

#selecting only 5 relevant variables
clim_data <- r[[c(1,4, 7, 12, 15)]]
names(clim_data) <- c("Mean_Temp", "Temp_Seas", "Temp_Range", "Prec", "Prec_Seas")

#read data with coordinates corrected for WorldClim. In this case some of the cell centre coordinates where in the ocean as the cells were in the edge of land. I manually checked for this cells and decide to correct the coordinates using one of the cell coordinates (one of the hexagon points) in land. 

#Joining the data used for the models which includes amount of data per cell, coordinates, mean nucleotide diversity and estimated "especies richness". 
cell_correlation <- left_join(cell_centre_df_more10, data_GD10_table, by = "cell")
cell_correlation <- left_join (cell_correlation, mean_by_cell_final10, by = "cell")

write_csv(cell_correlation, "cell_correlation_parasitoids.csv")

#Getting the coordinates for extracting climatic values 

lats <- cell_correlation$lat
lons <- cell_correlation$long 

coords <- data.frame(x=lons,y=lats)
#transforming to spatial points
points <- SpatialPoints(coords, proj4string = r@crs)
#extracting values (this is done with raster, will need to check how to move to terra)
values <- extract(r,points)

df_climatic <- cbind.data.frame(coordinates(points),values)
colnames(df_climatic) <- c("long", "lat", "Mean_Temp", "Temp_Seas", "Temp_Range", "Prec", "Prec_Seas")

#joining all data
data_parasitoids_models <- left_join(cell_correlation, df_climatic, by = c("lat", "long"))


#adding "species richness" estimates using the point estimate calculated at 70% of coverage. 
data_parasitoids_models <- left_join(data_parasitoids_models, point_estimate, by = "cell")
names(data_parasitoids_models)
colnames(data_parasitoids_models) <- c("cell", "mean_nuc.div", "sequences","bins", "long", "lat", "Mean_Temp", "Temp_Seas","Temp_Range", "Prec",              "Prec_Seas",  "Coverage","Species Richness" ,"Lower Confidence Level" , "Upper Confidence Level" , "Climatic_Zone" )

write_csv(data_parasitoids_models, "data_parasitoids_models.csv" )

data_parasitoids_models <- read.csv("data_parasitoids_models_JULY2023.csv")

#Remove the cells that were considered as outliers because those have the higher amount of data an influenced the species richness estimation. 
data_parasitoids_models <- filter(!is.na(data_parasitoids_models$`Species Richness`))

#Linear model relationships with variables. Temperatures are in C *10 in WorldClim

data_parasitoids_models$Mean_Temp <- data_parasitoids_models$Mean_Temp/10

plot((data_parasitoids_models$Mean_Temp), sqrt(data_parasitoids_models$mean_by_cell_pi10), xlab = "Mean Temperature (°C)", ylab = "Mean Nucleotide Diversity per Cell", pch = 19, col = "dark grey")
abline(lm(sqrt(data_parasitoids_models$mean_by_cell_pi10) ~ data_parasitoids_models$Mean_Temp), col = "red", lwd = 3)


plot((abs(data_parasitoids_models$lat)), sqrt(data_parasitoids_models$mean_by_cell_pi10), xlab = "Absolute Latitude", ylab = "Mean Nucleotide Diversity per Cell", pch = 19, col = "dark grey",  cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
abline(lm( sqrt(data_parasitoids_models$mean_by_cell_pi10) ~ abs(data_parasitoids_models$lat)), col = "red", lwd = 3)

lm_L_Nuc2 <- lm (sqrt(data_parasitoids_models_F$mean_by_cell_pi10) ~ abs(data_parasitoids_models_F$lat))


summary(lm_L_Nuc2)
plot(lm_L_Nuc2)

#Removing outliers 

data_parasitoids_models_F <- data_parasitoids_models[c(-15, -37), ]

plot(log(data_parasitoids_models$bins), sqrt(data_parasitoids_models$mean_by_cell_pi10), xlab = "Log(Number of Molecular Clusters)", ylab = "Mean Nucleotide Diversity per Cell", pch = 19, col = "dark grey",  cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5 )
abline(lm(sqrt(data_parasitoids_models$mean_by_cell_pi10) ~ log(data_parasitoids_models$bins)), col = "red", lwd = 3)

lm_bins_nucdiv <- lm(sqrt(data_parasitoids_models_F$mean_by_cell_pi10) ~ log(data_parasitoids_models_F$bins))
summary(lm_bins_nucdiv)

plot(lm_bins_nucdiv)

# get list of residuals 
res <- resid(lm_L_Nuc2)
res

# produce residual vs. fitted plot
plot(fitted(lm_L_Nuc2), res , xlab = "Fitted values", ylab = "Residuals", pch = 19, col = "dark grey")

# add a horizontal line at 0 
abline(0,0)

# create Q-Q plot for residuals
qqnorm(res)

# add a straight diagonal line 
# to the plot
qqline(res)

#Trying to generate figure 1 by combining all the figures 
install.packages("cowplot")
library(cowplot)
install.packages("gridGraphics")
library(gridGraphics)

all_data_combined <- read.csv("all_data_combined.csv")

# Initiate a ggplot
b <- ggplot(all_data_combined, aes(x = lat.y, y = mean_by_cell_pi10))+
  geom_point(aes(color = family))+
  geom_smooth(aes(color = family))+
  scale_colour_viridis(discrete = TRUE)+
  labs(colour='Family')+
  labs(y="Mean Nucleotide Diversity per Cell", x = "Latitude") +
  theme_bw() 

b + theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 20))  +
  theme(legend.title = element_text(size = 14)) +
  theme(legend.text = element_text(size = 12))




plot_grid(mapmore10.3_parasitoids, map_cell_centre, b, labels = 'AUTO')


#Adding unique number of locations per BIN to the model analysis

data_perBIN10_all <- left_join(data_perBIN10_all, distance_bin_df_braco, by ="bin_uri")

lm3 <- lm(data_perBIN10_all$Pi_per_site ~ data_perBIN10_all$unique_locations.y)
summary(lm3)

lm4 <- lm(sqrt(data_perBIN10_all$Pi_per_site) ~ log(data_perBIN10_all$unique_locations.y))
summary(lm4)

plot(log(data_perBIN10_all$unique_locations.y), sqrt(data_perBIN10_all$Pi_per_site))
abline(lm(sqrt(data_perBIN10_all$Pi_per_site) ~ log(data_perBIN10_all$unique_locations.y)), col = "red", lwd = 3)
plot(lm4)


