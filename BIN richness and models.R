#BIN richness calculation and accumulation curves. 

#The parasitoids data frame has all the information for Ichneumonidae + Braconidae and will be used to calculate species richness per cell for both groups. This include all data filtered for presence of BINs, lat/lon, sequence quality and outlier detection. 

parasitoids_SR <- data_geo 

#creating a grid c ell with same spatial resolution and adding the cell number ot the data frame. 
gridSR <- dgconstruct(res = 5)

parasitoids_SR$cell <- dgGEO_to_SEQNUM(gridSR, parasitoids_SR$lon, parasitoids_SR$lat)$seqnum

#count records per BINs per cell 
total_record_parasitoids <- parasitoids_SR %>%
  dplyr::count(cell,bin_uri)

bin_cell_SR <- total_record_parasitoids %>%
  group_by(cell) %>%
  dplyr::count(length(bin_uri))

bin_cell_SR <- bin_cell_SR [, c(1,3)]


#Adding count(number of BINs per Cell to the parasitoids_SR data set)
parasitoids_SR <- left_join(parasitoids_SR, bin_cell_SR, by =c("cell"))

#Cells with less than 3 BINs were filtered to keep same filtering as with the Nuc div analysis. In this case we don't filter by the number of sequences per BIN as singletons and doubletons are relevant in the abundance based estimate. 

parasitoids_SR <- filter(parasitoids_SR, n >= 3)

#creating community object to use with the iNEXT package. This represents the number of sequences per BIn in each cell. 

parasitoids_BIN_percell <- parasitoids_SR %>%
  group_by(cell,bin_uri) %>%
  dplyr::summarize(count = length((bin_uri)))

data_perBINcomm <- pivot_wider (parasitoids_BIN_percell, names_from = bin_uri , values_from = count)

data_perBINcomm <- data_perBINcomm %>%
  remove_rownames() %>%
  column_to_rownames(var = 'cell')

data_perBINcomm[is.na(data_perBINcomm)] <- 0


install_github('JohnsonHsieh/iNEXT')
library(iNEXT)

parasitoids_BIN_percell <- data.frame(parasitoids_BIN_percell)

data_abundance <- split(parasitoids_BIN_percell$count, parasitoids_BIN_percell$cell)

#We decided to remove the outliers cells with the highest amount of data. (cells 1297, 199). This is because the estimation of species richness having these two data rich cells was affected. 

cells_models2 <- c("10", "1000", "1009", "1018" ,"104" , "105",  "1085", "113" , "114"  , "1222"  , "126" , "1279",  "131" , "132" , "141" , "142" , "1481" ,"150" , "151" , "154" , "163" , "170",  "1708", "171" , "1755" ,"181" , "186",  "189" , "19"  , "190" , "195",  "1951", "1977" , "204",  "207",  "208" , "2092", "2110", "2111", "2122", "2129", "2131", "2133", "2164", "217",  "226" , "2342", "2351", "248" , "256", "257",  "36" ,  "39" ,  "40" ,  "496" , "514",  "523" , "532",  "559", "59",   "618",  "636",  "67" ,  "712" , "717",  "727" , "757",  "76",   "77",   "78",   "782",  "8",    "86",   "87" ,  "881" , "94",   "954" , "970",  "982",  "991",  "992", "124", "122")

#selecting just the cells included. 

data_abundance2 <- data_abundance [cells_models2]

# set a series of sample sizes (m) for R/E computation
m <- c(10, 50, 100, 150, 200, 250, 300, 400, 450, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500)

#richness estimation usign Hill Numbers

rich_estimates_parasitoids <- iNEXT(data_abundance2, q=0, datatype="abundance")
head(rich_estimates_parasitoids[2])

data_iNEXT <- rich_estimates_parasitoids$iNextEst$coverage_based

data_iNEXTf <- data_iNEXT %>%
  filter(Method == "Observed")

df_SR <- fortify(rich_estimates_parasitoids, type = 1)
head(df_SR)

names(df_SR)
names(df_SR) <- c("datatype", "plottype", "cell", "Method", "Order.q", "x", "y", "y.lwr", "y.upr")

#Adding latitudinal information 

df_SR2 <- left_join(df_SR, cell_lat_info, by = "cell")

#Plotting rearefaction curves 
df_SR2.point <- df_SR2[which(df_SR2$Method=="Observed"),]
df_SR2.line <- df_SR2[which(df_SR2$Method!="Observed"),]
df_SR2.line$Method <- factor(df_SR2.line$Method, c("Rarefaction", "Extrapolation"), c("interpolation", "extrapolation"))


#Climatic Band, trying adding the Climatic zones, as this is consistent with the previous analysis
data_CZ_SR <- df_SR2 

data_CZ_SR$rndCoord.lat <- RoundCoordinates((data_CZ_SR$lat), res = "fine", latlong = "lat")

data_CZ_SR$rndCoord.lon <- RoundCoordinates((data_CZ_SR$lon), res = "fine", latlong = "lon")

data_CZ_SR$ClimateZ <- LookupCZ(data_CZ_SR, res = "fine")


write.csv(data_CZ_SR, "data_CZ_SR_update.csv")

#Data after correcting the coordinates of the ones in the Ocean. This is because the rounding issue with the kgc package. I used a manual correction to include the coordinates in land and search the climatic zone using the same code. Then I added the Climatic_BAD column. 

data_CZ_SR2 <- read.csv("data_CZ_SR_update.csv")

CZ_SR <- data_CZ_SR2 [, c(3,12, 14)]

CZ_SR$cell <- as.character(CZ_SR$cell)

data_CZ_SR2 <- left_join (df_SR2 , CZ_SR, by = "cell")

data_CZ_SR2 <- data_CZ_SR2[ , -1]

#Rarefaction curves colour coding for the Climatic Zone 
data_CZ_SR2$Climatic_BAND <- factor(data_CZ_SR2$Climatic_BAND , levels = c("Tropical ", "Arid/Dry", "Warm/Mild Temperate", "Continental", "Polar/Tundra"))

inext_graph <-  ggplot(data_CZ_SR2, aes(x=x, y=y, group = cell)) +
  geom_line()+
  #geom_point(aes(shape=cell), size=5, data=df2.point) +
  #geom_line(aes(linetype = Method), lwd=1.5, data=df2.line) +
  geom_ribbon(aes(ymin=y.lwr, ymax=y.upr,
                 fill=Climatic_BAND), alpha=0.6) +
  labs(x="Number of records", y="Molecular clusters relative diversity") +
  theme_bw()+
  scale_fill_viridis(direction = 1, discrete = TRUE, alpha=0.6)+
  #scale_fill_viridis(discrete = TRUE)+
  labs(fill='Climatic Zones')
  #theme(legend.position = "none") 

inext_graph
inext_graph +
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 20))  +
  theme(legend.title = element_text(size = 14))


# The lapply function applies the same logic across elements in a list. Selecting value of species richess at the same coverage value of 1000
point_1000 <- lapply(rich_estimates_parasitoids$iNextEst, function(x) { x[ x$m == 1000, ] })
# Turn the output into a dataframe
point_df <- point_1000 [[1]]
# Extract the strata
point_df$Strata <- gsub("\\..*","",rownames(point_df))

q.vals <- c("test", "why")
point_df$Assemblage <- as.factor(point_df$Assemblage)
# Make a nice ggplot!

names(point_df)
names(point_df) <- c("cell","m","Method","Order.q", "qD", "qD.LCL", "qD.UCL", "SC", "SC.LCL", "SC.UCL","Strata")

df_SR3 <- data_CZ_SR2m[, c(5, 12, 13, 14, 20)] #cell, mean_bycellPi, Climatic Zone, lat long
df_SR3 <- distinct(df_SR3)

df_SR3$cell <- as.character(df_SR3$cell)

point_df <- left_join(point_df, df_SR3, by = "cell")

#another way of estimating diversity for a given coverage. level = 0.70 means 70% of coverage

point_estimate <- estimateD(data_abundance2, q = 0, datatype="abundance", base="coverage", level=0.70, conf=0.95)

names(point_estimate) <- c("cell","m","Method","Order.q", "SC", "qD", "qD.LCL", "qD.UCL")

point_estimate <- left_join(point_estimate, CZ_SR, by = "cell")
point_estimate <- distinct(point_estimate)

point_estimate$Climatic_BAND <- factor(point_estimate$Climatic_BAND , levels = c("Tropical ", "Arid/Dry", "Warm/Mild Temperate", "Continental", "Polar/Tundra"))

write.csv (point_estimate, "point_estimate.csv")

#Graphic of the "species richness" per cell color coded for the climatic zone. 

SR_plot <- ggplot(point_estimate, aes(x=cell, y=qD, colour=Climatic_BAND, shape = Method)) + 
  theme_bw() +
  #scale_x_discrete(breaks=c("1","2"),labels= c("1","2")) +
  geom_errorbar(aes(ymin=qD.LCL, ymax=qD.UCL), width=.01) +
  labs(y="Estimated richness", x = "Cell") +
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 45, size = 5))+
  scale_colour_viridis(direction = 1, discrete = TRUE)+
  #scale_colour_viridis(discrete = TRUE)+
  labs(colour='Climatic Zones', shape = 'Method')

SR_plot

SR_plot + theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 20))  +
  theme(legend.title = element_text(size = 14))

#Linear model of the mean nucleotide diversity values per cell and the "species richness" estimate.

plot(sqrt(point_estimate$qD), sqrt(point_estimate$mean_by_cell_pi10), xlab = "Estimated richness per Cell", ylab = "Mean Nucleotide Diversity per Cell", pch = 19, col = "dark grey", cex.lab=1.5, cex.axis=1.5, cex.main=1.5, cex.sub=1.5)
abline(lm(sqrt(point_df$mean_by_cell_pi10) ~ sqrt(point_df$qD)), col = "red", lwd = 3)


ggplot(point_estimate, aes(x =sqrt(point_estimate$qD), y = sqrt(point_estimate$mean_by_cell_pi10))) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  theme_bw()


#Saving values of Species richness estimates using iNEXT.
sps_riche_parasitoids <- rbind(rich_estimates_parasitoids$AsyEst)
sps_riche_parasitoids <- filter(sps_riche_parasitoids, Diversity == 'Species richness')
names(sps_riche_parasitoids)
colnames(sps_riche_parasitoids) <- c("cell", "Diversity", "Observed",  "Estimator", "s.e.", "LCL", "UCL")

write_csv(sps_riche_parasitoids , "Species_div_parasitoids,csv")

#To explore the potential determinants of intraspecific genetic diversity, we used linear models to test the influence of geographic and environmental variables. We modelled the nucleotide diversity values per cell as a function of the latitude and climatic variables. Using the grid cell centre coordinates, we retrieved bioclimatic variables from WorldClim. Particularly Annual Mean Temperature, Temperature Seasonality, Temperature Annual Range, Annual Precipitation and Precipitation Seasonality. These variables represent annual trends and were used as predictors of the nucleotide diversity values.Regression analyses were conducted using geographic and environmental factors as predictor variables. 

## Bring in predictor variables

a<-abs(data_parasitoids_models$lat)
b<-(data_parasitoids_models$Mean_Temp/10)
c<-(data_parasitoids_models$Temp_Seas/10)
d<-(data_parasitoids_models$Temp_Range)
e<-(data_parasitoids_models$Prec)
f<-(data_parasitoids_models$Prec_Seas)
g<-(data_parasitoids_models$Estimator)


variablenames <- data.frame(names(data_parasitoids_models[,6:12]))
colnames(variablenames) <- "varname"

## response

y<-sqrt(data_parasitoids_models$mean_by_cell_pi10) # square root of the mean nuc div per cell 

## Calculate Variance Inflation Factor to decide which predictors to drop
## remove predictors stepwise
library(car)
full <- lm(y ~ a + b + c + d + e + f + g )
vif(full)

## First for response variable y
## remove c, which has the highest VIF
reduced1 <- lm(y ~ a+ b + d + e + f + g)
vif(reduced1)
## All models now have VIF <= 10 

## Full model used for the response variable
full.nd <- y ~ a + b + d + e + f + g 

#Another approach to model averaging: use MuMIn. I like this because it is simple and highly repeatable (i.e., does not rely on many custom functions)

library(MuMIn)

options(na.action = "na.fail")

## Generate full set of possible models
nd1 <- dredge(lm(full.nd))

## Rel importance
#Warning message:'importance' is deprecated. Use 'sw' instead.

import.nd <- sw(model.avg(nd1))

## Get best model and coefficients

summary(get.models(nd1, 1)[[1]])

nd1_f <- subset(nd1, delta < 3)

par(mar = c(3,5,6,4))

plot(nd1_f, labAsExpr = TRUE)



#PLOT the best model 

plot((data_parasitoids_models$Temp_Range/10), sqrt(data_parasitoids_models$mean_by_cell_pi10), xlab = "Temperature Range (°C)", ylab = "Mean Nucleotide Diversity per Cell", pch = 19, col = "dark grey")
abline(lm(sqrt(data_parasitoids_models$mean_by_cell_pi10) ~  y+1), col = "red", lwd = 3)

y <- data_parasitoids_models$Temp_Range/10
lm_TR_nucdiv <- lm(sqrt(data_parasitoids_models$mean_by_cell_pi10) ~ (data_parasitoids_models$Temp_Range) +1)

summary(lm_TR_nucdiv)

ggplot(data_parasitoids_models, aes(x = Temp_Range/10, y = sqrt(data_parasitoids_models$mean_by_cell_pi10))) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red") +
  theme_bw()
  
                   
## Plot output
#########################3
names(data_parasitoids_models)
namesforplot <- c("A.Latitude", "Mean T." ,"T. Range", "Prec" ,"Prec S." ,"E. Richness" )

pdf("rel_import_nd1.pdf", width=3, height=6)
par(mfrow=c(3,1), mar=c(5,5,2,1))

barplot(rev(as.vector(import.nd)), horiz=T, names=rev(namesforplot[as.numeric(as.factor(names(import.nd)))]), las=2, xlim=c(0,1.2), cex.names=0.75, axes=F, xlab="Relative importance (AICc)", main="Mean Nucleotide Diversity", col="grey")
axis(1, at=seq(0,1,0.2))
text(rep(1.05, 12),(12:1*1.2)-0.6, c("+","-","+","-","-","+"), cex=1.25)
text(1.1,(12*1.2)-0.6,"*", cex=1.25, col="red")

dev.off()



residuals_nd1 <- get.models(nd1, 1)[[1]]

summary(get.models(nd1, 2)[[1]])

residuals_nd1.2 <- get.models(nd1, 2)[[1]]

summary(get.models(nd1, 2)[[1]])

residuals_nd1.3 <- get.models(nd1, 3)[[1]]

residuals_nd1.4 <- get.models(nd1, 4)[[1]]

summary(get.models(nd1, 4)[[1]])


residuals_nd1.5 <- get.models(nd1, 5)[[1]]

summary(get.models(nd1, 5)[[1]])

w <- (data_parasitoids_models$Temp_Range/10 + abs(data_parasitoids_models$lat) +1)

plot((data_parasitoids_models$Temp_Range/10 + abs(data_parasitoids_models$lat) +1), sqrt(data_parasitoids_models$mean_by_cell_pi10), xlab = "Temperature Range (°C)", ylab = "Mean Nucleotide Diversity per Cell", pch = 19, col = "dark grey")
abline(lm(sqrt(data_parasitoids_models$mean_by_cell_pi10) ~  w), col = "red", lwd = 3)






#Modified t.test for spatial autocorrelated data in the relationship of nuc div and sps richness.

cell_correlation2 <- left_join(cell_correlation2, sps_riche_parasitoids, by = "cell")

cell_correlation2 <- cell_correlation2 [, c(-7, -8, -10, -11, -12)]

cell_correlation_fil <- filter(cell_correlation2, bins < 105) 

coords_t_lm <- as.matrix(cell_correlation2[, c(5,6)])

correlation_nucdiv_spsr <- modified.ttest(cell_correlation2$Estimator, sqrt(cell_correlation2$mean_by_cell_pi10), coords_t_lm)

print(correlation_nucdiv_spsr)


