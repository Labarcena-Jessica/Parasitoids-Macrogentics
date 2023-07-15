install.packages("kgc")
library(kgc)
library(dplyr)
library(ape)
library(pegas)

#Climatic Zone analysis. Here we use the package kgc to assign records to the climatic zone based in the coordinates. For this, we first constructed rarefaction curves per climatic zone to explore the completeness of each regions. We use the data that wasn't filtered a minimum number of records per BIN, this is the same data used for estimations of species "BIN" richness. 

data_SR <- data_geo

data_SR <- data_SR[, c("recordID","bin_uri","phylum_taxID", "country","phylum_name","class_taxID","class_name","order_taxID","order_name","family_taxID","family_name","subfamily_taxID", "subfamily_name","genus_taxID","genus_name", "species_taxID","species_name","nucleotides","lat","lon")]

#Using the package kgc to assign climatic zone classification to the coordinates. 
data_CZ <- data_SR

data_CZ$rndCoord.lat <- RoundCoordinates((data_CZ$lat), res = "fine", latlong = "lat")

data_CZ$rndCoord.lon <- RoundCoordinates((data_CZ$lon), res = "fine", latlong = "lon")
  
data_CZ$ClimateZ <- LookupCZ(data_CZ , res = 'fine')

data_CZ$ClimateZ <- as.character(data_CZ$ClimateZ)

#Checking for records missing information and filtering those out. 

sum(data_CZ$ClimateZ == "Climate Zone info missing")

data_CZ_missing <- dplyr::filter(data_CZ, data_CZ$ClimateZ == "Climate Zone info missing")

data_CZ <- dplyr::filter(data_CZ, data_CZ$ClimateZ != "Climate Zone info missing")

#Running script again to search the ones missing information and saving the ones that were assigned. 
data_CZ_missing$rndCoord.lat <- RoundCoordinates((data_CZ_missing$lat), res = "fine", latlong = "lat")

data_CZ_missing$rndCoord.lon <- RoundCoordinates((data_CZ_missing$lon), res = "fine", latlong = "lon")

data_CZ_missing$ClimateZ <- LookupCZ(data_CZ_missing, res = "fine")

data_CZ_missing2 <- filter(data_CZ_missing, data_CZ_missing$ClimateZ == "Ocean")

data_CZ_missing <- filter(data_CZ_missing, data_CZ_missing$ClimateZ != "Ocean")

data_CZ_missing2 <- data_CZ_missing2[, c(-21, -22, -23, -24)]
colnames(data_CZ_missing2) <- c("recordID" , "bin_uri"  , "phylum_taxID", "country" , "phylum_name" , "class_taxID" ,"class_name", "order_taxID","order_name", "family_taxID",    "family_name", "subfamily_taxID", "subfamily_name",  "genus_taxID","genus_name", "species_taxID","species_name",    "nucleotides", "rndCoord.lat","rndCoord.lon"  )

data_CZ_missing2$ClimateZ <- LookupCZ(data_CZ_missing2, res = "fine")

data_CZ <- full_join(data_CZ, data_CZ_missing) 

#The package allows to check uncertainty associated to the classification. I noticed some Tundra records and others in the tropic that seemed conflicting after producing a map. I recommend running this step and also producing a map to check for potential errors. I decided to manually correct some Polar records with high uncertainty and also the ones in Costa Rica. High uncertainty was defined as > 0.5. The package provides a new column with the more likely assignation. If the possible Climatic Zone was from the same zone (ej Aw instead of Am) I didn't changed it. I only corrected the ones whit high uncertainty and a complete change of climatic zone (E to D). 

data_CZ2 <- data.frame(data_CZ, CZUncertainty(data_CZ))

#creating an intermediate variable to avoid having to run again the previous step. 
data_CZ2.1 <- data_CZ2

#First, keeping the records that don't need revision. 

data_CZ2_f <- data_CZ2.1 %>%
  filter(!(uncertainty > 0.5 & country =="Costa Rica"))

data_CZ2_f <- data_CZ2_f %>%
  filter(!(uncertainty > 0.5 & ClimateZ =="ET"))

#Second, filtering the ones that need to be checked. 

data_CZ2_f1 <- data_CZ2.1 %>%
  filter((uncertainty > 0.5 & country =="Costa Rica"))

data_CZ2_f2 <- data_CZ2.1 %>%
  filter((uncertainty > 0.5 & ClimateZ =="ET"))

#write csv file and check manually in Excel. 

write_csv(data_CZ2_f1, "data_CZ2_filter.csv")
write_csv(data_CZ2_f2, "data_CZ2_filter2.csv")

#read files corrected and join in one data frame. 

data_CZf1 <- read.csv("data_CZ2_filter.csv")
data_CZf2 <- read.csv("data_CZ2_filter2.csv")

data_CZ_corrected <- bind_rows(data_CZ2_f, data_CZf1, data_CZf2)

#updating CZ object for downstream analysis. 

data_CZ <- data_CZ_corrected
  
#subsetting the data in the tropical area 
tropical <- dplyr::filter(data_CZ, grepl("^A",ClimateZ))
class(tropical)

#count the number of records per BIN 
bin_tropical <-tropical %>%
  group_by(bin_uri) %>%
  tally(length(recordID))

#filtering BINs with less than 10 records
bin_tropical10 <- filter(bin_tropical, n>= 10 )

tropical10 <- left_join(tropical,bin_tropical10, by= "bin_uri")

tropical10 <- tropical10 %>%
  filter(!is.na(n))


ggplot(data = bin_tropical)+
  geom_histogram(mapping = aes(x = n), color = "black", alpha = 0.2, binwidth = 5)+
  labs(x = "Number of Sequences per BIN", y = "Frequency", title = "Frequency Histogram of the Number of Sequences per BIN")

#Continental
continental <- dplyr::filter(data_CZ, grepl("^D",ClimateZ))

##count the number of records per BIN 
bin_continental <-continental %>%
  group_by(bin_uri) %>%
  tally(length(recordID))

#filtering BINs with less than 10 records
bin_continental10 <- filter(bin_continental, n>= 10 )

continental10 <- left_join(continental,bin_continental10, by= "bin_uri")

continental10 <- continental10 %>%
  filter(!is.na(n))

#Repeating steps for the other Climatic Zones. This code is repeated for each zone because I had to explore each subset independently. 
#Dry Zone
dry <- dplyr::filter(data_CZ, grepl("^B",ClimateZ))
bin_dry <-dry %>%
  group_by(bin_uri) %>%
  tally(length(recordID))

bin_dry10 <- filter(bin_dry, n>= 10 )
dry10 <- left_join(dry,bin_dry10, by= "bin_uri")
dry10 <- dry10 %>%
  filter(!is.na(n))

#Temperate Zone

tempe <- dplyr::filter(data_CZ, grepl("^C",ClimateZ))

bin_tempe <- tempe %>%
  group_by(bin_uri) %>%
  tally(length(recordID))

bin_tempe10 <- filter(bin_tempe, n>= 10 )
tempe10 <- left_join(tempe, bin_tempe10,by= "bin_uri")
tempe10 <- tempe10 %>%
  filter(!is.na(n))

#Polar Zone
polar <- dplyr::filter(data_CZ, grepl("^E",ClimateZ))

bin_polar <- polar%>%
  group_by(bin_uri) %>%
  tally(length(recordID))

bin_polar10 <- filter(bin_polar, n>= 10 )
polar10 <- left_join(polar, bin_polar10,by= "bin_uri")
polar10 <- polar10 %>%
  filter(!is.na(n))


#Exploring the completeness of sampling in the climatic zones as BINs will be used as a molecular proxy for species identities. First, creating a community object to use the function rarecurve. Hera I used the records per BIN count object created in the stpes before. 

tropical_comm <- pivot_wider (bin_tropical, names_from = bin_uri, values_from = n)

tropical_rc <- rarecurve(tropical_comm)


continental_com <- pivot_wider(bin_continental, names_from = bin_uri, values_from = n)

continetal_rc <- rarecurve(continental_com)


dry_com <- pivot_wider(bin_dry, names_from = bin_uri, values_from = n)

dry_rc <- rarecurve(dry_com)


tempe_com <- pivot_wider(bin_tempe, names_from = bin_uri, values_from = n)

tempe_rc <- rarecurve(tempe_com)


polar_com <- pivot_wider(bin_polar, names_from = bin_uri, values_from = n)

polar_rc <- rarecurve(polar_com)

#Running the function for the output object of the function rarecurve for all the zones in one graphic. 

out_comm <- as_tibble_rc(c( tropical_rc, dry_rc, tempe_rc, continetal_rc,polar_rc))

#Plotting the results using ggplot
mytheme = list(
  theme_classic()+
    theme(panel.background = element_blank(),strip.background = element_rect(colour=NA, fill=NA),panel.border = element_rect(fill = NA, color = "black"),
          axis.text=element_text(face="bold"),axis.title = element_text(face="bold"),plot.title = element_text(face = "bold", hjust = 0.5,size=13))
)

rarefaction_CZ <- ggplot(data = out_comm, mapping = aes(x = Sample_size, y = Species, color = Site))+
  geom_line()+ geom_path(size = 1.5)+
  labs( y = "Molecular Clusters Counts")+
  labs( x = "Sequences")+
  labs(title = "Rarefaction Curves")+
  labs( c(  "Tropical",  "Dry",  "Warm","Continental", "Polar"))+
  guides(color=guide_legend("Climatic Zones"))+
  mytheme

rarefaction_CZ  + scale_colour_viridis_d(direction =1, labels = c( "Tropical/Equatorial Zone", "Arid/Dry Zone",  "Warm/Mild Temperate", "Continental" ,"Polar Zone"))

#amount of data per CZ for the nucletoide diversity analysis. 
tropical_data <- get_data_list(tropical10)
tempe_data <- get_data_list(tempe10)
dry_data <- get_data_list(dry10)
continental_data <- get_data_list(continental10)
polar_data <- get_data_list(polar10)

data_CZ_10 <- rbind(tropical10, continental10, dry10, tempe10, polar10)
class(data_CZ_10$ClimateZ)


#Transforming the data to sf object for plotting with tmap. 
data_CZ_10$ClimateZ <- lapply(data_CZ_10$ClimateZ, FUN= str_extract, pattern=  "^.{1}")
sf_records_CZ <- st_as_sf(data_CZ_10, coords = c("lon", "lat"), crs = 4326)
sf_records_CZ$ClimateZ <- unlist(sf_records_CZ$ClimateZ)


#Map scatter plot of the data available per Climatic Zone. Molecular Clusters with 10 or more records. 

map_CZ_records <- tm_shape(world, projection="+proj=robin")+
  tm_polygons(col="grey",border.col= "grey")+
  tm_shape(sf_records_CZ) +
  tm_bubbles(col = 'ClimateZ' ,border.col = "grey", size = 0.15, palette = "viridis")+
  tm_layout (frame = FALSE)
  
map_CZ_records

#Saving file
tmap_save(map_CZ_records, filename = "CZ_records_parasitoids_July2023")

#Saving amount of data per Climatic Zones 
seq_perCZ <- rbind(continental_data, tropical_data, dry_data, tempe_data, polar_data)
seq_perCZ$ClimaticZ <- c("continental", "tropical", "dry", "tempe", "polar")

write_csv(seq_perCZ, "data_CZ_all.csv")

#Running alignments per BIN per Climatic Zone. 
tropical_list <- get_bin_list(tropical10)

tropical_alignment <- dnaStringSet_op(tropical_list)


#Calculating Nucleotide Diversity using pegas. Initially in this project I tested the performance of PopGenome and pegas nucleotide diversity calculation and the values are the same as both packages use the same function. 

pegas_nucDiv <- function(l){
  alignment <- l
  al_bin <- lapply(alignment, as.DNAbin)
  nuc_div_bin <- lapply(al_bin, nuc.div)
  return(nuc_div_bin)
}

#Calculating mean nucleotide diversity per BIN
nuc_div_tropical <- pegas_nucDiv(tropical_alignment)

#mean value per Climatic Zone. 
mean_list_cell <- function(l)
{
  l_values <- as.numeric(paste(unlist(l)))
  mean_values <- mean(l_values)
  return(mean_values)
  
}

mean_nuc_div_tropical <- mean_list_cell(nuc_div_tropical)

#Continental
continental_list <- get_bin_list(continental10)

continental_alignment <- dnaStringSet_op(continental_list)

nuc_div_continental <- pegas_nucDiv(continental_alignment)

mean_nuc_div_continental <- mean_list_cell(nuc_div_continental)

#Dry Climate
dry_list <- get_bin_list(dry10)

dry_alignment <- dnaStringSet_op(dry_list)

nuc_div_dry <- pegas_nucDiv(dry_alignment)

mean_nuc_div_dry <- mean_list_cell(nuc_div_dry)

#Temperate Climate 
tempe_list <- get_bin_list(tempe10)

tempe_alignment <- dnaStringSet_op(tempe_list)

nuc_div_tempe <- pegas_nucDiv(tempe_alignment)

mean_nuc_div_tempe <- mean_list_cell(nuc_div_tempe)

#Polar climate 
polar_list <- get_bin_list(polar10)

polar_alignment <- dnaStringSet_op(polar_list)

nuc_div_polar <- pegas_nucDiv(polar_alignment)

mean_nuc_div_polar <- mean_list_cell(nuc_div_polar)



#Getting values of nucleotide diversity per molecular cluster per Climatic Zone in a table and also developing a violin plot of the distribution of values and mean. 
#Saving values per CZ in a data frame

nuc_div_polar_df <- as.data.frame(do.call(rbind, nuc_div_polar))
nuc_div_continental_df <- as.data.frame(do.call(rbind, nuc_div_continental))
nuc_div_tropical_df <- as.data.frame(do.call(rbind, nuc_div_tropical))
nuc_div_dry_df <- as.data.frame(do.call(rbind, nuc_div_dry))
nuc_div_tempe_df <- as.data.frame(do.call(rbind, nuc_div_tempe))

data_CZ_nucDiv <- c(nuc_div_tropical_df,nuc_div_dry_df, nuc_div_tempe_df, nuc_div_continental_df, nuc_div_polar_df)

names(data_CZ_nucDiv) <- c("Tropical", "Arid/Dry", "Warm/Mild Temperate", "Continental","Polar")

data_CZ_nucDiv <- as.data.frame(do.call(rbind, data_CZ_nucDiv))

data_CZ_nucDiv <-rownames_to_column(data_CZ_nucDiv, "CZ")

#Transforming data frame to longer form for graphic 
data_CZ_nucDiv_L <- data_CZ_nucDiv %>% pivot_longer(cols= 2:2014,
                    names_to='Climatic Zone',
                    values_to='Nuc.Div')

data_CZ_nucDiv_L$CZ <- factor(data_CZ_nucDiv_L$CZ, levels = c("Tropical", "Arid/Dry", "Warm/Mild Temperate", "Continental", "Polar"))

data_CZ_nucDiv_L$Nuc.Div <- sqrt(data_CZ_nucDiv_L$Nuc.Div)

#Violing plot
viol_plot <- ggplot(data_CZ_nucDiv_L, aes(x=CZ, y=Nuc.Div, fill=CZ)) +
  geom_violin(width=1.4) +
  geom_boxplot(width=0.1, color="grey", alpha=0.2) +
  scale_fill_viridis(discrete = TRUE) +
  stat_summary(fun="mean", geom="point", shape=20, size=8, color="red", fill="red")+
  ylab("Mean Nucleotide Diversity per Molecular Cluster") +
  xlab("Climatic Zone")+
  theme_bw() 

viol_plot + theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 20))  

viol_plot 


#Rarefaction analysis 
#The lower number of BINs and sequences among the 4 Climatic Zones are 162 and 4939 respectively. To adress the efect of the amount of data in the values calculated I decided to filter BINs with more than 150 records in the continental/tropical/temperate area data set to do the rarefaction analysis from this dataset. The polar area has 150 but this I decided to try the one with 162. In the Dry Climatic Zone there is a BIN wit 815 records and then all the others are 150 and less. That's why I decided to filter by 150.  

#Filterin BINs with more than 150 records 
tropical150 <- filter(tropical10, n < 150)
continental150 <- filter(continental10, n <150)
tempe150 <- filter (tempe10, n < 150)

#Running the alignments and mean Nucleotide Diversity calculation. This section needs to be run for each CZ and the results need to be saved independently. 

tropical150_list <- get_bin_list(tropical150)

tropical150_alignment <- names_DNAStringSet(tropical150_list) # dnaStringSet_op


#Loop for rarefaction 
nuc_div_tropical150_df <- NULL

nuc_div_tropical162 <- NULL

mean_nuc_div_tropical162 <- NULL

  for (i in 1:100){
    
    sample_alignment<- sample(tropical150_alignment, 162, replace = TRUE) # sample equal number of BINS like in the tropical area
    nuc_div_tropical162[[i]] <- pegas_nucDiv(sample_alignment)
    mean_nuc_div_tropical162[[i]] <- mean_list_cell(nuc_div_tropical162[[i]])
    
    nuc_div_tropical150_df[[i]] <- as.data.frame(mean_nuc_div_tropical162[[i]])#saving the results as a data frame 

  }

nuc_div_tropical162_df <- as.data.frame(matrix(unlist(nuc_div_tropical150_df)), nrow=length(nuc_div_tropical150_df))


rare_CZ <- cbind(nuc_div_continental162_df, nuc_div_tropical162_df, nuc_div_tempe162_df )
colnames(rare_CZ) <- c("Continental", "Tropical", "Warm/Mild Temperate")

rare_CZ <-  pivot_longer(rare_CZ, cols=c('Continental', 'Tropical', 'Warm/Mild Temperate'), names_to='Climatic Zone',
                    values_to='Mean Nucleotide Diversity') 
rare_dry <- data.frame("Arid/Dry Zone", mean_nuc_div_dry )
names(rare_dry) <- c('Climatic Zone', 'Mean Nucleotide Diversity')

rare_CZ <- rbind(rare_CZ, rare_dry)

#Saving values of nucleotide diversity per Climatic Zone

CZ_vals <- c(mean_nuc_div_continental, mean_nuc_div_dry, mean_nuc_div_tempe, mean_nuc_div_tropical, mean_nuc_div_polar)

CZ_vals_df <- data.frame(name=c("Continental", "Arid/Dry Zone","Warm/Mild Temperate", "Tropical", "Polar/Tundra"), value=CZ_vals)
names(CZ_vals_df) <- c('Climatic Zone', 'Mean Nucleotide Diversity')

#Boxplot of the rarefaction analysis. 
rare_CZ_boxplot <-
  ggplot( rare_CZ, aes(x=`Climatic Zone`, y=`Mean Nucleotide Diversity`, fill=`Climatic Zone`)) +
  geom_boxplot() +
  scale_fill_viridis(direction = -1, discrete = TRUE, alpha=0.6, limits = c("Continental" , "Warm/Mild Temperate", "Arid/Dry Zone", "Tropical")) +
  geom_jitter(color="black", size=0.4, alpha=0.9) +
  theme(axis.text = element_text(size = 20))+
  theme_bw()
  #theme_ipsum() 
  #theme(
   # legend.position="none",
   # plot.title = element_text(size=11)) +
  #ggtitle("A boxplot with jitter") +
  #xlab("")
rare_CZ_boxplot


rare_CZ_boxplot <- rare_CZ_boxplot + geom_point(data=CZ_vals_df, aes(x=`Climatic Zone`, y=`Mean Nucleotide Diversity`, fill=`Climatic Zone`), color="red", size=3,  shape=17)

rare_CZ_boxplot2 <- rare_CZ_boxplot + theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 20))  +
  theme(legend.title = element_text(size = 14)) 

rare_CZ_boxplot2

