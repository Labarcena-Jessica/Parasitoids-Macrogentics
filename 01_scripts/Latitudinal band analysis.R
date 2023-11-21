######
#Latitudinal Scale diversity analysis.

library(tmap)
library(terra)
library(viridis)
library(sp)
library(tidyverse)

#Grouping the data per latitudinal band. I will have a list of dataframes of the data grouped by increments of 10 degrees for BINs with more than 10 sequences. 

LatBand_data10 <- split(data_df_GD10, cut(data_df_GD10$lat, seq(-50, 90, by=10)))

data_length <- length(LatBand_data10)

#Changing the class of the data in the list as data frame. It was previously a grouped data frame. 

LatBand_data10  <- lapply(LatBand_data10 , as.data.frame)

#Get the amount of information per latitudinal band

get_data_list <- function(df){
  
  x <- df
  bins <- length(unique(x$bin_uri))
  seqs <- length(x$recordID)
  data_final <- data.frame(bins, seqs)
  return(data_final)
}

data_perLatBand10 <- lapply(LatBand_data10 , get_data_list)

data_perLatBand_final10 <- as.data.frame(do.call(rbind, data_perLatBand10))

#Iteration using lapply through each latitudinal band

bin_list_latitude_all10 <- lapply(LatBand_data10, get_bin_list)

#Retrieving the BIN_id for all the records per latitudinal band/Cell. It will be used to name the FASTA files with the alignments. 

BIN_id_lat10 <- lapply(bin_list_latitude_all10, function(x) { sapply(x, function(y)(unique(y$bin_uri)))})

#Retrieving the number of records per latitudinal band/Cell

bin_number_lat_test_all10 <- lapply(bin_list_latitude_all10, function(x) length(x))

#Iterating over each latitudinal band/Cell

Alignment_Latitudinal_all10 <- lapply(bin_list_latitude_all10, dnaStringSet_op)

write_rds(Alignment_Latitudinal_all10, file = "Alignment_Latitudinal_all.rds")

#Writing alignments per BIN as a FASTA file.This function uses as arguments an alignment per BIN and a list of the BIN ID for naming the file. The output FASTA file will be named with the BIN ID. 

class(Alignment_Latitudinal_all)


write_fasta <- function(x,q){
  data_popG_lat <- x
    align_number2 <- length(x)
    
  BINs_list <- q

  for(i in seq(from = 1, to = align_number2, by= 1)) {
    data_popG_lat[[i]] = as(x[[i]], "DNAStringSet")
    
  }
  for (i in seq(from = 1, to = align_number2, by= 1) ) {
    writeXStringSet(data_popG_lat [[i]], file = paste0(BINs_list[i],'.fa'))
    
    }
}

#Iteration over each latitudinal band/Cell. The following loop creates a directory for each cell/latitudinal band and saves the FASTA files of each BIN. It allows an easy access to the alignments per latitudinal band/cell for the nucleotide diversity analysis. 

 for (i in seq( from = 1, to = data_length, by = 1)) {
   wse <- as.list(names(NucDiv_perLatBand_data10))
   work.dir <- getwd()
   newdir <- paste0(wse[[i]])
   dir.create(newdir)
   setwd((file.path(newdir))) 
   write_fasta(Alignment_Latitudinal_all10[[i]], BIN_id_lat10[[i]])
   setwd((file.path(work.dir)))
  }


#Reading the alignments per BIN in one cell using the function readData of the package PopGenome. The FASTA files need to be in one folder that is used as the argument of the readData function. That folder must be in the working directory. This produce an object of class GENOME. When typing GENOME.class, we get some information about the main methods provided by PopGenome and how to access the results.

#Using the function diversity.stats to calculate nucleotide and haplotype diversities. pi = TRUE to activate Nei’s calculation of pi. The authors suggest that the nucleotide diversities have to be divided by the number of sites to obtain diversities per site. This shouldn't be a problem for us because after running coil we will have the sequences framed and trimmed to the same length of the barcode region. Nevertheless, it will be discussed. 

#After running PopGenome, I have the values of Pi per BIN en each latitudinal band/cell. I want that information in a data frame so I can then merge that with the rest of the metadata for each record per BIN.
 
#Accessing Nei’s calculation of pi in the latitudinal band small data set.

#Function to get the Pi value from the package PopGenome and save it in a data frame. Resumes all the previous steps in one function for easy manipulation of the data.

Nuc_div_function <- function(directory){
  z <- directory
  x <- readData(directory)
  y <- diversity.stats(x, pi = TRUE)
  Nuc_div_df <- as.data.frame(y@Pi, optional = TRUE)
  Nuc_div_df <- rownames_to_column(Nuc_div_df)
  Nuc_div_df$cell <- paste0(z)
  Nuc_div_df$NSites <- x@n.sites
  colnames(Nuc_div_df) <- c("bin_uri", "Pi", "cell", "N_sites")
  return(Nuc_div_df)
}
 

 ############
#Reading the fasta files to use PopGenome. 
 
 files_lat_fasta <- as.list(list.dirs("/Users/jessica/Documents/MSc Research Project/June2021/April2022/Latitudinal_band", full.names = FALSE, recursive = FALSE))
 
 #Reading all the files per Latitudinal Band and estimating nucleotide diversity.

 Nuc_div_df_lat <- lapply(files_lat_fasta, Nuc_div_function) 

 Nuc_div_df_lat <- lapply(Nuc_div_df_lat, get_mean_Pi)
 
#Getting mean value per latitudinal band   
 mean_by_lat_pi <- foreach (i = 1:length(Nuc_div_df_lat)) %do%
   mean_by_cell(Nuc_div_df_lat[[i]])
 
#Getting the latitudinal band information and saving a csv with mean values of nucleotide diversity per latitudinal band 
 lat_number_pi <- sapply( Nuc_div_df_lat , function(x) (unique(x$cell)))
 
  mean_by_lat_final <- do.call(rbind, Map(data.frame, A=mean_by_lat_pi, B=lat_number_pi))
 colnames(mean_by_lat_final) <- c("mean_by_lat_pi", "lat_band")
  
 write_csv(mean_by_lat_final, "mean_pi_by_lat.csv")
  
#Creating a raster file with the mean values of nucleotide diversity per latitudinal band to create a map. 
 
 #Create a matrix of the values per latitudinal band to create raster file. This needs to be in the order of the latitudinal bands.
 matrix_lat <- c(mean_by_lat_final$mean_by_lat_pi)
 matrix_lat <- matrix_lat[c(14, 13, 12,  11,  10,  9,  8,7,  6,1, 2, 3, 4, 5)]
 matrix_lat <- c(matrix_lat, 0, 0, 0, 0)

lat_matrix_values <- matrix (rep(matrix_lat, each= 36), nrow = 14, ncol = 36, byrow = TRUE)

 rast_lat <- rast(nrow = 14, ncol = 36, 
             ymin = -50, ymax= 90,
             vals = lat_matrix_values)
 rast_lat
 rast_lat2 <- terra::disagg (rast_lat, fact = 15)
 
 class(rast_lat2)
 plot(rast_lat2)
 
 # Load land polygon, we will use it to mask out the marine environment.  
 #This file needs to be provided to reproduce the map
 land <- st_read('ne_50m_land/ne_50m_land.shp')
 
 landRast <- rasterize(vect(land), rast_lat2)
 
 #Mask areas with no values
 topo <- terra::mask(rast_lat2, mask = landRast)
 
 plot(topo) 
 crs <- "+proj=robin"
 
 #Re-project raster
 topo <- project(topo, crs)
 
 #mapping 
 tmap_options(max.raster = c(plot = 1e9, view = 1e9))
 
 lat_mean_map <- tm_shape(world, projection="+proj=robin")+ tm_polygons(col="grey",border.col= "grey")+ 
   tm_shape(topo, border.col = NULL)+ tm_raster(palette = "viridis", title = "Nucleotide Diversity")+tm_layout(legend.position = c("left","bottom"), legend.bg.color = "white", frame = FALSE)
 
 
 lat_mean_map 
 
 tmap_save(lat_mean_map, filename = "Nuc_div_latBand.tiff")
 
 data_perLatBand_final <- tibble::rownames_to_column(data_perLatBand_final, "lat_band")

 data_all_LatBand <- left_join(data_perLatBand_final, mean_by_lat_final, by = "lat_band") 



