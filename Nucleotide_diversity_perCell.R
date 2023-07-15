###### This script can be used to calculate nucleotide diversity per BIN (molecular cluster) at a grid cell scale using a dataset of georeferenced nucletoide sequences. In this example we use BOLD data. 

library(foreach)
library (tidyverse)
install.packages('devtools')
library(devtools)
install_github('r-barnes/dggridR')
library(dggrid1R)
install.packages("PopGenome")
library(PopGenome)
library(sf)
library(tmap)
library(viridis)

#Creating a hexagonal cell grid of the world map using the package ddgrid. The function dgconstruct construct a discrete global grid system (dggs) object. The resolution can be modified to select the best cell area and spacing. I used the default parameters for this test. The topology can be modified, too; the default is Hexagon. It is posible to test the amount of data using different resolution of the hexagon. In our case res = 6 was the most appropiate one. 

grid6 <- dgconstruct(res = 6)

#Get the corresponding grid cells for each record (lat-long pair)

data_df_GD <- data_df_GD %>%
  mutate(cell = dgGEO_to_SEQNUM(grid6, data_df_GD$lon, data_df_GD$lat)$seqnum)

#Determining number of cells

length(unique(data_df_GD$cell))

#Get the number of records per BIN in each cell, the mean value, and the number of BINs per cell 

record_count <- data_df_GD %>%
  count(cell, bin_uri)

#adding the number of records per BIN per cell to the data frame
data_df_GD2 <- left_join(data_df_GD, record_count, by=c("cell", "bin_uri"))

#Number of BINs per cell 
bin_cell <- data_df_GD2 %>%
  group_by(cell) %>%
  summarize(count = length(unique(bin_uri)))

max(record_count$n)
mean(record_count$n)
min(record_count$n)
length(record_count$cell)
sum(record_count$n ==1)
one_record <- as.list(record_count$bin_uri[record_count$n == 1])

#Filtering records of only one sequence per BIN in a cell.

data_df_GD2 <- filter(data_df_GD2, n != 1) #This n represents sequences/per BIN/per Cell
 
#Saving to file this data frame as it is the one ready for downstream analyses. 
 
write_tsv(data_df_GD2, "data_df_GD2.tsv")

#Getting the total number of records per cell 
   
total_record_count <- data_df_GD2 %>%
  count(cell)

#Summarizing all values per cell

cell_values_table <- left_join(total_record_count, bin_cell, by = "cell")
colnames(cell_values_table) <- c("cell", "sequences", "bins")

write_tsv(cell_values_table, "cell_value_table.tsv")


#Get the grid cell boundaries for cells that have records using the dggrid R package. The next steps are recommended in the package's vignette.

grid <- dgcellstogrid(grid6, data_df_GD2$cell,frame=TRUE,wrapcells=TRUE)

cellcenters <- dgSEQNUM_to_GEO(grid6, in_seqnum = data_df_GD2$cell)


#Update the grid cells' properties to include the number of records in each cell
grid <- merge(grid,total_record_count,by.x="cell",by.y="cell")

class(grid)

#Get polygons for each country of the world
countries <- map_data("world")


#Creating a list of data frames. Each data frame will have the records of one cell. Having this information in a list will be important to handle a hundred of cells. 
cell_list <- lapply(unique(data_df_GD2$cell), function(x) data_df_GD2[data_df_GD2$cell == x,])

#Number of cells in the dataset
cell_number <- length(unique(data_df_GD2$cell))
cell_number

#Changing the format to a data framed as it was as a grouped data frame.
df_cell <- lapply(cell_list, as.data.frame)

#############
#Load functions needed for this section
#Function that groups the data per BIN in one Latitudinal Band/cell, making a list that includes all the records per BIN in one latitudinal band/Cell.

get_bin_list <- function(df){
  
  df_f <- df
  
  #Creating a list of data frames with the records per BIN.
  bin_list <- lapply(unique(df_f$bin_uri), function(x) df_f[df_f$bin_uri == x,])
  
  return(bin_list)
}

#Function to name the records DNAStringset object with the Record ID and perform the Multiple Sequence Alignment (MSA) for each BINs in one latitudinal band/Cell. Also performs the MSA. It takes as an argument a list of 

dnaStringSet_op <- function (l){
  bin_data <- l
  
  number_data <- length(bin_data)
  
  bin_list_recordID <- lapply(bin_data, function(x) (x$recordID))
  
  dnaStringSet1<- lapply(bin_data, function(x) DNAStringSet(x$nucleotides))#here if the column where the nucleotide sequences are stored should be changed accordingly.
  
  for (i in seq(from = 1, to = number_data, by = 1)){
    names(dnaStringSet1[[i]]) <- bin_list_recordID[[i]]}
  
  #Multiple sequence alignment using the muscle package. Running a multiple sequence alignment on each element of the dnaStringSet1 list. Using diags = TRUE with Muscle command to speed up alignment of each BIN.
  
  Alignment_df_cell <- foreach(i=1:number_data) %do%
    muscle::muscle(dnaStringSet1[[i]], maxiters = 3, diags = TRUE, gapopen = -3000)
  
  return( Alignment_df_cell)
  
}


#Steps using the function get_bin_list to group the data per BIN in one cell.
bin_list_df_cell <- lapply(df_cell, get_bin_list)

length(bin_list_df_cell)

#Retrieving the BIN_id for all the records per cell. It will be used to name the FASTA files with the alignments. 

BIN_id_cell <- lapply(bin_list_df_cell, function(x) { sapply(x, function(y)(unique(y$bin_uri)))})

#Retrieving the number of records per cell

bin_number_df_cell <- lapply(bin_list_df_cell, function(x) length(x))

#Iterating over each Cell. 
#Function to name the records DNAStringset object with the Record ID and perform the Multiple Sequence Alignment for the BINs in one latitudinal band/Cell. Also performs the MSA.

Alignment_df_cell <- lapply(bin_list_df_cell, dnaStringSet_op)

#Saving the alignment to file. 
base::saveRDS(Alignment_df_cell, file = "Alignment_df_cell.rds")

#The following loop creates a directory for each cell/latitudinal band and saves the FASTA files of each BIN. It allows an easy access to the alignments per latitudinal band/cell for the nucleotide diversity analysis. 

#Function write_fasta

#Writing alignments per BIN as a FASTA file.This function uses as arguments an alignment per BIN and a list of the BIN ID for naming the file. The output FASTA file will be named with the BIN ID. 


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



cell_id <- lapply(bin_list_df_cell, function(x){sapply(x, function(y) unique(y$cell))})
cell_id <- lapply(cell_id, function(x) x[[1]])

for (i in seq( from = 1, to = cell_number, by = 1)) {
  cell_id <- lapply(bin_list_df_cell, function(x){sapply(x, function(y)       unique(y$cell))})
  cell_id <- lapply(cell_id, function(x) x[[1]])
  work.dir <- getwd()
  newdir <- paste0(cell_id[[i]])
  dir.create(newdir)
  setwd((file.path(newdir))) 
  write_fasta(Alignment_df_cell[[i]], BIN_id_cell[[i]])
  setwd((file.path(work.dir)))
}

######

#Reading the alignments per BIN in one cell using the function readData of the package PopGenome. The FASTA files need to be in one folder that is used as the argument of the readData function. That folder must be in the working directory. This produce an object of class GENOME. When typing GENOME.class, we get some information about the main methods provided by PopGenome and how to access the results.

#Change the working directory accordingly. 
files_cell_fasta <- as.list(list.dirs("/Users/jessica/Documents/MSc Research Project/June2021/cell_analysis_final", full.names = FALSE, recursive = FALSE))

class(files_cell_fasta)

#Using the function diversity.stats to calculate nucleotide and haplotype diversities. pi = TRUE to activate Nei’s calculation of pi. The authors suggest that the nucleotide diversities have to be divided by the number of sites to obtain diversities per site. This shouldn't be a problem for us because after running coil we will have the sequences framed and trimmed to the same length of the barcode region.


#Function to get the Pi value from the package PopGenome and save it in a data frame.
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

#Iterate over all the folders with the fasta files consisting of each alignment per BIN
Nuc_div_df_cell <- lapply(files_cell_fasta, Nuc_div_function)

#Function to calculate Nucleotide diversity per site by dividing Pi by the number of sites. 

class(Nuc_div_df_cell[[1]])

get_mean_Pi <- function(df){
  x <- df
  x$Pi_per_site <- x$Pi/x$N_sites
  return(x)
}

#Iteartion over each BIN per cell.
Nuc_div_df_cell <- lapply(Nuc_div_df_cell, get_mean_Pi)


#The function below calculates the mean nucleotide diversity per cell. 
mean_by_cell <- function(df){
  mean_cell <- mean(df$Pi_per_site)
    return(mean_cell)
}

#Calculating mean nucleotide diversity per cell 
mean_by_cell_pi <- foreach (i = 1:length(Nuc_div_df_cell)) %do%
  mean_by_cell(Nuc_div_df_cell[[i]])

#getting the cell number                              
cell_number_pi <- sapply( Nuc_div_df_cell , function(x) (unique(x$cell)))       

#Saving values per cell 
mean_by_cell_final <- do.call(rbind, Map(data.frame, A=mean_by_cell_pi, B=cell_number_pi))
colnames(mean_by_cell_final) <- c("mean_by_cell_pi", "cell")

class(mean_by_cell_final$mean_by_cell_pi)

#Update the grid cells' properties to include the mean nucleotide diversity in each cell

grid_df <- left_join(grid,mean_by_cell_final,by = "cell") 

#Plot of the cells using tmap package. For that I transform the grid object to an sf object. 

sf_grid <- st_as_sf( grid_df, coords = c("long", "lat"), crs = 4326) %>%
  group_by(group) %>%
  dplyr::summarise(geometry = st_combine(geometry))%>% 
  st_cast("POLYGON") 


#adding the lat long columns as with the transformation you obtain the geometry column 
sf_grid <- left_join(sf_grid, grid_df, by= "group", keep = FALSE)

class(sf_grid)
plot(sf_grid)


data("World")

world <- World %>% 
  filter(sovereignt != "Antarctica")

#mapping using the Robin projection for the World Map. 
map_nuc_div <- tm_shape(world, projection="+proj=robin")+
  tm_polygons(col="grey",border.col= "grey")+
  tm_shape(sf_grid, projection="+proj=robin") + tm_polygons(col = "mean_by_cell_pi",  palette = "viridis", n = 4, border.col = NULL, border.alpha = 0.4, title = "Nucleotide Diversity per cell")+
  tm_layout(legend.position = c("left","bottom"),
            legend.bg.color = "white", frame = FALSE)
map_nuc_div

tmap_save(map_nuc_div, filename = "Mean_Nuc_cell.tiff")


###########
#graphics summarizing data per cell 
cell_values_sum <- left_join(cell_values_table, mean_by_cell_final, by = "cell")
cell_values_sum <- left_join(cell_values_sum, cell_centre_df, by = "cell")


png('cell_sequences_bins.png', width = 7, height = 7, units = 'in', res = 400, bg = 'transparent')

ggplot(data = cell_values_sum, aes(x=bins, y=sqrt(mean_by_cell_pi), size = sequences)) +
  geom_point(alpha=0.7) +
  scale_size(range = c(1.4, 19), name="Number of Sequences") +
  scale_color_viridis(discrete=FALSE ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))+
xlab("Number of BINs per Cell") + ylab("Mean Nucleotide Diversity per Cell")+
labs(color = "Mean Nucleotide Diversity")
dev.off()

Mean_pi_percell <- ggplot(data = cell_values_sum)+
  #geom_bar(mapping = aes(x = cell, y = mean_by_cell_pi), stat = 'identity')+
  geom_point(mapping = aes(x = sequences, y = sqrt(mean_by_cell_pi)))+
  labs(x = "Sequences/Cell", y = "Mean Nucleotide diversity/Cell", title = "Mean Nucleotide Diversity per cell")+
  #scale_x_discrete(breaks = Mean_percell_graph$cell[c(T,F,F)]) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  labs(size = "Sequences per Cell")
  #stat_smooth(method = "lm", formula = y ~ x)

cell_values_sum_filt <- cell_values_sum %>%
  filter(sequences < 1000)

Mean_pi_percell <- ggplot(data = cell_values_sum_filt)+
  #geom_bar(mapping = aes(x = cell, y = mean_by_cell_pi), stat = 'identity')+
  geom_point(mapping = aes(x = (sequences), y = sqrt(mean_by_cell_pi)))+
  labs(x = "Sequences/Cell", y = "Mean Nucleotide diversity/Cell", title = "Mean Nucleotide Diversity per cell < 1000 seq")+
  #scale_x_discrete(breaks = Mean_percell_graph$cell[c(T,F,F)]) +
  #theme(axis.text.x = element_text(angle = 90, vjust = 0.5))+
  labs(size = "Sequences per Cell")
#stat_smooth(method = "lm", formula = y ~ x)

Mean_pi_percell


png('cell_sequences_pi.png', width = 7, height = 7, units = 'in', res = 400, bg = 'transparent')


