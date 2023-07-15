##################
#Spatial Re-analysis following Miraldo et al 2016. This section includes filtering steps for the amount of data per cell, rarefaction analysis and spatial autocorrelation. 

#Filtering to select only BINs that have 10 or more sequences to account for the effect of data. Run the GD script filtering for BINs with more than 10 sequences  

data_df_GD10 <- filter(data_df_GD2, n >= 10) #`seqs/bin/cell`

#Sumarizing the data after filtering step. 

total_record_count10 <- data_df_GD10 %>%
 dplyr::count(cell, bin_uri)

total <- total_record_count10 %>%
  group_by(cell) %>%
  dplyr::count(length(bin_uri))
total <- total [, c(1,3)]

colnames(total) <- c("cell", "bins")

record_cell10 <- data_df_GD10 %>%
  group_by(cell) %>%
  dplyr::summarize(count = length(unique(recordID)))

colnames(record_cell10) <- c("cell", 'sequences')

data_GD10_table <- left_join(record_cell10, total, by = "cell")

write_csv(data_GD10_table, "cell_info_dataGD10.csv")

#Adding the number of BINs and sequences per BIN per cell to the data frame. 
data_df_GD10 <- data_df_GD10 [, c(-22, -23)] #remove previous columns
data_df_GD10 <- left_join(data_df_GD10, data_GD10_table, by = "cell")

which(data_GD10_table$bins < 3) #Cells with less than 3 BINs

#Updating grid for cell that have BINs with 10 or more sequences.

grid_more10 <- dgcellstogrid(grid,data_df_GD10$cell,frame=TRUE,wrapcells=TRUE)
length(unique(grid_more10$cell))

#Creating a list of data frames. Each data frame will have the records of one cell. Having this information in a list will be important to handle a hundred of cells. 

cell_list10 <- lapply(unique(data_df_GD10$cell), function(x) data_df_GD10[data_df_GD10$cell == x,])


#Number of cells in the dataset
cell_number10 <- length(unique(data_df_GD10$cell))

df_cell_more10 <- lapply(cell_list10, as.data.frame)

#List of BINs per cell to run alignment
bin_list_df_cell_more10 <- lapply(df_cell_more10, get_bin_list)


#Retrieving the BIN_id for all the records per latitudinal band/Cell. It will be used to name the FASTA files with the alignments. 

BIN_id_cell_more10 <- lapply(bin_list_df_cell_more10, function(x) { sapply(x, function(y)(unique(y$bin_uri)))})

#Retrieving the number of records per latitudinal band/Cell

bin_number_df_cell_more10 <- lapply(bin_list_df_cell_more10, function(x) length(x))
 
#Function to name the records DNAStringset object with the Record ID and perform the Multiple Sequence Alignment for the BINs in one latitudinal band/Cell. Also performs the MSA.

Alignment_df_cell_more10 <- lapply(bin_list_df_cell_more10, dnaStringSet_op)

#Iteration over each latitudinal band/Cell. The following loop creates a directory for each cell/latitudinal band and saves the FASTA files of each BIN. It allows an easy access to the alignments per latitudinal band/cell for the nucleotide diversity analysis. 

cell_id10 <- lapply(bin_list_df_cell_more10, function(x){sapply(x, function(y)       unique(y$cell))})
cell_id10<- lapply(cell_id10, function(x) x[[1]])

for (i in seq( from = 1, to = cell_number10, by = 1)) {
  cell_id <- lapply(bin_list_df_cell_more10, function(x){sapply(x, function(y)       unique(y$cell))})
  cell_id5 <- lapply(cell_id10, function(x) x[[1]])
  work.dir <- getwd()
  newdir <- paste0(cell_id10[[i]])
  dir.create(newdir)
  setwd((file.path(newdir))) 
  write_fasta(Alignment_df_cell_more10[[i]], BIN_id_cell_more10[[i]])
  setwd((file.path(work.dir)))
}

#Reading fasta files for PoPGenome
files_cell_fasta10 <- as.list(list.dirs("/Users/jessica/Documents/MSc Research Project/June2021/April2022/Parasitoids/Parasitoids_redone", full.names = FALSE, recursive = FALSE))

Nuc_div_df_cell10 <- lapply(files_cell_fasta10, Nuc_div_function)

#Function to calculate Nucleotide diversity per site by dividing Pi by the number of sites. 

Nuc_div_df_cell10 <- lapply(Nuc_div_df_cell10, get_mean_Pi)

mean_by_cell_pi10 <- foreach (i = 1:length(Nuc_div_df_cell10)) %do%
  mean_by_cell(Nuc_div_df_cell10[[i]])

cell_number_pi10 <- sapply( Nuc_div_df_cell10 , function(x) (unique(x$cell)))                         
mean_by_cell_final10 <- do.call(rbind, Map(data.frame, A=mean_by_cell_pi10, B=cell_number_pi10))
colnames(mean_by_cell_final10) <- c("mean_by_cell_pi10", "cell")
class(mean_by_cell_final10$mean_by_cell_pi10)

grid_more10 <- left_join(grid_more10,mean_by_cell_final10,by = "cell")

data_GD10_table$cell <- as.character(data_GD10_table$cell)
grid_more_10.3 <- left_join(data_GD10_table, grid_more10, by = "cell")

#Filter cells with less than 3 BINs.
grid_more_10.3 <- grid_more_10.3 %>%
  filter(bins >= 3)
length(unique(grid_more_10.3$cell))

#Plot of the cells using tmap package
sf_grid_more10 <- st_as_sf( grid_more_10.3, coords = c("long", "lat"), crs = 4326) %>%
  group_by(group) %>%
  dplyr::summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 
class(sf_grid_more10)
plot(sf_grid_more10)

sf_grid_more10 <- left_join(sf_grid_more10, grid_more_10.3, by= "group", keep = FALSE)

library(tmap)
data("World")

world <- World %>% 
  filter(sovereignt != "Antarctica")


mapmore10.3_parasitoids2 <- tm_shape(world, projection="+proj=robin")+
  tm_polygons(col="grey",border.col= "grey")+
  tm_shape(sf_grid_more10.3, projection="+proj=robin") + tm_polygons(col = "mean_by_cell_pi10",  palette = "viridis", style = "quantile", n = 4, border.col = NULL, border.alpha = 0.4, title = "Nucleotide Diversity")+
  tm_layout(legend.position = c("left","bottom"),
            legend.bg.color = "white", frame = FALSE)
mapmore10.3_parasitoids2

tmap_save(mapmore10.3_parasitoids2, filename = "parasitoids_nuc.div_10.3.tiff")

#Correlation Analysis of Nucleotide diversity values using the whole dataset whit a minimum of 2 sequences per BIN and the one filtering for 10 as a minimum number of sequences. 

mean_cell_correlation10 <- inner_join(mean_by_cell_final, mean_by_cell_final10, by ="cell")

#square transformation
mean_cell_correlation10$mean_by_cell_pi <- sqrt(mean_cell_correlation10$mean_by_cell_pi)

mean_cell_correlation10$mean_by_cell_pi10 <- sqrt(mean_cell_correlation10$mean_by_cell_pi10)

cor_lat<- cor.test(mean_cell_correlation10$mean_by_cell_pi, mean_cell_correlation10$mean_by_cell_pi10)

install.packages("ggpubr")
library(ggpubr)

ggscatter(mean_cell_correlation10, x = "mean_by_cell_pi10", "mean_by_cell_pi", add = "reg.line", conf.int = TRUE, cor.coef = TRUE, cor.method = "pearson", xlab = "Nucleotide diversity per cell/>10 sequences per BIN", ylab = "Nucleotide diversity per cell/>2 sequences per BIN")

plot(mean_cell_correlation10$mean_by_cell_pi, mean_cell_correlation10$mean_by_cell_pi10, pch = 19, col = "dark grey", xlab = "Nucleotide diversity per cell/>2 sequences per BIN", ylab = "Nucleotide diversity per cell/>10 sequences per BIN")

abline(lm(mean_cell_correlation10$mean_by_cell_pi10 ~ mean_cell_correlation10$mean_by_cell_pi), col = "red", lwd = 3)


#####
#Rarefaction analysis doing a master alignment per BIN/per Cell and sampling sequences from the alignment.  
#Convert alignment to a DNAstringset
list_string <- lapply(Alignment_df_cell_more10, function(x) { sapply(x, function(y) (DNAStringSet(y)))})

length(list_string)

#Create a list of data frame with each unique bin as one element of the list
#list_vec_seq <- lapply(list_string[[1]], function(x) data.frame("sequences" = x, "processid" = names(x), row.names = NULL))


list_vec_seq <- NULL

 for(i in 1:length(list_string)){
  var <- list_string[[i]]
   list_vec_seq[[i]] <-lapply(var, function(x) data.frame("sequences" = x, "processid" = names(x), row.names = NULL))
 }

#Function 2: Function to name the sequences using the recordID and generate DNAstringset
names_DNAS<- function(df){
  data <- df #the argument is a data frame
  
  data_processid <- unique(data$processid) #get the processid for naming the sequences. Since I will be sampling with replacement use the unique processid
  
  dnaStringSet_data <-  DNAStringSet(data$sequence) #transforming nucleotide sequences to a DNAString format
  
  names(dnaStringSet_data) <- data_processid #name the sequences with the processid
  
  return(dnaStringSet_data)
  
}

#Defining number of iterations

library(pegas)
N <- 1000
nuc_div_df <- NULL

bin_number <- lapply(bin_list_df_cell10, function(x) { sapply(x, function(y)(unique(y$bin_uri)))}) #saving BIN number for the final file
bin_number <- lapply(bin_number, as.data.frame)

rarefaction_func <- function(l){
  
  bin_list_df_cell10_test <- l
  
  for (i in 1:N){
    
    sample10 <- lapply(bin_list_df_cell10_test, slice_sample, n = 10, replace = TRUE) # sample 10 records per BIN
    data_processid <- lapply(sample10, function(x) (x$processid))
    data_processid <- lapply(data_processid, as.data.frame)
    alignment_all <- lapply(sample10, names_DNAS_string_Mehra)
    al_bin <- lapply(alignment_all, as.DNAbin) #transforming to DNAbin format
    nuc_div_bin <- lapply(al_bin, nuc.div)#calculating nucleotide diversity using pegas
    nuc_div_df[[i]] <- as.data.frame(matrix(unlist(nuc_div_bin)), nrow=length(nuc_div_bin),byrow=TRUE)#saving the results as a data frame 
    final_df <- bind_cols(nuc_div_df)# adding the population number to the data frame 
    
  }
  
  name_col <- as.character(1:1000) #Naming the columns in the iterations dataset values
  colnames(final_df) <- c( name_col)# This needs to be adjusted depending on the number of iterations. Maybe we don't need the column name to be this specific, it is just an idea.
  
  return(final_df)
}


#Running the iteration for all the cells with parallel computing  
s <- system.time({
  rarefaction_list_cell <- parallel::mclapply(list_vec_seq, rarefaction_func, mc.cores = 6)
})
rarefaction_list_cell <- lapply(list_vec_seq, rarefaction_func)


base::saveRDS(rarefaction_list_cell, "rarefaction_list_cell.Rds")

rarefaction_list_cell <- readRDS("rarefaction_list_cell.Rds")

#Naming each data frame in the list with the number of the cell 
names(rarefaction_list_cell) <- cell_id10

#Getting the mean by row, which is the mean of the iterations/per BIN/per cell

mean_by_row_rare <- function(df){
  y <- df
  name_col <- as.character(1:1000) #Naming the columns in the iterations dataset values
  colnames(y) <- c( name_col)
  w <- y %>% mutate( m = rowMeans(y [  , c(1:1000)] ))
  return(w)
}

#Getting the mean by column. That means the value per cell for each iteration. 
mean_by_column <- function(df){
  
  y <- df 
  w <- colMeans(y)
    return(w)
}

r <- rarefaction_list_cell[[1]]
r2 <- colMeans(r)

rarefaction_list_cell_mean <- lapply(rarefaction_list_cell, mean_by_row_rare)


#Mean by column 
rarefaction_list_cell_mean_column <- lapply(rarefaction_list_cell, mean_by_column)

#Mean of the rarefaction 
rarefaction_list_cell_mean2 <- lapply(rarefaction_list_cell_mean_column, mean)

rarefaction_mean2 <- bind_rows(rarefaction_list_cell_mean_column)

#Using transpose (t) function and wrapping around a data frame to transform columns to rows. Now each row has the 1000 iterations of the one cell, and each column the mean per cell for one iteration using 10 sequences per BIN

rarefaction_mean2 <- as.data.frame(t(rarefaction_mean2))

rarefaction_mean2<- rownames_to_column(rarefaction_mean2)
 colnames(rarefaction_mean2) <- c("cell", 1:1000)
 
 #Transforming the data for plotting.
 
 rarefaction_mean2_longer <- rarefaction_mean2 %>%
   pivot_longer(!cell)
 
 rarefaction_mean2_longer$cell <- as.factor(rarefaction_mean2_longer$cell)
 
 plot(rarefaction_mean2_longer$cell, rarefaction_mean2_longer$value, pch = 19, col = "dark grey", xlab = "cell", ylab= "Nucleotide diversity 10 seq/BIN")
 
 
 #Transforming the list with the mean values of the 1000 iterations.
 rarefaction_list_cell_mean2 <- as.tibble(rarefaction_list_cell_mean2)
 
 rarefaction_list_cell_mean2 <- as.data.frame(t( rarefaction_list_cell_mean2))
 
 rarefaction_list_cell_mean2 <- rownames_to_column(rarefaction_list_cell_mean2)
 colnames(rarefaction_list_cell_mean2) <- c("cell", "mean")
 
 #Adding the mean of the iterations and plotting the regression line of the mean of the iterations and the mean of the nucleotide diversity per cell using BINs with >=10 sequences.  
 rarefaction_mean2_longer <- left_join(rarefaction_mean2_longer, rarefaction_list_cell_mean2, by = "cell")
 
corre_rare <- left_join(rarefaction_mean2_longer, mean_by_cell_final10, by = "cell")

#removing values of zero for the correlation analysis with the log transformation. I do the log transformation as in the Miraldo et al 2016 paper.
corre_rare <- corre_rare %>%
  filter(mean != 0)

plot(log(corre_rare$mean_by_cell_pi10), log(corre_rare$value), pch = 19, col = "dark grey", xlab = "log Nucleotide Diversity >= 10 seq/BIN", ylab= "log Nucleotide Diversity of 10 random sequences")

abline(lm(log(corre_rare$mean_by_cell_pi10) ~ log(corre_rare$mean)), col = "red", lwd = 3)

# Spatial autocorrelation
# Code to test for spatial autocorrelation in the data.

library(spdep)
library(spatialreg)
#library(ncf)
install.packages("pgirmess")
library(pgirmess)

#### 
#Spatial Autocorrelogram using the package pgirmess
cell_data10 <- left_join (cell_centre_df_more10, mean_by_cell_final10, by="cell" )
coords <- cell_data10[, c(2,3)] #using cell center of cells. 

values_nucdiv <- cell_data10 [,4]

co1 <- pgirmess::correlog(coords,values_nucdiv, method = "Moran" )
class(co1)
#co1 <- as.data.frame(co1)
plot(co1)

#Graphic
plot(co1[,1], co1[,2], type='b', cex=0.8, ylab= 'Moran I', xlab='Distance class (km)', ylim=range(c(co1[,2], co1[,2])))

#Highlight significant values. 
points(co1[,1], co1[,2],col = ifelse(co1[,3] < 0.05, "red", "black"), cex=0.8, pch = 16)

abline(h = 0, lty = 2)

# Correlogram: A plot of Moran's I against spatial distance classes
## Quantifies the degree of spatial autocorrelation across distances


#Build a neighbors list based in regions with contiguous boundaries. Using the sf polygons for the grid cells.
install.packages("spdep")
library(spdep)

sf_grid_more10.3 <- st_as_sf(grid_more_10.3, coords = c("long", "lat"), crs = 4326) %>%
  group_by(group) %>%
  dplyr::summarise(geometry = st_combine(geometry)) %>%
  st_cast("POLYGON") 
class(sf_grid_more10.3)
plot(sf_grid_more10.3)

#Steps to add values of nucleotide diversity and lat long to the sf object for the neighbor analysis. 
grid_fil <- grid_more_10.3 [, c(9, 10)]

sf_grid_more10.3 <- left_join(sf_grid_more10.3, grid_fil , by ="group" )
sf_grid_more10.3 <- left_join(sf_grid_more10.3, cell_correlation, by = "mean_by_cell_pi10")

data_moran <- sf_grid_more10[, c(8,9)]
data_moran <- distinct(sf_grid_more10.3)
data_moran <- left_join(data_moran, cell_centre_df, by = "cell")
data_moran$mean_by_cell_pi10 <- sqrt(data_moran$mean_by_cell_pi10)


#neighbors #Removing cells 1297 ans 199 as those were not included in the modelling(This step is optional and comes from downstream analysis)
data_moran_f <- data_moran  %>%
  filter(!(cell == "1297" | cell == "199"))
neigh <- poly2nb(data_moran_f$geometry)

#Steps described in the package vignette. 
w_obj <- nb2listw(neigh, style = 'B', zero.policy = TRUE)
 
m_test <- moran.test(data_moran$mean_by_cell_pi10, w_obj, zero.policy = TRUE, randomisation = TRUE)

print(m_test)

moran.plot(data_moran$mean_by_cell_pi10, w_obj, zero.policy = TRUE)

moran.mc(data_moran$mean_by_cell_pi10, w_obj, zero.policy = TRUE, nsim = 10000)

#autocorrelation of the residuals. This step is after running the models but I wanted to add it to the same section of Spatial Autocorrelation.  
Moran_lm <- lm.morantest(residuals_nd1.3, w_obj, zero.policy = TRUE)

print(Moran_lm)

moran.mc(residuals_nd1.2$residuals, w_obj, zero.policy = TRUE, nsim = 10000)

#Putting together cell centre and cell values for mapping the mean nucleotide diversity value and number of BINs

sf_data_N <- sf_grid_more10 [, c(3, 4, 11)]
sf_data_N <- as.data.frame(sf_data_N)
sf_data_N <- sf_data_N [, -4]

sf_data_N <- sf_data_N [!duplicated(sf_data_N ), ]
sf_data_N <- left_join(sf_data_N, cell_centre_df_more10, by = "cell")
sf_data_N <- st_as_sf(sf_data_N, coords = c("long", "lat"))
total$cell <- as.character(total$cell)
sf_data_N <- left_join (sf_data_N, total, by = "cell")

map_cell_centre <- tm_shape(world, projection="+proj=robin")+
  tm_polygons(col="grey",border.col= "grey")+
  tm_shape(sf_data_N)+
  tm_bubbles("bins", col = "mean_by_cell_pi10", palette = "viridis", n = 5, border.col = "black", border.alpha = .5, title.size="Number Molecular clusters", title.col="Mean Nucleotide diversity", scale = 5)+
  tm_layout (frame = FALSE)

map_cell_centre

tmap_save(map_cell_centre, filename = "Bubble map Nuc_cell.tiff")


