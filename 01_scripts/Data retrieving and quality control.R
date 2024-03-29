install.packages("tidyverse")
install.packages("vegan")
install.packages("foreach")
library(foreach)
library (tidyverse)
library(vegan)
library(devtools) 
devtools::install_github("CNuge/coil", build_vignettes = TRUE)
library(coil)
library(parallel)

#Installing the package Biostring for manipulation of nucleotide sequences 

#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")

library(Biostrings)

################
#Data exploration. Reading the data. Set the working directory to where the data is saved.  

data_bold <- read_delim("Copy of Dirk_Hym_20231116.txt")
names(data_bold)

#Filtering for lat/long presence. We need all the records to have this information

data_bold <- data_bold %>%
  filter(!is.na(lat) & !is.na(long))

sum(is.na(data_bold$long))

#To this point there are 661 635 records for Ichneumonidae, Braconidae and Formicidae. 

#Creating a separate file for each family 

unique(data_bold$family)

ichneumonidae <- data_bold %>%
  filter(family == "Ichneumonidae")

braconidae <- data_bold %>%
  filter(family == "Braconidae")

formicidae <- data_bold %>%
  filter(family == "Formicidae")

parasitoids <- full_join(ichneumonidae, braconidae)

#Creating generic variable to make script reproducible. Just changing the input file it can run completely. 

data <- parasitoids

###

#Determining the number of records missing BIN and filtering those records. The file is now called Ichneumonidae_BIN because it has all the records bearing BINs. 

sum(is.na(data$GUID))

data_BIN <- data %>%
  filter(!is.na(GUID))

#Checking how many unique BINs are present in the data set as this will be used as a proxy for species. There are 43093 BINs

length(unique(data_BIN$GUID))

#Determining the number of records missing a sequence (0, all the records have a nucleotide sequence)

sum(is.na(data_BIN$nucraw))

#Removing records missing nucleotide sequence 

data_BIN <- data_BIN %>%
  filter(!is.na(nucraw))


############
#Exploring the distribution of sequence length 
max(data_BIN$nucraw_length)
min(data_BIN$nucraw_length)

#Histogram to visualize the distribution of the sequence length

ggplot(data = data_BIN)+
  geom_histogram(mapping = aes(x = nucraw_length), color = "black", alpha = 0.2, binwidth = 10)+
  labs(x = "Nucleotide Sequence Length", y = "Frequency", title = "Frequency Histogram")

#Filtering by sequence length. Removing sequences too short and too long. Keeping all the sequences between 500 to 700 nucleotides. 

data_BIN <- data_BIN %>%
  filter (nucraw_length >= 500 & nucraw_length <= 700)

#For the analysis of Nucleotide Diversity I need to filter all the records that have 1 sequence per BIN. 

data_GD <- data_BIN %>%
  group_by(GUID) %>%
  filter(length(processid) != 1)

#We have 334123 records for parasitoids that have more than one sequences per BIN.

# Determining the number of sequences per BIN and summarizing that into a table

Sequences_bin <- data_GD %>%
  group_by(GUID) %>%
  summarize(count = length(processid))

#writing filtered data to file. 
write_tsv(data_GD, "data_GD_NOV_2023.tsv")

############
#The next step is to check the quality of the nucleotide sequences before the alignment. The R package coil contains functions for placing COI-5P barcode sequences into a common reading frame, translating DNA sequences to amino acids and for assessing the likelihood that a given barcode sequence includes an insertion or deletion error.


#Some records have non-IUPAC symbols WYSRDV in the nucleotide sequences (here some authors suggest removing the sequences but we decided to place an N. If this character is generating any problem in the sequence the open reading frame will be shifted and filtered by coil). I further investigated those sequences as this character is not accepted to run coil.

length(grep(data_GD$nucraw, pattern = "[Y, S, R, D, V, K, M, H, W, y, s, r, d, v, k, m, h, w]\\d+"))

seq_error3 <- vector('list')

for (i in 1:length(data_GD$nucraw)) {
  ind <- grep(data_GD$nucraw, pattern = "[y, s, r, d, v, k, m, h, w]")
  seq_error3 <- data_GD[ind,]
}

#seq_error2 <- seq_error [[1]]

seq_error2$nucraw[[2]]

#Removing sequences. There are 549 with problems 

data_GD <- data_GD %>% 
  filter(!processid %in% seq_error3$processid)

#Code to run coil after removing problematic characters.

data_GD$processid <- as.character(data_GD$processid)

genetic_code_coil <- which_trans_table("Ichneumonidae")

full_coi5p_df = data.frame(matrix(ncol = 9, nrow = 0),stringsAsFactors = FALSE)

colnames(full_coi5p_df) = c("name", "raw", "framed", "was_trimmed", "align_report", "aaSeq", "aaScore", "indel_likely", "stop_codons")


output_data_all = parallel::mclapply(1:length(data_GD$processid), function(i){
  coi5p_pipe(data_GD$nucraw[i], 
             name = data_GD$processid[i], 
             trans_table = genetic_code_coil, 
             triple_translate = TRUE)}, mc.cores = 5 )


head(output_data_all)
tail(output_data_all)
length(output_data_all)

full_coi5p_df <- rbind(full_coi5p_df, flatten_coi5p(output_data_all))

write_tsv(full_coi5p_df, "coil_results_date.tsv")

#Summarizing the results of coil.Remove sequences with errors. 

sum(full_coi5p_df$stop_codons == FALSE)
sum(full_coi5p_df$stop_codons == TRUE)
sum(full_coi5p_df$was_trimmed == TRUE)

#I am going to place de framed and trimmed sequences in the working data frame (data_df_GD).

#Changing the column name in full_coi5p_df to join the frame sequences to the working data frame. 

colnames (full_coi5p_df) [1] <- "recordID"

#check class and fix if needed to be the same
class(data_df_GD$recordID)
class(full_coi5p_df$recordID)


data_df_GD <- left_join(data_df_GD, full_coi5p_df[c(1,3)], by = "recordID" ) 

data_df_GD<- as.data.frame(data_df_GD)

#Adding the sequence length to a column and plotting the distribution
data_df_GD <- data_df_GD %>%
  mutate(Sequence_Length_framed = str_count (data_df_GD$framed))

max(data_df_GD$Sequence_Length_framed)
min(data_df_GD$Sequence_Length_framed)

ggplot(data = data_df_GD)+
  geom_histogram(mapping = aes(x = Sequence_Length_framed), color = "black", alpha = 0.2, binwidth = 10)+
  labs(x = "Nucleotide Sequence Length", y = "Frequency", title = "Frequency Histogram")

###########
#A Barcode Index Number can have multiple records and sequences. A single representative sequence (centroid) will be selected from each BIN for further analysis following the methodology described by (Orton et al., 2019).This will simplify the analysis as I will select one sequence per BIN. A Centroid sequence is defined as the BIN sequence with minimum average pairwise distance to all other sequences in a given BIN.

bin_list_df_cs <- (data_df_GD [, c("recordID", "bin_uri", "species_name", "framed")])

bin_list_cs <- lapply(unique(bin_list_df_cs$bin_uri), function(x) bin_list_df_cs[bin_list_df_cs$bin_uri == x,])

#Taking information regarding each BIN including number of recordIDs within each BIN.
bin_list_df_cs <- bin_list_df_cs %>%
  group_by(bin_uri) %>%
  mutate(bin_size = (count = length(recordID)))

#I already filtered the records that have only 1 sequence per BIN so I will need to determine the centroid sequences for all the dataset. 

# Also need to find the number of unique BINs in bin_list_df
bin_numberCentroid <- unique(bin_list_df_cs$bin_uri)
bin_numberCentroid <- length(bin_numberCentroid)

# Extract record id from each BIN.
bin_list_RecordId_cs <- sapply( bin_list_cs , function(x) ( x$recordID) )

# Convert all of the sequences in the bin_list_cs to dnaStringSet format for the alignment step.
dnaStringSet1 <-
  sapply( bin_list, function(x) DNAStringSet(x$framed) )

class(dnaStringSet1)

# name DNAStringSet with the record ids.
for (i in seq(from = 1, to = bin_numberCentroid, by = 1)) {
  names(dnaStringSet1[[i]]) <- bin_list_RecordId_cs[[i]]
}

# Multiple sequence alignment using the muscle package. Running a multiple sequence alignment on each element of the dnaStringSet1 list. Using diags = TRUE with Muscle command to speed up alignment of each BIN.

BiocManager::install("muscle")
library(muscle)

alignment_framed <- foreach(i=1:bin_numberCentroid) %do%
  muscle::muscle(dnaStringSet1[[i]], maxiters = 3, diags = TRUE, gapopen = -3000)

class(alignment_framed)

# I can then convert each alignment to DNAbin format.
install.packages("ape")
library (ape)

dna_BINCentroid <-
  foreach(i=1:bin_numberCentroid) %do% as.DNAbin(alignment_framed[[i]])

# Then, I perform genetic distance determination with the TN93 model on each DNAbin list.

genetic_distance_centroid <- foreach(i=1:bin_numberCentroid) %do%
  dist.dna(dna_BINCentroid[[i]], model = "TN93", as.matrix = TRUE,
           pairwise.deletion = TRUE)

# The centroid sequence can be determined from the distance matrix alone. It is the sequence in a bin with minimum average pairwise distance to all other sequences in its BIN.

centroid_seq <- foreach(i=1:bin_numberCentroid) %do%
  which.min(rowSums(genetic_distance_centroid[[i]]))

centroid_seq <- unlist(centroid_seq)
centroid_seq <- names(centroid_seq)
centroid_seq <- as.numeric(centroid_seq)

# Subset bin_list_df by the record ids on this list.

bin_list_df_centroid <- subset(bin_list_df_cs, recordID %in% centroid_seq)

class(bin_list_df_centroid$framed)
class(bin_list_df_centroid)

bin_list_df_centroid <- ungroup(bin_list_df_centroid)

bin_list_df_centroid <- bin_list_df_centroid %>% 
  mutate(Sequence_Length = str_count(bin_list_df_centroid$framed))

max(bin_list_df_centroid$Sequence_Length)
min(bin_list_df_centroid$Sequence_Length)

#Outlier detection using distance matrix 

outlier_DNAStringset <- DNAStringSet(bin_list_df_centroid_cs$framed)
names(outlier_DNAStringset) <- bin_list_df_centroid$bin_uri
class(outlier_DNAStringset)

#Running Muscle for the centroid sequences

outlier_align <- DNAStringSet(muscle::muscle(outlier_DNAStringset, log = "logMUSCLE.tx", maxiters = 3, diags = TRUE, verbose = T, gapopen = -3000))

#Transforming the data to DNAbin format

outlier_align <- as.DNAbin(outlier_align)

class(outlier_align)

# I am going to generate the distance matrix with the model TN93.

distanceMatrix_outlier <- dist.dna(x = outlier_align, model = "TN93", as.matrix = TRUE, pairwise.deletion = TRUE)

#All BINs were checked to ensure that their genetic distances did not deviate greatly from the typical range of divergence for the BIN dataset. This was performed by calculating the interquartile range (IQR) of the divergence values of the sequences and identifying sequences that exceeded the upper threshold (quartile 3 + (1.5 × IQR)). These checks ensured the accuracy of the taxonomic assignments and reduced the potential for human error as a result of gross contamination or misidentification (May 2020, SI.) 
#Code adapted from Jacqueline May in github

# Use the upper threshold of the IQR to detect outliers.

lowerQuantile <- quantile(distanceMatrix_outlier)[2]
upperQuantile <- quantile(distanceMatrix_outlier)[4]
iqr <- upperQuantile - lowerQuantile
upperThreshold <- (iqr * 1.5) + upperQuantile
# Remove 0 values so that these are not considered (when a species is compared to itself - the diagonal values).
distanceMatrix_outlier[distanceMatrix_outlier == 0] <- NA
# Convert to datatable.
dfOutliers <- as.data.table(distanceMatrix_outlier, keep.rownames = T)
# Rename the "rn" column (row names).
setnames(dfOutliers, "rn", "species_name")
# Identify BINs with no relatives within "typical" range of genetic divergence (i.e. all of their genetic distances are greater than 1.5 x IQR upper threshold.)
dfOutliers <- dfOutliers[, outlier := apply(.SD, 1, function(x) all(x > upperThreshold, na.rm = T))][outlier == TRUE]

# The results show that there are no outliers on the dataset. 
