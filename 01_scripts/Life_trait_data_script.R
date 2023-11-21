#We gathered information on life history traits to examine the ecological component of the anomalous diversity pattern in parasitic wasps. Idiobiosis (0) versus koinobiosis (1) classification at the genus level was sourced from recent phylogenomic studies of Ichneumonidae 66 and Braconidae 67 where this data was used for ancestral state reconstruction. Additionally, we used a data set comprising information for 382 Ichneumonoidea parasitoid species (26 braconid and 25 ichneumonid subfamilies) (Traynor, 2004). Following 67, the biological traits were extrapolated from congeneric species using at least one member of the genus and considering no conflicts among datasets. This considers that host use is almost always conserved within genera, which reduces the possibility of finding genera with records of both koinobionts and idiobionts 67,68. The final dataset comprised idiobioisis and koinobiosis classification at the genus level for 1492 molecular clusters. 

#Life trait data summary for Ichenumonoidea. The files will be uploaded in order to reproduce the code. 

#Reading data Braconidae 
phylog_data <- read.csv("phylo_braco.csv")

#Filtering step
phylog_data <- phylog_data[-(369:421), ]

#Braconbidae and Ichneumonidae, thesis Traynor, 2004.  
traynor_data <- read.csv ("traynor_data.csv")
traynor_data <- traynor_data[, c(1, 5, 6)]

#Selecting just Genus name from the species name clasification. 
traynor_data$Genus <- word(traynor_data$Taxonomy, 1)

#Filtering records missing information. 
traynor_data <- traynor_data %>% 
  filter(Genus != "Family" & Genus != "Tribe" & Genus != "Subfamily")
traynor_data <- traynor_data[-(395:1919), ]

#Joining datasets
life_trait_data <- full_join(traynor_data, phylog_data, by = "Genus")
life_trait_data <- life_trait_data[, c(3, 4, 10)]

#Changing class of the K-I classification
life_trait_data$Idiobiont..0..Koinobiont..1.<- as.numeric(life_trait_data$Idiobiont..0..Koinobiont..1.)
life_trait_data$Idiobiont..0..Koinobiont..1.<- life_trait_data$Idiobiont..0..Koinobiont..1.+1

#Reading Ichneumonidae data 
phylo_ichne <- read.csv("phylo_ichne.csv")

#String manipulation to keep genus name 
phylo_ichne$Taxon <- str_replace(phylo_ichne$Taxon, pattern = "^*_", replacement = "")
phylo_ichne$Taxon <- sub(".*?_", "", phylo_ichne$Taxon)
phylo_ichne$Taxon <- sub("_sp.*", "", phylo_ichne$Taxon)
phylo_ichne$Taxon <- sub("_.*", "", phylo_ichne$Taxon)

#Changing class of K-I calssification
phylo_ichne$Character.1 <- as.numeric(phylo_ichne$Character.1)

#Removing records withouth clasification. 
phylo_ichne <- phylo_ichne %>%
  filter(!is.na(Character.1))

#Removing duplicates
phylo_ichne <- phylo_ichne[!duplicated(phylo_ichne[ , c("Taxon","Character.1")]),]

#String manipulation to keep genus information equal among datasets
life_trait_data$Genus <- sub("^cf.", "", life_trait_data$Genus)
life_trait_data$Genus <- sub("\\s*\\([^\\)]+\\)", "", life_trait_data$Genus)
life_trait_data$Genus <- sub("\\[.*?\\]", "", life_trait_data$Genus)

life_trait_data <- life_trait_data[!with(life_trait_data,is.na(X5)& is.na(Idiobiont..0..Koinobiont..1.)),]

write_csv(life_trait_data, "life_trait_data.csv")

#Changing column names to merge
names(phylo_ichne) <- c("Genus", "Character")

#Standardizing classiifcation to  `1` = "idiobiont",`2` = "koinobiont"
phylo_ichne$Character <- phylo_ichne$Character + 1


all_data_LT <- full_join(life_trait_data, phylo_ichne, by = "Genus")

all_data_LT <- distinct(all_data_LT)

all_data_LT <- all_data_LT %>% 
  mutate(across(Character.x, ~ ifelse(is.na(Character.x), Character.y, Character.x)))

all_data_LT <- all_data_LT[, c(1, 2)]
names(all_data_LT) <- c("Genus", "Character")

##Getting coordinates for the BINS
data_df_GD_LT <- data_df_GD[, c(2,15, 19, 20)]

colnames(data_df_GD_LT ) <- c("bin_uri", "Genus", "lat", "long")

#Here I use values of nucleotide diversity per BIN estimated at the latitudinal band level to explore relationship of the values with the K-I classification of the genus. 

Nuc_div_df_trait <- bind_rows(Nuc_div_df_lat)
Nuc_div_df_trait$bin_uri <- sub(".fa", "", Nuc_div_df_trait$bin_uri)

data_genus <- data_df_GD[, c(2, 15)]

data_genus <- drop_na(data_genus)
names(data_genus) <- c("bin_uri", "Genus")

trait_data <- left_join(Nuc_div_df_trait, data_genus, by = "bin_uri")

trait_data <- drop_na(trait_data)
trait_data <- distinct(trait_data)

trait_data2 <- left_join(all_data_LT, trait_data, by = "Genus")
trait_data2 <- drop_na(trait_data2)

write_csv(trait_data2, "trait_data_all.csv")

trait_data_coord <- left_join(trait_data2, data_df_GD_LT, by = c("Genus", "bin_uri"))
trait_data_coord <- distinct(trait_data_coord)

#Graphic of BINs Genus K-I classification 
trait_data_coord <- trait_data_coord %>% group_by(bin_uri, Pi_per_site, cell, Character) %>% sample_n(size = 1)
trait_data_coord$Character <- as.factor(trait_data_coord$Character)

#Linear Model of latitude and nucleotide diversity with life trait effect 

b <- ggplot(trait_data_coord, aes(x = abs(lat), y = sqrt(Pi_per_site)))
# Scatter plot with regression line
b + geom_point()+
  geom_smooth(aes(color = Character), 
              method = "lm", fullrange = TRUE) +
  facet_wrap(~Character ~ ., labeller = as_labeller(Character_names)) +
  labs( x = "Absolute value of Latitude", y = "Nucleotide Diversity per Molecular Cluster", color = "Life Trait") +
  scale_color_manual(labels = c("idiobiont", "koinobiont"), values = c("#C7E020FF", "#FDE725FF")) +
  #scale_color_manual(values = c("#440154FF", "#FDE725FF"), name = "Life Trait")+
  scale_fill_manual(values = c("#C7E020FF", "#FDE725FF"))+
  ggpubr::stat_cor(label.x = 3)+
  ylab("Nucleotide Diversity per BIN")+
  xlab("Absolute value of Latitude")+
  theme_bw()+
  theme(axis.text = element_text(size = 14)) +
  theme(axis.title = element_text(size = 20)) +
  theme(strip.text = element_text(size = 14))

trait_data_coord_K <- trait_data_coord %>%
  filter(Character == 2)
mod_LT <- lm(sqrt(trait_data_coord_K$Pi_per_site) ~ abs(trait_data_coord_K$lat))

summary(mod_LT)

plot(mod_LT)

Character_names <- c(
  `1` = "idiobiont",
  `2` = "koinobiont"
)


##trying a Sankey diagram for a representation of the location of the records. 
install.packages("riverplot")
library(riverplot)

trait_data2$Lat <- sapply(trait_data2$cell, str_extract, "-?[0-9.]+")

trait_data2 <- dplyr::filter(trait_data2, !grepl("Malai",Genus))

unique(trait_data2$Lat )

node_names <- c("70", "60", "50", "40", "30", "20", "10", "0", "-10", "-20", "-30", "-40", "-50", "1", "2")
nodes <- data.frame(name = node_names)
node_x <- c(0.1, 0.1 , 0.1, 0.1 , 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.99, 0.99)
node_y <- c(0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 0.1, 0.70)

links <- data.frame(source = match(trait_data2$Lat, node_names) - 1, target = match(trait_data2$Character, node_names) - 1, value = trait_data2$Pi_per_site)

#Plot
ps <- plot_ly(type='sankey',
        orientation = "h",
        arrangement = "fixed",
        
        node = list(
          label = node_names,
          x = node_x, 
          y = node_y, 
          color = c("#440154FF", "#481B6DFF", "#46337EFF", "#3F4889FF","#365C8DFF", "#2E6E8EFF", "#277F8EFF", "#21908CFF", "#1FA187FF", "#2DB27DFF", "#4AC16DFF", "#71CF57FF", "#9FDA3AFF", "#CFE11CFF", "#FDE725FF"),
          pad = 20,
          thinkness = 0.5,
          line = list(width = 0.05)
          ),
        
        link = list(
          source = links$source,
          target = links$target,
          value = links$value))

layout <- list(
  font = list(size = 15), 
  title = "Trait data"
)
ps <- layout(ps, font=layout$font, title=layout$title)
ps



