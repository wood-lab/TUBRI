### This script is meant to analyse the impact of multiple stressors on myxozoa abundance
### Written by Daki Diaz-Morales
### Revised August 17, 2024



### Preparation of data frames----

## Import data

full_dataset <- read_csv("data/processed/Full_dataset_with_psite_life_history_info_2024.08.17.csv")

## Make a subset for myxozoans

full_dataset_myxo <- subset(full_dataset, Parasite_taxonomic_group == "Myxozoa")

## Make Parasite_genus a factor

full_dataset_myxo$Parasite_genus <- as.factor(full_dataset_myxo$Parasite_genus)



## Make a subset for data before the clean water act

full_dataset_myxo_bca <- subset(full_dataset_myxo, YearCollected < 1973)



### Before clean water act----

##What is the impact of distance from pulp mill on the abundance of myxozoans?

#Set a theme for all your plots
apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

# Plot myxozoan abundance against distance from pulp mill (latitude)

ggplot(full_dataset_myxo_bca,aes(Latitude,psite_count, color = Parasite_genus,group = Parasite_genus))+
  geom_point(size=4)+
  facet_wrap("Fish_sp.x")+
  geom_vline(xintercept=30.76, linetype="dashed", color = "black", size=0.5)+
  apatheme


