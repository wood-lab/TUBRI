#Script for Imani to select random monogeneans from Notropis atherinoides
#Revised July 12, 2024
#Authors: Imani Jones, Dakeishla Diaz, Chelsea Wood


NOTATH <- read_csv("data/processed/Notropis_atherinoides_processed_machine_readable.csv")


#Make a subset of only monogeneans

NOTATH_mon <- subset(NOTATH, psite_spp) 


# importing required libraries 
library("dplyr") 

# pick 3 samples of each from data frame 
NOTATH_select <- NOTATH %>% group_by(col1) %>% sample_n(3) 
print("Modified DataFrame") 
print (data_mod)