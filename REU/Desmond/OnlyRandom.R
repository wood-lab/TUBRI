#Attempting to create body size graphs including ONLY randomly-selected fishes

setwd("C:/Users/Test 1/OneDrive/Documents/WOODLAB/TUBRI")

#reading in ictalurus data

#processed data
library(readr)
ictpun <- read_csv("data/processed/Ictalurus_punctatus_processed_machine_readable_UPDATED_2024.07.10.csv")
View(ictpun)

#rawdata
View(ictpunraw)

#plotting ALL
plot(ictpunraw$StandardLength_mm~ictpunraw$YearCollected)

