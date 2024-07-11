#Attempting to create body size graphs including ONLY randomly-selected fishes

setwd("C:/Users/Test 1/OneDrive/Documents/WOODLAB/TUBRI")

#ICTALURUS PUNCTATUS

#processed data
library(readr)
ictpun <- read_csv("data/processed/Ictalurus_punctatus_processed_machine_readable_UPDATED_2024.07.10.csv")
View(ictpun)

#rawdata
View(ictpunraw)

#plotting ALL
plot(ictpunraw$StandardLength_mm~ictpunraw$YearCollected)
#BOO! no pattern.


#which rows contain fish where selection was influenced by size?
#13, 17, 35, 39, 44, 46, 56, 57, 65, 67, 69, 70, 74, 79, 80, 87
ictpuntrim <- ictpunraw[-c(13, 17, 35, 39, 44, 46, 56, 57, 65, 67, 69, 70, 74, 79, 80, 87), ]
View(ictpuntrim)

#plot only randomly selected fishes
plot(ictpuntrim$StandardLength_mm~ictpuntrim$YearCollected)
summary(lm(ictpuntrim$StandardLength_mm~ictpuntrim$YearCollected))
abline(lm(ictpuntrim$StandardLength_mm~ictpuntrim$YearCollected))

#NOTROPIS ATHERNOIDES

#plotting ALL
plot(notropraw$StandardLength_mm~notropraw$YearCollected)

#which rows contain size-selected fish?
#