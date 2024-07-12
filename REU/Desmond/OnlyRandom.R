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
#35, 36, 37, 54, 75, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 96, 99, 100, 101, 102, 106, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 121, 123, 125, 126, 127, 128, 129, 130, 131, 132, 133, 136, 137, 138, 139, 141, 142, 144, 145, 146, 147, 148, 149, 151, 152, 154, 155, 156, 158, 159, 160, 161, 162, 164, 165, 168, 169, 171, 172, 173, 174, 175, 176, 177, 178, 179, 181, 182, 183, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208
notroptrim <- notropraw[-c(35, 36, 37, 54, 75, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 96, 99, 100, 101, 102, 106, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 121, 123, 125, 126, 127, 128, 129, 130, 131, 132, 133, 136, 137, 138, 139, 141, 142, 144, 145, 146, 147, 148, 149, 151, 152, 154, 155, 156, 158, 159, 160, 161, 162, 164, 165, 168, 169, 171, 172, 173, 174, 175, 176, 177, 178, 179, 181, 182, 183, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208), ]
View(notroptrim)

#plot only randomly selected fishes
plot(notroptrim$StandardLength_mm~notroptrim$YearCollected)
summary(lm(notroptrim$StandardLength_mm~notroptrim$YearCollected))
abline(lm(notroptrim$StandardLength_mm~notroptrim$YearCollected))
#WOO!