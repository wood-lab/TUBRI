## This script adds nutrients and element concentrations to the myxozoan data
## Created by Daki Diaz Morales (diazdakeishla@gmail.com)


#### Load packages-----
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyverse)
library(FactoMineR)
library(missMDA)

#### Load data and data exploration-----
## Set a theme for downstream plots
apatheme= theme_bw(base_size = 18,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## Import parasite dataset
full_dataset <- read_csv("data/processed/Full_dataset_with_psite_life_history_info_2024.11.21.csv")

# What is the number of fish we dissected per fish species
fish_dissected <- full_dataset %>% group_by(Fish_sp.x) %>% 
  summarise(number_dissected=length(unique(IndividualFishID)),
            .groups = 'drop') 

# What is the size range of fish we dissected per fish species
fish_oneperfish <- full_dataset %>% group_by(Fish_sp.x,IndividualFishID) %>% 
  summarise(TL_mm=mean(TotalLength_mm),
            .groups = 'drop') # Remove multiple entries of a single fish individual

fish_sizerange <- fish_oneperfish %>% group_by(Fish_sp.x) %>% 
  summarise(mean_TL=mean(TL_mm),
            min_TL=min(TL_mm),
            max_TL=max(TL_mm),
            .groups = 'drop') # Summarize data on size

## Create column for presence and absence of myxozoans, useful for gallbladder-infecting myxozoans
full_dataset$psite_presence <- ifelse(full_dataset$psite_count > 0,               
                                      1,    # what if condition is TRUE
                                      0)       # what if condition is FALSE

#### Add temperature data----

## This was done in the mother data frame before importing the mother data frame here but the code is added for reference

# Import physical data (which includes water temperature and streamflow)
physicalUSGS <- read_csv("data/Physicochemical/physical_resultphyschem.csv")

# Read USGS-gage site metadata
siteinfo <-  read_csv("data/Physicochemical/station_info.csv")

# Add info on sampling sites (latitude, longitude, CI, etc.) to the physical data

physicalUSGS_withmetadata <- merge(physicalUSGS, siteinfo, by.x = "MonitoringLocationIdentifier", by.y = "MonitoringLocationIdentifier", all.x = TRUE)

# Make a subset of only water temperature
wtemp <- subset(physicalUSGS_withmetadata, Measure=="Temperature, water")

## Check if there is a significant difference between control and impact 

# Summarize by taking the mean, max, and median temperature per 'year'

wtemp_summ <- wtemp %>% group_by(YearCollected,CI) %>% 
  summarise(meanTemp=mean(Result),medianTemp=median(Result),maxTemp=max(Result),
            .groups = 'drop') # For water temperature there is also data for the impact category (although for less years). However, there is no significant difference between control and impact so we take the mean per year ignoring CI. See analysis (LM) relevant to this data above We have to keep in mind that the gage at the impact category is way downstream. So there might be some warming frmo the pulp mill that is probably not captured by the gage.


wtemp_month <- wtemp %>% group_by(MonthCollected,CI) %>% 
  summarise(meanTemp=mean(Result),medianTemp=median(Result),maxTemp=max(Result),
            .groups = 'drop') # For water temperature there is also data for the impact category (although for less years). However, there is no significant difference between control and impact so we take the mean per year ignoring CI. See analysis (LM) relevant to this data above We have to keep in mind that the gage at the impact category is way downstream. So there might be some warming frmo the pulp mill that is probably not captured by the gage.

# Data overview
mean(wtemp_summ$meanTemp)
max(wtemp_summ$maxTemp)



# wtemp change throughout time

wtemp_summ$CI <- as.factor(wtemp_summ$CI)

m1 <- lm(meanTemp~YearCollected*CI, data=wtemp_summ)
m2 <- lm(meanTemp~YearCollected, data=wtemp_summ)

summary(m1) 
summary(m2) 

tab_model(m1)
tab_model(m2)# The estimate of YearCollected is the same/consistent for both models (an increase of 0.09 degC/year). No difference between control and impact.

#### The script below assigns to each fish individual the mean temperature (mean_temperature) that the fish experienced over one year before its collection

# Add a new column to store the mean temperatures
full_dataset$mean_temperature <- NA

# Combine year, month, and day columns into a date column
full_dataset$collection_date <- make_date(full_dataset$YearCollected, full_dataset_with_LH$MonthCollected, 
                                          full_dataset$DayCollected)

wtemp$measurement_date <- make_date(wtemp$YearCollected, wtemp$MonthCollected, 
                                    wtemp$DayCollected)

# Loop through each row of 'full_dataset_with_LH'
for (i in 1:nrow(full_dataset_with_LH)) {
  # Extract the collection date for the individual
  collection_date <- as.Date(full_dataset_with_LH$collection_date[i])
  
  # Define the start and end date for one year before the collection date
  start_date <- collection_date - 365
  end_date <- collection_date
  
  # Filter 'wtemp' to get the temperature data within that one-year range
  temp_data <- wtemp[wtemp$measurement_date >= start_date & wtemp$measurement_date < end_date, "Result"]
  
  # Calculate the mean temperature and store it in 'full_dataset_with_LH'
  if (length(temp_data) > 0) {
    full_dataset_with_LH$mean_temperature[i] <- mean(temp_data, na.rm = TRUE)
  } else {
    full_dataset_with_LH$mean_temperature[i] <- NA  # If no data is found, store NA
  }
}

#### Add nutrients----
## Import data for nutrients
nutrients <- read_csv("data/Physicochemical/nutrients_resultphyschem.csv")

# Read USGS-gauge site metadata
siteinfo <-  read_csv("data/Physicochemical/station_info.csv")

# Add info on sampling sites (latitude, longitude, CI, etc.)
nutrients_withmetadata <- merge(nutrients, siteinfo, by.x = "MonitoringLocationIdentifier", by.y = "MonitoringLocationIdentifier", all.x = TRUE)

## Nitrogen species: nitrogen available to plants occurs in the form of 
# ammonia (NH3), ammonium (NH4+), nitrate (NO3−), nitrite (NO2-), and 
# organic nitrogen. Therefore, we used the mixed nitrogen measurement for our analyses.

# Take only mixed forms of nitrogen in mg-N/L
mixN_unitall <- subset(nutrients_withmetadata, Measure == "Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3)")

## The script below assigns to each fish individual the mean nutrients (mean_nitrogen) that the fish experienced over one year before its collection

# Add a new column to store the mean nutrient in the form of mixed nitrogen forms
full_dataset$mean_nitrogen <- NA

mix_N_mgL <- subset(mixN_unitall, Unit == "mg/l")

# Combine year, month, and day columns into a date column
full_dataset$collection_date <- make_date(full_dataset$YearCollected, full_dataset$MonthCollected, 
                                          full_dataset$DayCollected)

mix_N_mgL$measurement_date <- make_date(mix_N_mgL$YearCollected, mix_N_mgL$MonthCollected, 
                                        mix_N_mgL$DayCollected)

# Loop through each row of 'full_dataset'
for (i in 1:nrow(full_dataset)) {
  # Extract the collection date for the individual
  collection_date <- as.Date(full_dataset$collection_date[i])
  
  # Define the start and end date for one year before the collection date
  start_date <- collection_date - 365
  end_date <- collection_date
  
  # Filter 'mix_N_mgL' to get the nitrogen data within that one-year range
  nutrient_data <- mix_N_mgL[mix_N_mgL$measurement_date >= start_date & mix_N_mgL$measurement_date < end_date, "Result"]
  
  # Calculate the mean nitrogen and store it in 'full_dataset'
  if (length(nutrient_data) > 0) {
    full_dataset$mean_nitrogen[i] <- mean(nutrient_data, na.rm = TRUE)
  } else {
    full_dataset$mean_nitrogen[i] <- NA  # If no data is found, store NA
  }
}

#### Create daughter data frames----
## Take only myxozoan parasites
full_dataset_myxo <- subset(full_dataset, Parasite_taxonomic_group == "Myxozoa")

## Take the log of the size to later include it in the analysis as offset
full_dataset_myxo$logTL_mm <- log(full_dataset_myxo$TotalLength_mm)

## Remove columns with myxozoans that aren't confirmed to be myxozoans and therefore are not identified
full_dataset_myxo <- subset(full_dataset_myxo, !is.na(Parasite_genus)) # These were thought to be myxozoans at the beginning but they are not, therefore no genus was assigned
full_dataset_myxo <- subset(full_dataset_myxo, psite_spp.x!="MYX.GM") # According to Stephen Atkinson (OSU), this is likely a contamination and also it is in the gall bladder so does not count for abundance
full_dataset_myxo <- subset(full_dataset_myxo, psite_spp.x!="MYX.GEOM") # According to Stephen Atkinson (OSU), this is not a myxozoan

## Export the dataset for downstream analyses
write.csv(full_dataset_myxo, 
          file = "/Users/dakeishladiaz-morales/TUBRI/Manuscripts/Myxozoans//data/myxozoans_time.csv", 
          row.names = FALSE)

#### Add element concentrations----

## We add the element concentrations at this stage because we have to reduce 
## dimensionality with a Principal Component Analysis and I want to restrict it 
## to the data between 1972-1995, which is the data we will use for this analysis.

## Crease subset for myxozoans that can be counted (pseudo-abundance).
myxo_count <- subset(full_dataset_myxo, Parasite_genus == "Myxobolus"|Parasite_genus == "Henneguya"|Parasite_genus == "Unicauda"|Parasite_genus == "Thelohanellus")

## Restrict to year and control for the interactions analysis

# Restrict to the years that we are going to work with
myxo_count_y <- subset(myxo_count, 
                       YearCollected > 1972 &
                         YearCollected < 1995 )

# Restrict to control because we have environmental data for that area of the river
myxo_count_yc <- subset(myxo_count_y, 
                        CI == "control" )

### Elements
## Import data for metals
inorganics <- read_csv("data/Physicochemical/inorganics_resultphyschem.csv")

# Add info on sampling sites (latitude, longitude, CI, etc.)
inorganics_withmetadata <- merge(inorganics, siteinfo, by.x = "MonitoringLocationIdentifier", by.y = "MonitoringLocationIdentifier", all.x = TRUE)

# Summarize data table to year, day, and month collected
inorganics_means_ymd <- inorganics_withmetadata %>% group_by(YearCollected,MonthCollected,DayCollected,CI,Unit,Measure) %>% 
  summarise(Result=mean(Result),
            .groups = 'drop') 

# make the table wide
elements_na_wide <- spread(inorganics_means_ymd, Measure,Result)

## take only the years and control sites
# Restrict to the years that we are going to work with
elements_y <- subset(elements_na_wide, 
                       YearCollected > 1972 &
                         YearCollected < 1995 )

# Restrict to control because we have environmental data for that area of the river
elements_yc <- subset(elements_y, 
                        CI == "control" )

# Make a subset for element concentration in water

elements_water_w <- subset(elements_yc, Unit=="ug/l"|
                             Unit=="mg/l")


# Combine year, month, and day columns into a date column
elements_water_w$measurement_date <- make_date(elements_water_w$YearCollected, 
                                               elements_water_w$MonthCollected, 
                                               elements_water_w$DayCollected)

# Initialize new columns in myxo_count to store mean values for each variable

myxo_count_yc$As <- NA
myxo_count_yc$Ba <- NA
myxo_count_yc$Cd <- NA
myxo_count_yc$Cr <- NA
myxo_count_yc$Cu <- NA
myxo_count_yc$Pb <- NA
myxo_count_yc$Fe <- NA
myxo_count_yc$Mg <- NA
myxo_count_yc$Mn <- NA
myxo_count_yc$Hg <- NA
myxo_count_yc$Ni <- NA
myxo_count_yc$Zn <- NA

# Loop through each row of 'myxo_count_yc'
for (i in 1:nrow(myxo_count_yc)) {
  # Extract the collection date for the individual
  collection_date <- as.Date(myxo_count_yc$collection_date[i])
  
  # Define the start and end date for one year before the collection date
  start_date <- collection_date - 365
  end_date <- collection_date
  
  # Filter "elements_water_w" to get the data within that one-year range
  filtered_data <- elements_water_w[elements_water_w$measurement_date >= start_date & elements_water_w$measurement_date < end_date, ]
  
  # Calculate the mean for each variable and store it in 'myxo_count_yc'
  if (nrow(filtered_data) > 0) {
    myxo_count_yc$As[i] <- mean(filtered_data$Arsenic, na.rm = TRUE)
    myxo_count_yc$Ba[i] <- mean(filtered_data$Barium, na.rm = TRUE)
    myxo_count_yc$Cd[i] <- mean(filtered_data$Cadmium, na.rm = TRUE)
    myxo_count_yc$Cr[i] <- mean(filtered_data$Chromium, na.rm = TRUE)
    myxo_count_yc$Cu[i] <- mean(filtered_data$Copper, na.rm = TRUE)
    myxo_count_yc$Pb[i] <- mean(filtered_data$Lead, na.rm = TRUE)
    myxo_count_yc$Fe[i] <- mean(filtered_data$Iron, na.rm = TRUE)
    myxo_count_yc$Mg[i] <- mean(filtered_data$Magnesium, na.rm = TRUE)
    myxo_count_yc$Mn[i] <- mean(filtered_data$Manganese, na.rm = TRUE)
    myxo_count_yc$Hg[i] <- mean(filtered_data$Mercury, na.rm = TRUE)
    myxo_count_yc$Ni[i] <- mean(filtered_data$Nickel, na.rm = TRUE)
    myxo_count_yc$Zn[i] <- mean(filtered_data$Zinc, na.rm = TRUE)
    
  } else {
    myxo_count_yc$As[i] <- NA
    myxo_count_yc$Ba[i] <- NA
    myxo_count_yc$Cd[i] <- NA
    myxo_count_yc$Cr[i] <- NA
    myxo_count_yc$Cu[i] <- NA
    myxo_count_yc$Pb[i] <- NA
    myxo_count_yc$Fe[i] <- NA
    myxo_count_yc$Mg[i] <- NA
    myxo_count_yc$Mn[i] <- NA
    myxo_count_yc$Hg[i] <- NA
    myxo_count_yc$Ni[i] <- NA
    myxo_count_yc$Zn[i] <- NA
  }
}

## Remove non-numeric variables 

elements_na_wide_num <- myxo_count_yc[-c(1:39)]
View(elements_na_wide_num)

# scale variables
elements_na_wide_scaled <- scale(elements_na_wide_num)

# Estimate the number of dimensions for the Principal Component Analysis by cross-validation
# The number of components which leads to the smallest mean square error of prediction (MSEP) is retained. 
# For the Kfold cross-validation, pNA percentage of missing values is inserted and predicted with a PCA model using ncp.min to ncp.max dimensions.
nb <- estim_ncpPCA(elements_na_wide_scaled,method.cv = "Kfold", verbose = FALSE)
nb$ncp #5

plot(0:5, nb$criterion, xlab = "nb dim", ylab = "MSEP")

# Because we have a lot of NAs and PCAs can't handle that, we have to impute values based on our data
# The (regularized) iterative PCA algorithm (by default) first consists imputing missing values with initial values such as the mean of the variable.
# It is adviced to use the regularized version of the algorithm to avoid the overfitting problems which are very frequent when there are many missing values. 
res.comp <- imputePCA(elements_na_wide_scaled, ncp = 5) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set

imp <- cbind.data.frame(res.comp$completeObs)

res.pca <- PCA(imp, quanti.sup = 1, quali.sup = 12, ncp = 5, graph=FALSE) # generate PCA

plot(res.pca, choix="var")+apatheme# visualize PCA

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Results/elements_pca/PCA_elements.pdf",width=8,height=5,colormodel="rgb")
plot(res.pca, choix="var")+apatheme# visualize PCA
dev.off()

# Extract contribution of variables to the components
contributions_pca <- res.pca$var$contrib
loadings_pca <- res.pca$var$coord

# Save contributions
capture.output(contributions_pca, file = "Manuscripts/Myxozoans/Results/elements_pca/contributions.txt")

# Add principal components to your data frame
score <- as_tibble(factoextra::get_pca_ind(res.pca)$coord) #extract individual scores to be used in glm

myxo_count_ycm <- cbind(myxo_count_yc, score[1:2]) # merge scores with original data

# Rename columns
myxo_count_ycm <- myxo_count_ycm %>%
  rename(
    Elements_PC1 = Dim.1,
    Elements_PC2 = Dim.2
  )

## Export the dataset for downstream analyses
write.csv(myxo_count_ycm, 
          file = "/Users/dakeishladiaz-morales/TUBRI/Manuscripts/Myxozoans//data/myxozoans_multiplestressors.csv", 
          row.names = FALSE)

