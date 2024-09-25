### This script is meant to add predictor data to the parasite data set
# Daki Diaz-Morales (dmdiaz@uw.edu; diazdakeishla@gmail.com)
# September 2024


## The data used in this script on streamflow, temperature, metals (and other elements), nutrients and organic pollutants was downloaded from https://www.waterqualitydata.us.
## The code for the gage sites is USGS-02489500 for the upstream gage (30.79324276, -89.8209072) and USGS-02490193 for the downstream gage (30.70583333, -89.8463889)


### Load necessary packages and datasets-----

library(readxl)
library(readr)
library(dplyr)
library(sjPlot)

## Import parasite dataset

full_dataset <- read_csv("data/processed/Full_dataset_with_psite_life_history_info_2024.09.25.csv")

# Import physical data (which includes water temperature and streamflow)
physicalUSGS <- read_csv("data/Physicochemical/physical_resultphyschem.csv")

# Read USGS-gage site metadata
siteinfo <-  read_csv("data/Physicochemical/station_info.csv")

# Add info on sampling sites (latitude, longitude, CI, etc.) to the physical data

physicalUSGS_withmetadata <- merge(physicalUSGS, siteinfo, by.x = "MonitoringLocationIdentifier", by.y = "MonitoringLocationIdentifier", all.x = TRUE)

### Adding streamflow data----

## Unfortunately streamflow data is only available for the control site. However we assume that the streamflow should be continuous downstream the rivers and do not expect big differences in this parameter.

# Make a subset only for streamflow in m^3/sec
streamflow <- subset(physicalUSGS_withmetadata, Measure=="Stream flow")
streamflow_ms <- subset(streamflow, Unit=="m3/sec")

#### The script below assigns to each fish individual the mean streamflow (mean_streamflow) that the fish experienced over one year before its collection

# Add a new column to store the mean temperatures
full_dataset$mean_streamflow <- NA

library(lubridate)

# Combine year, month, and day columns into a date column
full_dataset$collection_date <- make_date(full_dataset$YearCollected, full_dataset$MonthCollected, 
                                          full_dataset$DayCollected)


streamflow_ms$measurement_date <- make_date(streamflow_ms$YearCollected, streamflow_ms$MonthCollected, 
                                            streamflow_ms$DayCollected)


# Loop through each row of 'full_dataset'
for (i in 1:nrow(full_dataset)) {
  # Extract the collection date for the individual
  collection_date <- as.Date(full_dataset$collection_date[i])
  
  # Define the start and end date for one year before the collection date
  start_date <- collection_date - 365
  end_date <- collection_date
  
  # Filter 'streamflow_ms' to get the streamflow data within that one-year range
  streamflow_data <- streamflow_ms[streamflow_ms$measurement_date >= start_date & streamflow_ms$measurement_date < end_date, "Result"]
  
  # Calculate the mean streamflow and store it in 'full_dataset'
  if (length(streamflow_data) > 0) {
    full_dataset$mean_streamflow[i] <- mean(streamflow_data, na.rm = TRUE)
  } else {
    full_dataset$mean_streamflow[i] <- NA  # If no data is found, store NA
  }
}

# View updated 'full_dataset'
print(full_dataset)

#### The script below assigns to each fish individual the mean streamflow per year regardless of when during that year the fish was collected

# Summarize by taking the mean streamflow per 'year'

# streamflow_mean <- streamflow_ms %>% group_by(YearCollected) %>% 
  # summarise(meanFlow=mean(Result),
            # .groups = 'drop') #For streamflow there is no data available for the impact category but we asume the same flow will be experienced throughout the relatively short strecht of the river we are working with.

# Merge with parasite data

# Full_dataset_streamflow <- merge(full_dataset, streamflow_mean, by.x = "YearCollected", by.y = "YearCollected", all.x = TRUE)


### Adding water temperature data----
# Make a subset of only water temperature
wtemp <- subset(physicalUSGS_withmetadata, Measure=="Temperature, water")

# Check if there is a significant difference between control and impact 

# Summarize by taking the mean, max, and median temperature per 'year'

wtemp_summ <- wtemp %>% group_by(YearCollected,CI) %>% 
  summarise(meanTemp=mean(Result),medianTemp=median(Result),maxTemp=max(Result),
            .groups = 'drop') # For water temperature there is also data for the impact category (although for less years). However, there is no significant difference between control and impact so we take the mean per year ignoring CI. See analysis (LM) relevant to this data above We have to keep in mind that the gage at the impact category is way downstream. So there might be some warming frmo the pulp mill that is probably not captured by the gage.

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

library(lubridate)

# Combine year, month, and day columns into a date column
full_dataset$collection_date <- make_date(full_dataset$YearCollected, full_dataset$MonthCollected, 
                                          full_dataset$DayCollected)


wtemp$measurement_date <- make_date(wtemp$YearCollected, wtemp$MonthCollected, 
                                    wtemp$DayCollected)


# Loop through each row of 'full_dataset'
for (i in 1:nrow(full_dataset)) {
  # Extract the collection date for the individual
  collection_date <- as.Date(full_dataset$collection_date[i])
  
  # Define the start and end date for one year before the collection date
  start_date <- collection_date - 365
  end_date <- collection_date
  
  # Filter 'wtemp' to get the temperature data within that one-year range
  temp_data <- wtemp[wtemp$measurement_date >= start_date & wtemp$measurement_date < end_date, "Result"]
  
  # Calculate the mean temperature and store it in 'full_dataset'
  if (length(temp_data) > 0) {
    full_dataset$mean_temperature[i] <- mean(temp_data, na.rm = TRUE)
  } else {
    full_dataset$mean_temperature[i] <- NA  # If no data is found, store NA
  }
}

# View updated 'full_dataset'
print(full_dataset)

#### The script below assigns to each fish individual the mean temperature per year regardless of when during that year the fish was collected

# Summarize by taking the mean, max, and median temperature per 'year'

wtemp_summ <- wtemp %>% group_by(YearCollected) %>% 
  summarise(meanTemp=mean(Result),medianTemp=median(Result),maxTemp=max(Result),
            .groups = 'drop') # For water temperature there is also data for the impact category (although for less years). However, there is no significant difference between control and impact so we take the mean per year ignoring CI. See analysis (LM) relevant to this data above We have to keep in mind that the gage at the impact category is way downstream. So there might be some warming frmo the pulp mill that is probably not captured by the gage.

# Merge with parasite data

Full_dataset_physical <- merge(Full_dataset_streamflow, wtemp_summ, by.x = "YearCollected", by.y = "YearCollected", all.x = TRUE)

### Season----
## Create a season column to be later used as a random factor in the models

full_dataset$MonthCollected <- as.numeric(full_dataset$MonthCollected)

full_dataset <- full_dataset %>% 
  mutate(season = case_when(  MonthCollected >= 1 &
                                MonthCollected <= 2 
                              ~ "winter",
                              MonthCollected >= 3 &
                                MonthCollected <= 5 ~ "spring",
                              MonthCollected >= 6 &
                                MonthCollected <= 8 ~ "summer",
                              MonthCollected >= 9 &
                                MonthCollected <= 11 ~ "fall",
                              MonthCollected == 12 ~ "winter",
  ))

full_dataset$season <- as.factor(full_dataset$season)

#### Export the final sheet----

write.csv(full_dataset, file="data/processed/Full_dataset_physical_2024.09.25.csv")


