## This script adds nutrients and element concentrations to the myxozoan data
## Created by Daki Diaz Morales (diazdakeishla@gmail.com)


#### Load packages-----
library(readr)
library(dplyr)
library(lubridate)
library(ggplot2)

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

# Make a subset of only water temperature
wtemp <- subset(physicalUSGS_withmetadata, Measure=="Temperature, water")

## Check if there is a significant difference between control and impact 

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
full_dataset_with_LH$mean_temperature <- NA

# Combine year, month, and day columns into a date column
full_dataset_with_LH$collection_date <- make_date(full_dataset_with_LH$YearCollected, full_dataset_with_LH$MonthCollected, 
                                                  full_dataset_with_LH$DayCollected)

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

#### Summary stats-----
## Remove columns with myxozoans that aren't confirmed to be myxozoans and therefore are not identified
full_dataset_myxo <- subset(full_dataset_myxo, !is.na(Parasite_genus)) # These were thought to be myxozoans at the beginning but they are not, therefore no genus was assigned
full_dataset_myxo <- subset(full_dataset_myxo, psite_spp.x!="MYX.GM") # According to Stephen Atkinson (OSU), this is likely a contamination and also it is in the gall bladder so does not count for abundance
full_dataset_myxo <- subset(full_dataset_myxo, psite_spp.x!="MYX.GEOM") # According to Stephen Atkinson (OSU), this is not a myxozoan

## Summary stats
prevalence_myxo <- subset(full_dataset_myxo, !is.na(psite_presence))
abundance_myxo <- subset(full_dataset_myxo, !is.na(psite_count))

## Which parasites have a prevalence of more than 5%
summary_prevalence <- prevalence_myxo %>% group_by(Fish_sp.x,psite_spp.x) %>% 
  summarise(Prevalence=mean(psite_presence),
            .groups = 'drop') 

prevalence_higherthan5percent <- subset(summary_prevalence, Prevalence > 0.05)

# The parasites with more than 5% of prevalence are: MYX.CM, MYX.F, MYX.G, MYX.SWTY, MYX.TAIL, MYX.SP, MYX.GO, MYX.THEL, MYXO.SBAD

### FIGURE 1 ###

## Make plot showing the prevalence of infection of fish per parasite
# What is the number of fish we dissected per fish species
fish_dissected <- full_dataset_myxo %>% group_by(Fish_sp.x) %>% 
  summarise(number_dissected=length(unique(IndividualFishID)),
            .groups = 'drop') 

# Calculate prevalence by parasite genus and fish species
summary_psitegenus <- prevalence_myxo %>%
  group_by(Parasite_genus, Fish_sp.x,IndividualFishID) %>%            
  summarize(genus_numinfected = sum(psite_presence), .groups = "drop") # This part ensures that we have one fish per parasite genus

summary_psitegenus$psite_presence <- ifelse(summary_psitegenus$genus_numinfected > 0,# create binary data to calculate prevalence
                                            1,    # what if condition is TRUE
                                            0)       # what if condition is FALSE

# Calculate prevalence by parasite genus and fish species
summary_prevalence <- summary_psitegenus %>%
  group_by(Parasite_genus, Fish_sp.x) %>%             # Group by parasite genus and snail species
  summarize(Prevalence = (mean(psite_presence)*100), .groups = "drop")

# plot
ggplot(summary_prevalence, aes(x= Fish_sp.x,
                               y=Prevalence,color=Parasite_genus,
                               fill=Parasite_genus))+
  geom_bar(stat = "identity")+apatheme+
  ggtitle(NULL)+
  scale_color_manual(name = "Parasite genus:",
                     values = c("Chloromyxum" = "#543005", 
                                "Henneguya" = "#8c510a", 
                                "Myxidium" = "#dfc27d", 
                                "Myxobolus" = "#c7eae5",
                                "Thelohanellus" = "#35978f",
                                "Unicauda" = "#01665e"),
                     labels = c("Chloromyxum" = expression(italic("Chloromyxum")*" sp."), 
                                "Henneguya"= expression(italic("Henneguya")*" sp."), 
                                "Myxidium" = expression(italic("Myxidium")*" sp."), 
                                "Myxobolus" = expression(italic("Myxobolus")*" spp."),
                                "Thelohanellus" = expression(italic("Thelohanellus")*" sp."),
                                "Unicauda" = expression(italic("Unicauda")*" sp."))) +
  scale_fill_manual(name = "Parasite genus:",
                    values = c("Chloromyxum" = "#543005", 
                               "Henneguya" = "#8c510a", 
                               "Myxidium" = "#dfc27d", 
                               "Myxobolus" = "#c7eae5",
                               "Thelohanellus" = "#35978f",
                               "Unicauda" = "#01665e"), 
                    labels = c("Chloromyxum" = expression(italic("Chloromyxum")*" sp."), 
                               "Henneguya"= expression(italic("Henneguya")*" sp."), 
                               "Myxidium" = expression(italic("Myxidium")*" sp."), 
                               "Myxobolus" = expression(italic("Myxobolus")*" spp."),
                               "Thelohanellus" = expression(italic("Thelohanellus")*" sp."),
                               "Unicauda" = expression(italic("Unicauda")*" sp."))) +
  
  xlab("Fish species")+ylab("Prevalence of infection (%)")+
  scale_x_discrete(limits = c("Carpiodes velifer", 
                              "Pimephales vigilax", 
                              "Notropis atherinoides",
                              "Gambusia affinis",
                              "Hybognathus nuchalis",
                              "Ictalurus punctatus"),
                   labels = c("Carpiodes velifer"="Carpiodes velifer 
                             n = 180", 
                              "Pimephales vigilax"="Pimephales vigilax 
                             n = 193", 
                              "Notropis atherinoides"="Notropis atherinoides 
                             n = 206",
                              "Gambusia affinis"="Gambusia affinis 
                             n = 208",
                              "Hybognathus nuchalis"="Hybognathus nuchalis 
                             n = 221",
                              "Ictalurus punctatus"="Ictalurus punctatus 
                             n = 88")) + ylim(0,100)

ggsave(file="/Users/dakeishladiaz-morales/TUBRI/Manuscripts/Myxozoans/Figures/Figure1.png", width=240, height=180, dpi=1000, units = "mm")
ggsave(file="/Users/dakeishladiaz-morales/TUBRI/Manuscripts/Myxozoans/Figures/Figure1.pdf", width=240, height=180, dpi=1000, units = "mm")
