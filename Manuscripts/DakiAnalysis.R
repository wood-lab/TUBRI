
#### Load packages----

##Load packages

library(readxl)
library(stats)
library(DHARMa)
library(MASS)
library(sjPlot)
library(ggplot2)
library(ggeffects)
library(sjmisc)
library(sjlabelled)
library(lme4)
library(investr)
library(ggsci)
library(ggpubr)
library(multcomp)
library(effects)
library(glmmTMB)
library(itsadug)
library(mgcv)
library(readr)
library(dplyr)
library(MuMIn)
library(tidyverse)

#### Load parasite data----

## Import parasite dataset

full_dataset <- read_csv("data/processed/Full_dataset_with_psite_life_history_info_2024.11.21.csv")

# What is the number of fish we dissected per fish species
fish_dissected <- full_dataset %>% group_by(Fish_sp.x) %>% 
  summarise(number_dissected=length(unique(IndividualFishID)),
            .groups = 'drop') 

# What is the size range of fish we dissected per fish species
fish_oneperfish <- full_dataset %>% group_by(Fish_sp.x,IndividualFishID) %>% 
  summarise(TL_mm=mean(TotalLength_mm),
            .groups = 'drop') 

fish_sizerange <- fish_oneperfish %>% group_by(Fish_sp.x) %>% 
  summarise(mean_TL=mean(TL_mm),
            min_TL=min(TL_mm),
            max_TL=max(TL_mm),
            .groups = 'drop') 

full_dataset_na <- subset(full_dataset, !is.na(psite_count))

# Check number of trematode metacercariae counted --- this is for leopoldina. Not for this analysis

#trem_cest <- subset(full_dataset,Parasite_taxonomic_group == "Cestoda"|
                      #Parasite_taxonomic_group == "Trematoda")
  
#total_tremcest <- trem_cest %>% group_by(fish_psite_combo) %>% 
  #summarise(number_counted=sum(psite_count),
            #.groups = 'drop') 

#total_tremcest_high <- subset(total_tremcest,number_counted >= 20)

#full_dataset_na <- subset(full_dataset, !is.na(psite_count))


# What is the number of fish we dissected per fish species
parasites_counted <- full_dataset_na %>% group_by(Parasite_taxonomic_group) %>% 
  summarise(number_counted=sum(psite_count),
            .groups = 'drop') 

parasites_counted_all <- full_dataset_na %>%
  summarise(number_counted=sum(psite_count),
            .groups = 'drop') 

result <- full_dataset_na %>%
  group_by(Fish_sp.y, YearCollected, MonthCollected) %>%
  summarise(observations = n_distinct(IndividualFishID), .groups = "drop")

# View the result
View(result)

## Create column for presence and absegroup_by()## Create column for presence and absence of myxozoans
full_dataset$psite_presence <- ifelse(full_dataset$psite_count > 0,                # condition
                                           1,    # what if condition is TRUE
                                           0)       # what if condition is FALSE

# From Carpiodes the only other parasite that was present as well and might compete with MYX.G for space if CILI.FU 
cilimyxo <- subset(full_dataset,psite_spp.x=="MYX.G"|
                     psite_spp.x=="CILI.FU")

cili <- subset(full_dataset,psite_spp.x=="CILI.FU")


# Set a theme for all your plots
apatheme= theme_bw(base_size = 14,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

apatheme2= theme_bw(base_size = 24,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



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

library(lubridate)

# Parasite abundance against time

ggplot(mixN_unitall, aes(x= YearCollected,
                              y=Result))+
  geom_point()+apatheme+ylab("# myxozoan cysts/fish")+
  facet_wrap("Unit",scales = "free")+
  ggtitle("Cyst number per fish")

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

# View updated 'full_dataset'
View(full_dataset)

# rescale nutrients
full_dataset$scaledmean_nitrogen <- scale(full_dataset$mean_nitrogen)


#### Add specific conductance----

## Import physical data

physical <- read_csv("data/Physicochemical/physical_resultphyschem.csv")

# Read USGS-gauge site metadata

siteinfo <-  read_csv("data/Physicochemical/station_info.csv")

# Add info on sampling sites (latitude, longitude, CI, etc.)

physical_withmetadata <- merge(physical, siteinfo, by.x = "MonitoringLocationIdentifier", by.y = "MonitoringLocationIdentifier", all.x = TRUE)

# Take only data on specific conductance

spec_cond <- subset(physical_withmetadata, Measure == "Specific conductance")


## The script below assigns to each fish individual the mean specific conductance (mean_conductivity) that the fish experienced over one year before its collection

# Add a new column to store the mean nutrient in the form of mixed nitrogen forms
full_dataset$mean_conductivity <- NA

library(lubridate)

# Combine year, month, and day columns into a date column
full_dataset$collection_date <- make_date(full_dataset$YearCollected, full_dataset$MonthCollected, 
                                          full_dataset$DayCollected)


spec_cond$measurement_date <- make_date(spec_cond$YearCollected, spec_cond$MonthCollected, 
                                        spec_cond$DayCollected)


# Loop through each row of 'full_dataset'
for (i in 1:nrow(full_dataset)) {
  # Extract the collection date for the individual
  collection_date <- as.Date(full_dataset$collection_date[i])
  
  # Define the start and end date for one year before the collection date
  start_date <- collection_date - 365
  end_date <- collection_date
  
  # Filter 'spec_cond' to get the nitrogen data within that one-year range
  conduct_data <- spec_cond[spec_cond$measurement_date >= start_date & spec_cond$measurement_date < end_date, "Result"]
  
  # Calculate the mean nitrogen and store it in 'full_dataset'
  if (length(conduct_data) > 0) {
    full_dataset$mean_conductivity[i] <- mean(conduct_data, na.rm = TRUE)
  } else {
    full_dataset$mean_conductivity[i] <- NA  # If no data is found, store NA
  }
}

# View updated 'full_dataset'
View(full_dataset)

# rescale nutrients
full_dataset$scaledmean_conductivity <- scale(full_dataset$mean_conductivity)








# Some data only spans between 1973 and 1994. So we will perform the analyses 
# only for that time period and for the control because I cannot assume that the 
# downstram impact sites are going to be represented by this data.
mixN_year <- subset(mixN_unitall,  
                    YearCollected > 1972 &
                      YearCollected < 1995 )
mixN_control <- subset(mixN_year, CI == "control")
mixN <- subset(mixN_control, Unit_full == "mg/l")








#### Add turbidity----

## Import physical data

physical <- read_csv("data/Physicochemical/physical_resultphyschem.csv")

# Read USGS-gauge site metadata

siteinfo <-  read_csv("data/Physicochemical/station_info.csv")

# Add info on sampling sites (latitude, longitude, CI, etc.)

physical_withmetadata <- merge(physical, siteinfo, by.x = "MonitoringLocationIdentifier", by.y = "MonitoringLocationIdentifier", all.x = TRUE)

# Take only data on turbidity

turbidity <- subset(physical_withmetadata, Measure == "Turbidity")


## The script below assigns to each fish individual the mean turbidity that the fish experienced over one year before its collection

# Add a new column to store the mean turbidity
full_dataset$mean_turbidity <- NA

library(lubridate)

# Combine year, month, and day columns into a date column
full_dataset$collection_date <- make_date(full_dataset$YearCollected, full_dataset$MonthCollected, 
                                          full_dataset$DayCollected)


turbidity$measurement_date <- make_date(turbidity$YearCollected, turbidity$MonthCollected, 
                                        turbidity$DayCollected)


# Loop through each row of 'full_dataset'
for (i in 1:nrow(full_dataset)) {
  # Extract the collection date for the individual
  collection_date <- as.Date(full_dataset$collection_date[i])
  
  # Define the start and end date for one year before the collection date
  start_date <- collection_date - 365
  end_date <- collection_date
  
  # Filter 'turbidity' to get the turbidity data within that one-year range
  turb_data <- turbidity[turbidity$measurement_date >= start_date & turbidity$measurement_date < end_date, "Result"]
  
  # Calculate the mean turbidity and store it in 'full_dataset'
  if (length(turb_data) > 0) {
    full_dataset$mean_turbidity[i] <- mean(turb_data, na.rm = TRUE)
  } else {
    full_dataset$mean_turbidity[i] <- NA  # If no data is found, store NA
  }
}

# View updated 'full_dataset'
View(full_dataset)

# rescale nutrients
full_dataset$scaledmean_turbidity <- scale(full_dataset$mean_turbidity)

#### Create daughter data frames----

full_dataset_myxo <- subset(full_dataset, Parasite_taxonomic_group == "Myxozoa")


# What is the number of fish we dissected per fish species
fish_dissected <- full_dataset_myxo %>% group_by(Fish_sp.x) %>% 
  summarise(number_dissected=length(unique(IndividualFishID)),
            .groups = 'drop') 


## Take the log of the size to later include it in the analysis as offset

full_dataset_myxo$logTL_mm <- log(full_dataset_myxo$TotalLength_mm)

## Summary stats
prevalence_myxo <- subset(full_dataset_myxo, !is.na(psite_presence))
abundance_myxo <- subset(full_dataset_myxo, !is.na(psite_count))

## Which parasites have a prevalence of more than 5%
summary_prevalence <- prevalence_myxo %>% group_by(Fish_sp.x,psite_spp.x) %>% 
  summarise(Prevalence=mean(psite_presence),
            .groups = 'drop') 

prevalence_higherthan5percent <- subset(summary_prevalence, Prevalence > 0.05)

# The parasites with more than 5% of prevalence are

## Take a look at abundances
summary_abundance <- abundance_myxo %>% group_by(Fish_sp.x,psite_spp.x) %>% 
  summarise(Mean_abundance=mean(psite_count),
            Max_abundance=max(psite_count),
.groups = 'drop') 

## Remove columns with myxozoans that aren't confirmed to be myxozoans and therefore are not identified
full_dataset_myxo <- subset(full_dataset_myxo, !is.na(Parasite_genus))
full_dataset_myxo <- subset(full_dataset_myxo, psite_spp.x!="MYX.GM") # According to Stephen Atkinson (OSU), this is likely a contamination and also it is in the gall bladder so does not count for abundace
full_dataset_myxo <- subset(full_dataset_myxo, psite_spp.x!="MYX.GEOM") # According to Stephen Atkinson (OSU), this is not a myxozoan

# Parasite abundance against time

ggplot(full_dataset_myxo, aes(x= YearCollected,
                       y=psite_count,color=Parasite_genus))+
  geom_point()+apatheme+ylab("# myxozoan cysts/fish")+
  facet_wrap("Fish_sp.x")+
  ggtitle("Cyst number per fish")
  

ggplot(full_dataset_myxo, aes(x= YearCollected,
                            y=psite_presence,color=Parasite_genus,
                            fill=Parasite_genus))+
  geom_bar(stat = "identity")+apatheme+
  ggtitle("Prevalence of infection")+
  xlab("Year")+ylab("Prevalence of infection (%)")+
  facet_wrap("Fish_sp.x")

full_dataset_myxo

## Make plot showing the prevalence of infection of fish per parasite
# What is the number of fish we dissected per fish species
fish_dissected <- full_dataset_myxo %>% group_by(Fish_sp.x) %>% 
  summarise(number_dissected=length(unique(IndividualFishID)),
            .groups = 'drop') 

# Calculate prevalence by parasite genus and snail species
summary_psitegenus <- prevalence_myxo %>%
  group_by(Parasite_genus, Fish_sp.x,IndividualFishID) %>%             # Group by parasite genus and snail species
  summarize(genus_numinfected = sum(psite_presence), .groups = "drop")


summary_psitegenus$psite_presence <- ifelse(summary_psitegenus$genus_numinfected > 0,                # condition
                                           1,    # what if condition is TRUE
                                           0)       # what if condition is FALSE


stuff<-subset(full_dataset_myxo,IndividualFishID == "1169777_02")
View(stuff)


# Calculate prevalence by parasite genus and snail species
summary_prevalence <- summary_psitegenus %>%
  group_by(Parasite_genus, Fish_sp.x) %>%             # Group by parasite genus and snail species
  summarize(Prevalence = (mean(psite_presence)*100), .groups = "drop")



ggplot(summary_prevalence, aes(x= Fish_sp.x,
                              y=Prevalence,color=Parasite_genus,
                              fill=Parasite_genus))+
  geom_bar(stat = "identity")+apatheme2+
  ggtitle(NULL)+
  scale_color_manual(name = "Parasite:",
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
    scale_fill_manual(name = "Parasite:",
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
  
ggsave(file="Manuscripts/Figures/Year/overviewfish.png", width=240, height=180, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/overviewfish.pdf", width=240, height=180, dpi=1000, units = "mm")

# Decompose CARVEL

length(unique(full_carvel$IndividualFishID))

full_carvel <- subset(full_dataset_myxo,Fish_sp.x=="Carpiodes velifer")

full_carvel_inf <- subset(full_carvel,psite_count>0)

full_carvel_inf_myxobolus <- subset(full_carvel_inf,Parasite_genus== "Myxobolus")

length(unique(full_carvel_inf_myxobolus$IndividualFishID))

ggplot(full_carvel, aes(x= psite_spp.x,
                              y=psite_presence,color=Parasite_genus,
                              fill=Parasite_genus))+
  geom_bar(stat = "identity")+apatheme2+
  ggtitle(NULL)+
  scale_color_manual(name = "Parasite:",
                     values = c("Chloromyxum" = "#543005", 
                                "Myxobolus" = "#c7eae5"),
                     labels = c("Chloromyxum" = expression(italic("Chloromyxum")*" sp."), 
                                "Myxobolus" = expression(italic("Myxobolus")*" spp."))) +
  scale_fill_manual(name = "Parasite:",
                    values = c("Chloromyxum" = "#543005", 
                               "Myxobolus" = "#c7eae5"), 
                    labels = c("Chloromyxum" = expression(italic("Chloromyxum")*" sp."), 
                               "Myxobolus" = expression(italic("Myxobolus")*" spp."))) +
  xlab("Fish species")+ylab("Prevalence of infection (%)")


#Revise variable types
full_dataset_myxo$IndividualFishID <- as.factor(full_dataset_myxo$IndividualFishID)
full_dataset_myxo$CI <- as.factor(full_dataset_myxo$CI)
full_dataset_myxo$Fish_sp.x <- as.factor(full_dataset_myxo$Fish_sp.x)
full_dataset_myxo$Parasite_genus <- as.factor(full_dataset_myxo$Parasite_genus)
full_dataset_myxo$before_after <- as.factor(full_dataset_myxo$before_after)
full_dataset_myxo$season <- as.factor(full_dataset_myxo$season)
full_dataset_myxo$site <- as.factor(full_dataset_myxo$site)

## Revise reference levels

full_dataset_myxo$before_after <- relevel(full_dataset_myxo$before_after, ref = "before")
full_dataset_myxo$CI <- relevel(full_dataset_myxo$CI, ref = "control")

## Create subset per fish species

myx_cm <- subset(full_dataset_myxo, psite_spp.x == "MYX.CM")
myx_f <- subset(full_dataset_myxo, psite_spp.x == "MYX.F")
myx_g <- subset(full_dataset_myxo, psite_spp.x == "MYX.G")
myx_stwy <- subset(full_dataset_myxo, psite_spp.x == "MYX.SWTY")
myx_tail <- subset(full_dataset_myxo, psite_spp.x == "MYX.TAIL")
myx_sp <- subset(full_dataset_myxo, psite_spp.x == "MYX.SP")
myx_go <- subset(full_dataset_myxo, psite_spp.x == "MYX.GO")
myx_thel <- subset(full_dataset_myxo, psite_spp.x == "MYX.THEL")
myxo_sbad <- subset(full_dataset_myxo, psite_spp.x == "MYXO.SBAD")

## Crease subset for myxozoans that can be counted

myxo_count <- subset(full_dataset_myxo, Parasite_genus == "Myxobolus"|Parasite_genus == "Henneguya"|Parasite_genus == "Unicauda"|Parasite_genus == "Thelohanellus")

## Create subset of myxo_count per fish species

pimvig_count <- subset(myxo_count, Fish_sp.x == "Pimephales vigilax")
gamaff_count <- subset(myxo_count, Fish_sp.x == "Gambusia affinis")
ictpun_count <- subset(myxo_count, Fish_sp.x == "Ictalurus punctatus")
notath_count <- subset(myxo_count, Fish_sp.x == "Notropis atherinoides")
carvel_count <- subset(myxo_count, Fish_sp.x == "Carpiodes velifer")
hybnuc_count <- subset(myxo_count, Fish_sp.x == "Hybognathus nuchalis")

## Restrict to year and control for the interactions analysis

# Restrict to the years that we are going to work with
myxo_count_y <- subset(myxo_count, 
                             YearCollected > 1972 &
                               YearCollected < 1995 )

# Restrict to control because we have environmental data for that area of the river
myxo_count_yc <- subset(myxo_count_y, 
                                    CI == "control" )

## Create subset of myxo_count_ycm per fish species

pimvig_count_ycm <- subset(myxo_count_ycm, Fish_sp.x == "Pimephales vigilax")
gamaff_count_ycm <- subset(myxo_count_ycm, Fish_sp.x == "Gambusia affinis")
ictpun_count_ycm <- subset(myxo_count_ycm, Fish_sp.x == "Ictalurus punctatus")
notath_count_ycm <- subset(myxo_count_ycm, Fish_sp.x == "Notropis atherinoides")
carvel_count_ycm <- subset(myxo_count_ycm, Fish_sp.x == "Carpiodes velifer")
hybnuc_count_ycm <- subset(myxo_count_ycm, Fish_sp.x == "Hybognathus nuchalis")

# Set a theme for all your plots
apatheme= theme_bw(base_size = 14,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## Crease subset for myxozoans that cannot be counted

myxo_presence <- subset(full_dataset_myxo, Parasite_genus == "Chloromyxum"|Parasite_genus == "Myxidium")

## Create subset of myxo_presence per fish species

pimvig_presence <- subset(myxo_presence, Fish_sp.x == "Pimephales vigilax")
gamaff_presence <- subset(myxo_presence, Fish_sp.x == "Gambusia affinis")
ictpun_presence <- subset(myxo_presence, Fish_sp.x == "Ictalurus punctatus")
notath_presence <- subset(myxo_presence, Fish_sp.x == "Notropis atherinoides")
carvel_presence <- subset(myxo_presence, Fish_sp.x == "Carpiodes velifer")
hybnuc_presence <- subset(myxo_presence, Fish_sp.x == "Hybognathus nuchalis")

#### Add element concentrations----
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

# take only the years and control sites
# Make a subset for element concentration in water

elements_water_w <- subset(elements_na_wide, Unit=="ug/l"|
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

elements_na_wide_num <- myxo_count_yc[-c(1:44)]

# scale variables

elements_na_wide_scaled <- scale(elements_na_wide_num)

# Estimate the number of dimensions for the Principal Component Anal- ysis by cross-validation
# The number of components which leads to the smallest mean square error of prediction (MSEP) is retained. 
# For the Kfold cross-validation, pNA percentage of missing values is inserted and predicted with a PCA model using ncp.min to ncp.max dimensions.

library(missMDA)

#nb <- estim_ncpPCA(elements_na_wide_scaled,method.cv = "Kfold", verbose = FALSE)
nb$ncp #5

plot(0:5, nb$criterion, xlab = "nb dim", ylab = "MSEP")

# Because we have a lot of NAs and PCAs can't handle that, we have to impute values based on our data
# The (regularized) iterative PCA algorithm (by default) first consists imputing missing values with initial values such as the mean of the variable.
# It is adviced to use the regularized version of the algorithm to avoid the overfitting problems which are very frequent when there are many missing values. 

library(FactoMineR)

res.comp <- imputePCA(elements_na_wide_scaled, ncp = 2) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set

imp <- cbind.data.frame(res.comp$completeObs)

res.pca <- PCA(imp, quanti.sup = 1, quali.sup = 12, ncp = 2, graph=FALSE) # generate PCA

plot(res.pca, choix="var")+apatheme2# visualize PCA

# Extract contribution of variables to the components

contributions_pca <- res.pca$var$contrib
loadings_pca <- res.pca$var$coord

# Add principal components to your data frame
score <- as_tibble(factoextra::get_pca_ind(res.pca)$coord) #extract individual scores to be used in glm

myxo_count_ycm <- cbind(myxo_count_yc, score[1:2]) # merge scores with original data

myxo_count_ycm$Elements_PC1 <- myxo_count_ycm$Dim.1
myxo_count_ycm$Elements_PC2 <- myxo_count_ycm$Dim.2

#### Add element concentrations, PCA with N----
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

# take only the years and control sites
# Make a subset for element concentration in water

elements_water_w <- subset(elements_na_wide, Unit=="ug/l"|
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

elements_na_wide_num <- myxo_count_yc[-c(1:36,38:44)]

# scale variables

elements_na_wide_scaled <- scale(elements_na_wide_num)

# Estimate the number of dimensions for the Principal Component Anal- ysis by cross-validation
# The number of components which leads to the smallest mean square error of prediction (MSEP) is retained. 
# For the Kfold cross-validation, pNA percentage of missing values is inserted and predicted with a PCA model using ncp.min to ncp.max dimensions.

library(missMDA)

#nb <- estim_ncpPCA(elements_na_wide_scaled,method.cv = "Kfold", verbose = FALSE)
#nb$ncp #5

plot(0:5, nb$criterion, xlab = "nb dim", ylab = "MSEP")

# Because we have a lot of NAs and PCAs can't handle that, we have to impute values based on our data
# The (regularized) iterative PCA algorithm (by default) first consists imputing missing values with initial values such as the mean of the variable.
# It is adviced to use the regularized version of the algorithm to avoid the overfitting problems which are very frequent when there are many missing values. 

library(FactoMineR)

res.comp <- imputePCA(elements_na_wide_scaled, ncp = 2) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set

imp <- cbind.data.frame(res.comp$completeObs)

res.pca <- PCA(imp, quanti.sup = 1, quali.sup = 11, ncp = 2, graph=FALSE) # generate PCA

plot(res.pca, choix="var") # visualize PCA

score <- as_tibble(factoextra::get_pca_ind(res.pca)$coord) #extract individual scores to be used in glm

myxo_count_ycm <- cbind(myxo_count_yc, score[1:2]) # merge scores with original data


#### The effect of time on abundance data----

## Carpiodes velifer - MYX.G

# Water temperature per month

ggplot(myx_g, aes(x= as.numeric(YearCollected),
                  y=psite_count,color=season,fill=season))+
  geom_point()+apatheme+
  xlab("Year")+ylab("Parasite abundance (# pseudocysts/fish)")

ggplot(myx_g, aes(x= as.numeric(YearCollected),
                  y=psite_count))+
  geom_point()+apatheme+
  facet_wrap("season")+
  xlab("Year")+ylab("Parasite abundance (# pseudocysts/fish)")



# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_g,
              family = nbinom1(link="log")) # AIC = 1532

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_g,
              family = nbinom2(link="log")) # AIC = 1538


m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myx_g,
              family = nbinom1(link="sqrt")) # AIC = 1610

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myx_g,
              family = nbinom2(link="sqrt")) # AIC = 1538

summary(m2)
AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

# m2 has the best diagnostics and an acceptable AIC. It is 6 units higher than m1, but it is a good trade-off for better diagnostics 

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                        offset(logTL_mm)+
                        (1|CI/site)+
                        (1|season),
                      data = myx_g,
                      family = nbinom2(link="log"))

## Extract estimates and confidence intervals for myx_g

# Calculate confidence intervals and estimates
conf_int_myxg <- confint(myx_g_time)

# Combine estimates and confidence intervals into a data frame
results_myxg <- data.frame(
  Estimate = conf_int_myxg[2,3],
  Upper_CI = conf_int_myxg[2,2],
  Lower_CI = conf_int_myxg[2,1]
)

results_myxg$Fish_sp <- "Carpiodes velifer"
results_myxg$Parasite <- "MYX.G"

#Check over dispersion
performance::check_overdispersion(m2)

#Evaluate model
tab_model(m2)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
plot_model(m2,type = "re")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m2, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/Carvel_MyxG_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/Carvel_MyxG_Year.pdf", width=150, height=120, dpi=1000, units = "mm")


## Carpiodes velifer - MYX.F

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_f,
              family = nbinom1(link="log")) # 269

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_f,
              family = nbinom2(link="log")) # 269

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myx_f,
              family = nbinom1(link="sqrt")) # AIC = 300

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myx_f,
              family = nbinom2(link="sqrt")) # did not converge

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

# m2 has the best AIC and the best residuals
myx_f_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                        offset(logTL_mm)+
                        (1|CI/site)+
                        (1|season),
                      data = myx_f,
                      family = nbinom2(link="log"))


## Extract estimates and confidence intervals for myx_f

# Calculate confidence intervals and estimates
conf_int_myxf <- confint(myx_f_time)

# Combine estimates and confidence intervals into a data frame
results_myxf <- data.frame(
  Estimate = conf_int_myxf[2,3],
  Upper_CI = conf_int_myxf[2,2],
  Lower_CI = conf_int_myxf[2,1]
)

results_myxf$Fish_sp <- "Carpiodes velifer"
results_myxf$Parasite <- "MYX.F"

  
#Check over dispersion
performance::check_overdispersion(m2)

#Evaluate model
tab_model(m2)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m2, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/Carvel_MyxF_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/Carvel_MyxF_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

### Ictalurus punctatus
## Ictalurus puncatus - MYX.TAIL

# GLMM - which family is the best?

m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_tail,
              family = nbinom1(link="log")) # did not converge

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_tail,
              family = nbinom2(link="log")) # AIC = 90

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myx_tail,
              family = nbinom1(link="sqrt")) # AIC = 91

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myx_tail,
              family = nbinom2(link="sqrt")) # Did not converge

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m3,n=250)
s$scaledResiduals
plot(s)

# m3 has the best AIC and diagnostics

myx_tail_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                           offset(sqrt(TotalLength_mm))+
                           (1|CI/site)+
                           (1|season),
                         data = myx_tail,
                         family = nbinom1(link="sqrt"))

## Extract estimates and confidence intervals for myx_f

# Calculate confidence intervals and estimates
conf_int_myxtail <- confint(myx_tail_time)

# Combine estimates and confidence intervals into a data frame
results_myxtail <- data.frame(
  Estimate = conf_int_myxtail[2,3],
  Upper_CI = conf_int_myxtail[2,2],
  Lower_CI = conf_int_myxtail[2,1]
)

results_myxtail$Fish_sp <- "Ictalurus puncatus"
results_myxtail$Parasite <- "MYX.TAIL"


#Check over dispersion
performance::check_overdispersion(m2)

#Evaluate model
tab_model(m3)
plot_model(m3,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m3, terms= c("YearCollected")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/ICTPUN_MyxTail_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/ICTPUN_MyxTail_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

### Notropis atherinoides
## Notropis atherinoides - MYX.SP
# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_sp,
              family = nbinom1(link="log")) # AIC = 482

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_sp,
              family = nbinom2(link="log")) # AIC = 479

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myx_sp,
              family = nbinom1(link="sqrt")) # AIC = 491

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myx_sp,
              family = nbinom2(link="sqrt")) # AIC = 500

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

# m2 has the best AIC and diagnostics

myx_sp_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                         offset(logTL_mm)+
                         (1|CI/site)+
                         (1|season),
                       data = myx_sp,
                       family = nbinom2(link="log"))

## Extract estimates and confidence intervals for myx_f

# Calculate confidence intervals and estimates
conf_int_myxsp <- confint(myx_sp_time)

# Combine estimates and confidence intervals into a data frame
results_myxsp <- data.frame(
  Estimate = conf_int_myxsp[2,3],
  Upper_CI = conf_int_myxsp[2,2],
  Lower_CI = conf_int_myxsp[2,1]
)

results_myxsp$Fish_sp <- "Notropis atherinoides"
results_myxsp$Parasite <- "MYX.SP"


#Check over dispersion
performance::check_overdispersion(m2)

#Evaluate model
tab_model(m2)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m2, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/NOTATH_MyxSP_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/NOTATH_MyxSP_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

### Pimephales vigilax
## Pimephales vigilax - MYX.GO

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_go,
              family = nbinom1(link="log")) # AIC = 137

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_go,
              family = nbinom2(link="log")) # did not converge

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myx_go,
              family = nbinom1(link="sqrt")) # AIC = 151

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myx_go,
              family = nbinom2(link="sqrt")) # did not converge

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

#Check over dispersion
performance::check_overdispersion(m1)

# m1 has the best AIC and the best diagnostics, although there is a significant pattern in the residuals

myx_go_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                         offset(logTL_mm)+
                         (1|CI/site)+
                         (1|season),
                       data = myx_go,
                       family = nbinom1(link="log"))

## Extract estimates and confidence intervals for myx_f

# Calculate confidence intervals and estimates
conf_int_myxgo <- confint(myx_go_time)

# Combine estimates and confidence intervals into a data frame
results_myxgo <- data.frame(
  Estimate = conf_int_myxgo[2,3],
  Upper_CI = conf_int_myxgo[2,2],
  Lower_CI = conf_int_myxgo[2,1]
)

results_myxgo$Fish_sp <- "Pimephales vigilax"
results_myxgo$Parasite <- "MYX.GO"


#Evaluate model
tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m1, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/PIMVIG_MyxGO_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/PIMVIG_MyxGO_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")


### Pimephales vigilax
## Pimephales vigilax - MYX.THEL

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_thel,
              family = nbinom1(link="log")) # AIC = 279

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_thel,
              family = nbinom2(link="log")) # AIC = 274

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myx_thel,
              family = nbinom1(link="sqrt")) # AIC = 277

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myx_thel,
              family = nbinom2(link="sqrt")) # did not converge

AIC(m1,m2,m3,m4)

# Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

# Check over dispersion
performance::check_overdispersion(m1)

# m1 has the best diagnostics. It is 5 units higher in terms of AIC than m3, but the diagnostics of m3 are bad. 

myx_thel_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                           offset(logTL_mm)+
                           (1|CI/site)+
                           (1|season),
                         data = myx_thel,
                         family = nbinom1(link="log"))

## Extract estimates and confidence intervals for myx_f

# Calculate confidence intervals and estimates
conf_int_myxthel <- confint(myx_thel_time)

# Combine estimates and confidence intervals into a data frame
results_myxthel <- data.frame(
  Estimate = conf_int_myxthel[2,3],
  Upper_CI = conf_int_myxthel[2,2],
  Lower_CI = conf_int_myxthel[2,1]
)

results_myxthel$Fish_sp <- "Pimephales vigilax"
results_myxthel$Parasite <- "MYX.THEL"


#Evaluate model
tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m1, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/PIMVIG_MyxTHEL_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/PIMVIG_MyxTHEL_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

### Pimephales vigilax
## Pimephales vigilax - MYXO.SBAD

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myxo_sbad,
              family = nbinom1(link="log")) # AIC = 191

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myxo_sbad,
              family = nbinom2(link="log")) # AIC = 188

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myxo_sbad,
              family = nbinom1(link="sqrt")) # AIC = 187

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|CI/site)+
                (1|season),
              data = myxo_sbad,
              family = nbinom2(link="sqrt")) # did not converge

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

# m2 has the best AIC and the best diagnostics

myxo_sbad_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                            offset(logTL_mm)+
                            (1|CI/site)+
                            (1|season),
                          data = myxo_sbad,
                          family = nbinom2(link="log")) 

## Extract estimates and confidence intervals for myx_sbad

# Calculate confidence intervals and estimates
conf_int_myxsbad <- confint(myxo_sbad_time)

# Combine estimates and confidence intervals into a data frame
results_myxsbad <- data.frame(
  Estimate = conf_int_myxsbad[2,3],
  Upper_CI = conf_int_myxsbad[2,2],
  Lower_CI = conf_int_myxsbad[2,1]
)

results_myxsbad$Fish_sp <- "Pimephales vigilax"
results_myxsbad$Parasite <- "MYX.SBAD"


#Check over dispersion
performance::check_overdispersion(m2)

#Evaluate model
tab_model(m2)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m2, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/PIMVIG_MyxSBAD_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/PIMVIG_MyxSBAD_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")


## Combine all estimates and confidence intervals to make summary plot

## combine  data frames

results_combined <- bind_rows(results_myxg,results_myxf, results_myxtail,results_myxsp,results_myxthel,results_myxsbad,results_myxgo)

# What is the mean estimate across all fishes?
estimates_abundance_mean <- results_combined %>%  
  summarise(Est_mean=mean(Estimate),
            .groups = 'drop') 

# Data for horizontal lines: the mean estimate observed across all fish species. This helps you to see which fish is significantly deviating from the rest of the fish
# I took the values "-0.936" and "510" from the code lines 168-170.

#dummy2 <- data.frame(Predictor = c("psite_count", "CIimpact"), Z = c(-0.936, 510))

ggplot(dummy1, aes(x = D, y = Y)) + geom_point() + facet_grid(~X) + 
  geom_hline(data = dummy2, aes(yintercept = Z))

# Set a theme for all your plots
apatheme2= theme_bw(base_size = 24,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# Plot the estimates with confidence intervals
ggplot(results_combined, aes(x = Parasite, y = Estimate, ymin = Lower_CI, ymax = Upper_CI)) +
  geom_pointrange(size=0.6) +
  labs(x = "Host - parasite", y = "Effect of time") +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1.0,color="darkgrey")+apatheme2+coord_flip()+
  scale_x_discrete(labels = c("MYX.G" = expression(italic("Myxobolus")*" sp. 2 "*italic("- C. velifer ")*"gills"),
                              "MYX.F" = expression(italic("Myxobolus")*" sp. 1 "*italic("- C. velifer ")*"fins"),
                              "MYX.THEL" = expression(italic("Thelohanellus")*" sp."*italic("- P. vigilax ")*"gills"),
                              "MYX.SBAD" = expression(italic("Myxobolus")*" sp. 4"*italic("- P. vigilax ")*"fins"),
                              "MYX.GO" = expression(italic("Myxobolus")*" sp. 5 "*italic("- P. vigilax ")*"gills"),
                              "MYX.TAIL" = expression(italic("Henneguya")*" sp."*italic("- I. punctatus ")*"gills"),
                              "MYX.SP" = expression(italic("Myxobolus")*" sp. 3" *italic("- N. atherinoides ")*"fins")),
                   limits = c("MYX.GO","MYX.SBAD","MYX.THEL","MYX.SP","MYX.TAIL", "MYX.G","MYX.F"))


#geom_hline(data = dummy2, aes(yintercept = Z),linetype = "dashed")+coord_flip()+apatheme

#### The effect of time on prevalence data----

## Carpiodes velifer - MYX.CM

# GLMM - which family is the best?
m1 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_cm,
              family = binomial(link = "logit")) 

m2 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_cm,
              family = binomial(link = "probit")) 

# This one coverged and it's good
myxo_cm_time <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = myx_cm,
              family = binomial(link = "probit")) 


## Extract estimates and confidence intervals for myx_cm

# Calculate confidence intervals and estimates
conf_int_myxscm <- confint(myxo_cm_time)

# Combine estimates and confidence intervals into a data frame
results_myxscm <- data.frame(
  Estimate = conf_int_myxscm[2,3],
  Upper_CI = conf_int_myxscm[2,2],
  Lower_CI = conf_int_myxscm[2,1]
)

results_myxscm$Fish_sp <- "Carpiodes velifer"
results_myxscm$Parasite <- "MYX.CM"


#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

#Check over dispersion
performance::check_overdispersion(m2)

#Evaluate model
tab_model(m2)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m2, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Probability of myxozoan infection (%)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/Carvel_MyxCM_Presence_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/Carvel_MyxCM_Presence_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

## Gambusia affinis - MYX.STWY
gamaff_count_myxstwy <- subset(gamaff_presence,psite_spp.x == "MYX.STWY")

# GLMM - which family is the best?
m1 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = gamaff_count_myxstwy,
              family = binomial(link = "logit")) 

m2 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = gamaff_count_myxstwy,
              family = binomial(link = "probit")) 

m3 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = gamaff_count_myxstwy,
              family = binomial(link = "identity")) 

m4 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = gamaff_count_myxstwy,
              family = binomial(link = "inverse")) 

# GLMM - None of them converged so we try to remove CI from RE
m1 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = gamaff_count_myxstwy,
              family = binomial(link = "logit")) 

m2 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = gamaff_count_myxstwy,
              family = binomial(link = "probit")) 

m3 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = gamaff_count_myxstwy,
              family = binomial(link = "identity")) 

m4 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = gamaff_count_myxstwy,
              family = binomial(link = "inverse")) 

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

# The best is m1

myxo_stwy_time <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = gamaff_count_myxstwy,
              family = binomial(link = "logit")) 

## Extract estimates and confidence intervals for myx_cm

# Calculate confidence intervals and estimates
conf_int_myxstwy <- confint(myxo_stwy_time)

# Combine estimates and confidence intervals into a data frame
results_myxstwy <- data.frame(
  Estimate = conf_int_myxstwy[2,3],
  Upper_CI = conf_int_myxstwy[2,2],
  Lower_CI = conf_int_myxstwy[2,1]
)

results_myxstwy$Fish_sp <- "Gambusia affinis"
results_myxstwy$Parasite <- "MYX.STWY"

#Check over dispersion
performance::check_overdispersion(m1)

#Evaluate model
tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m1, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Probability of myxozoan infection (%)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/GAMAFF_MyxSTWY_Presence_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/GAMAFF_MyxSTWY_Presence_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

#### Non-point source - Only for MYX.G----

# Some data only spans between 1973 and 1994. So we will perform the analyses 
# only for that time period and for the control because I cannot assume that the 
# downstram impact sites are going to be represented by this data.

myx_g_yc <- subset(myxo_count_ycm, psite_spp.x == "MYX.G")



m1 <- glmmTMB(psite_count ~ 
                scale(mean_temperature)+
                Elements_PC1+
                Elements_PC2+
                scale(mean_nitrogen)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|YearCollected/season),
              data = myx_g_yc,
              family=nbinom1(link="sqrt"))

vif_values <- check_collinearity(m1)


m1 <- glmmTMB(psite_count ~ 
                scale(mean_temperature)*scale(mean_nitrogen)+
                scale(mean_temperature)*Elements_PC1+
                scale(mean_temperature)*Elements_PC2+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
               data = myx_g_yc,
              family=nbinom1(link="log")) #nbinom2 log bad residuals, nibinom1 log good residuals

m2 <- glmmTMB(psite_count ~ 
                scale(mean_temperature)*scale(mean_nitrogen)+
                scale(mean_temperature)*Elements_PC1+
                scale(mean_temperature)*Elements_PC2+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|YearCollected/season),
              data = myx_g_yc,
              family=nbinom2(link="sqrt")) #nbinom2 sqrt bad residuals, nbinom1 sqrt did not converge

m3 <- glmmTMB(psite_count ~ 
                scale(mean_temperature)*scale(mean_nitrogen)+
                scale(mean_temperature)*Elements_PC1+
                scale(mean_temperature)*Elements_PC2+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|YearCollected/season),
              data = myx_g_yc,
              family=nbinom1(link="sqrt")) #nbinom2 sqrt bad residuals, nbinom1 sqrt did not converge

vif_values <- check_collinearity(m2)

AIC(m1,m2,m3)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m2)

# The best model is 

ms_model <- glmmTMB(psite_count ~ 
                scale(mean_temperature)*scale(mean_nitrogen)+
                scale(mean_temperature)*Elements_PC1+
                scale(mean_temperature)*Elements_PC2+
                offset(sqrt(TotalLength_mm))+
                (1|site/CatalogNumber)+
                (1|YearCollected/season),
              data = myx_g_yc,
              family=nbinom2(link="sqrt")) #nbinom2 sqrt bad residuals, nbinom1 sqrt did not converge

ms_model <- glmmTMB(psite_count ~ 
                      scale(mean_temperature)*scale(mean_nitrogen)+
                      scale(mean_temperature)*Elements_PC1+
                      scale(mean_temperature)*Elements_PC2+
                      cili_fu+
                      offset(sqrt(TotalLength_mm))+
                      (1|site/CatalogNumber)+
                      (1|YearCollected/season),
                    data = myx_g_yc,
                    family=nbinom2(link="sqrt")) #nbinom2 sqrt bad residuals, nbinom1 sqrt did not converge

tab_model(ms_model)
plot_model(ms_model,type = "est",ci.lvl = 0.95,colors = c("black"))+apatheme2+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m2)

# Plot predictions

mydf <- ggpredict(ms_model, terms= c("mean_temperature[n=50]")) 

plot(mydf,show_data=TRUE,show_residuals=TRUE,color=c("#5aae61","#762a83"))+
  labs(x = "Temperature (°C)", y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme2

# Plot predictions

mydf <- ggpredict(ms_model, terms= c("mean_nitrogen[n=50]")) 

plot(mydf,show_data=TRUE,show_residuals=TRUE,color=c("#5aae61","#762a83"))+
  labs(x = "Nitrogen (mg/L)", y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme2

  
# Plot predictions

mydf <- ggpredict(ms_model, terms= c("Elements_PC1[n=50]")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,color=c("#5aae61"))+
labs(x = "Elements PC1", y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme2


# Plot predictions

mydf <- ggpredict(ms_model, terms= c("Elements_PC2[n=50]")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,color=c("#5aae61"))+
  labs(x = "Elements PC2", y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme2

# Plot predictions

mydf <- ggpredict(ms_model, terms= c("Elements_PC1[n=100]","mean_temperature")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,color=c("#74add1","#fee090","#d73027"))+
  labs(color="Temperature (°C)",x = "Elements PC1", y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme2

# Plot predictions

mydf <- ggpredict(ms_model, terms= c("Elements_PC2[n=100]","mean_temperature")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,color=c("#74add1","#fee090","#d73027"))+
  labs(color="Temperature (°C)",x = "Elements PC2", y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme2

# Plot interactions Elements PC1

library(interactions)

max(myx_g_yc$mean_temperature)

interact_plot(m3, pred = Elements_PC1, modx = mean_temperature,
              interval = TRUE,
              colors = c("#2166ac", "#d1e5f0", "#b2182b"),
              int.width = 0.95,
              partial.residuals = FALSE,
              plot.points=TRUE,
              facet.modx = FALSE,
              point.alpha  = 0.5,
              x.label = "Elements PC1",
              y.label = "Parasite abundance (# pseudocysts/fish)",
              legend.main = "Temperature")+
  apatheme2

# Plot interactions Elements PC2

library(interactions)

interact_plot(ms_model, pred = Elements_PC2, modx = mean_temperature,
              interval = FALSE,
              colors = c("#2166ac", "#d1e5f0", "#b2182b"),
              int.width = 0.95,
              partial.residuals = FALSE,
              plot.points=TRUE,
              facet.modx = FALSE,
              point.alpha  = 0.5,
              x.label = "Elements PC2",
              y.label = "Parasite abudance (# pseudocysts/fish)",
              legend.main = "Temperature")+
  apatheme2

loadings_pca <- as.data.frame(loadings_pca)
loadings_pca$Element <- c("Ba","Cd","Cr","Cu","Pb","Fe","Mg","Mn","Hg","Ni")
colnames(loadings_pca)[1]<-"Elements_PC1"
colnames(loadings_pca)[2]<-"Elements_PC2"

# Plot the loadings for PC1 and PC2
ggplot(loadings_pca, aes(x = Element, y = Elements_PC1)) +
  geom_bar(stat="identity",fill="darkgrey") +
  labs(x = "Element", y = "Loadings") +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1.0,color="black")+apatheme2+coord_flip()

ggplot(loadings_pca, aes(x = Element, y = Elements_PC2)) +
  geom_bar(stat="identity",fill="darkgrey") +
  labs(x = "Element", y = "Loadings") +
  geom_hline(yintercept = 0,linetype = "dashed",size = 1.0,color="black")+apatheme2+coord_flip()

# Seems that lead might be partially responsible for the changes in myxozoan abundace
model_lead1 <- glmer.nb(psite_count ~ 
                        scale(mean_temperature)*scale(Pb)+
                        offset(logTL_mm)+
                        (1|site/CatalogNumber)+
                        (1|YearCollected/season),
                      data = myx_g_yc) #nbinom2 sqrt bad residuals, nbinom1 sqrt did not converge

model_lead2 <- glmmTMB(psite_count ~ 
                      scale(mean_temperature)*scale(Pb)+
                      offset(sqrt(TotalLength_mm))+
                      (1|site/CatalogNumber)+
                      (1|YearCollected/season),
                    data = myx_g_yc,
                    family=nbinom2(link="sqrt")) #nbinom2 sqrt bad residuals, nbinom1 sqrt did not converge

model_lead3 <- glmmTMB(psite_count ~ 
                         scale(mean_temperature)*scale(Pb)+
                         offset(sqrt(TotalLength_mm))+
                         (1|site/CatalogNumber)+
                         (1|YearCollected/season),
                       data = myx_g_yc,
                       family=nbinom1(link="sqrt")) #nbinom2 sqrt bad residuals, nbinom1 sqrt did not converge

model_lead4 <- glmmTMB(psite_count ~ 
                         scale(mean_temperature)*scale(Pb)+
                         offset(logTL_mm)+
                         (1|site/CatalogNumber)+
                         (1|YearCollected/season),
                       data = myx_g_yc,
                       family=nbinom1(link="log")) #nbinom2 sqrt bad residuals, nbinom1 sqrt did not converge

model_lead5 <- glmmTMB(psite_count ~ 
                         scale(mean_temperature)*scale(Pb)+
                         offset(logTL_mm)+
                         (1|site/CatalogNumber)+
                         (1|YearCollected/season),
                       data = myx_g_yc,
                       family=nbinom2(link="log")) #nbinom2 sqrt bad residuals, nbinom1 sqrt did not converge

AIC(model_lead1,model_lead2,model_lead3,model_lead4,model_lead5)

#Evaluate residuals
s=simulateResiduals(fittedModel=model_lead5,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m2)

# The best model is 


tab_model(model_lead1)
plot_model(model_lead1,type = "est",ci.lvl = 0.95,colors = c("black"))+apatheme2+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m2)


# Plot predictions

mydf <- ggpredict(model_lead2, terms= c("Pb[n=50]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,color=c("#5aae61","#762a83"))+
  labs(x = "Elements PC2", y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme2

# Plot interactions Elements PC2

library(interactions)

interact_plot(model_lead2, pred = Pb, modx = mean_temperature,
              interval = FALSE,
              colors = c("#2166ac", "#d1e5f0", "#b2182b"),
              int.width = 0.95,
              partial.residuals = FALSE,
              plot.points=TRUE,
              facet.modx = FALSE,
              point.alpha  = 0.82,
              x.label = "Elements PC2",
              y.label = "Parasite abudance (# pseudocysts/fish)",
              legend.main = "Temperature")+
  apatheme2


### Impact of point-source pollution - CARVEL MYX.G only CI----

m1 <- glmer.nb(psite_count~CI*before_after+
                 offset(logTL_mm)+
                 (1|site/CatalogNumber)+
                 (1|YearCollected/season),
               data = myx_g)

### Which link function is the best?
m1 <- glmmTMB(psite_count ~ 
                CI*before_after+
                offset(logTL_mm)+
                (1|site/CatalogNumber)+
                (1|YearCollected/season),
              data = myx_g,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ 
                CI*before_after+
                offset(sqrt(TotalLength_mm))+
                (1|site/CatalogNumber)+
                (1|YearCollected/season),
              data = myx_g,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ 
                CI*before_after+
                offset(logTL_mm)+
                (1|site/CatalogNumber)+
                (1|YearCollected/season),
              data = myx_g,
              family = nbinom2(link="log")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ 
                CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|YearCollected/season),
              data = myx_g,
              family = nbinom2(link="sqrt")) # with log-link function, nbinom1 does not converge

AIC(m1,m2,m3,m4) # best is nbinom1(sqrt)

summary(m1)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m3)

tab_model(m3)
plot_model(m3,type = "est",title = "",color="black")+
  apatheme2+
  geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

r.squaredGLMM(m3)

# Post hoc for interaction
post_hoc_interaction <- emmeans(m3, ~ CI * before_after)
pairwise_results_interaction <- pairs(post_hoc_interaction)
summary(pairwise_results_interaction)

## visualize model

# Flow with CI and psite_genus
mydf <- ggpredict(m3, terms= c("before_after","CI")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter = 0.1,
     dodge = 0.3,colors=c("#35978f","#8c510a"))+xlab("Clean Water Act")+
  ylab("Parasite abundance (# of pseudocysts per fish)")+
  apatheme2+ggtitle(NULL)

# Plot the predictions using ggplot2
ggplot(mydf, aes(x = x, y = predicted, color = group, group = group)) +
  geom_line(size = 1,position = position_dodge(0.5)) +               # Line for predicted values
  geom_point(size = 3,position = position_dodge(0.5)) +              # Points for predicted values
  geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0,position = position_dodge(0.5), alpha = 1.0) + 
  labs(title = "Predicted Response by Factor Levels",
       x = "Clean Water Act",
       y = "Parasite abundance (# of pseudocysts per fish)") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  scale_color_manual(values = c("#35978f","#8c510a")) +apatheme2+
  # Raw data points with jitter for better visualization
  geom_jitter(data = myx_g, aes(x = interaction(CI, before_after), y = psite_count), 
              width = 0.2, height = 0, color = "black", alpha = 0.5)


# Define custom colors
custom_colors <- c("A_X" = "#FF5733",  # Red for A_X
                   "A_Y" = "#33FF57",  # Green for A_Y
                   "B_X" = "#3357FF",  # Blue for B_X
                   "B_Y" = "#F0E68C")  # Khaki for B_Y

# Create the plot
ggplot() +
  # Predicted values with line and point for each group
  geom_line(data = mydf, aes(x = x, y = predicted, color = group, group = group),
            size = 1, position = position_dodge(0.5)) +
  geom_point(data = mydf, aes(x = x, y = predicted, color = group, group = group),
             size = 3, position = position_dodge(0.5)) +
  # Error bars without caps
  geom_errorbar(data = mydf, aes(x = x, ymin = conf.low, ymax = conf.high, color = group, group = group),
                width = 0,  # Remove caps
                position = position_dodge(0.5), 
                size = 1) +  # Custom thickness of error bars
  # Labels and theme
  labs(x = "Clean Water Act",
       y = "Parasite abundance (# of pseudocysts per fish)") +
  apatheme2+
  scale_color_manual(values = c("#35978f","#8c510a"))  # Apply custom colors


### Impact of point-source pollution - CARVEL MYX.F only CI----

carvel_count_myxf <- subset(carvel_count,psite_spp.x == "MYX.F")


### Which link function is the best?
m1 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = carvel_count_myxf,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = carvel_count_myxf,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = carvel_count_myxf,
              family = nbinom2(link="log")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = carvel_count_myxf,
              family = nbinom2(link="sqrt")) # with log-link function, nbinom1 does not converge

AIC(m1,m2,m3,m4) # best is nbinom2(log)

summary(m3)

#Evaluate residuals
s=simulateResiduals(fittedModel=m3,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m3)

tab_model(m3)
plot_model(m3,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m3)

## visualize model

# Flow with CI and psite_genus
mydf <- ggpredict(m3, terms= c("CI")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlab("Streamflow (m3/sec)")+ylab("Parasite abundance (# of pseudocysts per fish)")


### Impact of point-source pollution - ICTPUN MYX.TAIL only CI----

### Which link function is the best?
m1 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|CatalogNumber)+
                (1|YearCollected/season),
              data = ictpun_count_myxtail,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ 
                CI+
                offset(sqrt(TotalLength_mm))+
                (1|CatalogNumber)+
                (1|YearCollected/season),
              data = ictpun_count_myxtail,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|CatalogNumber)+
                (1|YearCollected/season),
              data = ictpun_count_myxtail,
              family = nbinom2(link="log")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ 
                CI+
                offset(sqrt(TotalLength_mm))+
                (1|CatalogNumber)+
                (1|YearCollected/season),
              data = ictpun_count_myxtail,
              family = nbinom2(link="sqrt")) # with log-link function, nbinom1 does not converge

AIC(m1,m2,m3,m4) # best is nbinom2(log)

summary(m3)

#Evaluate residuals
s=simulateResiduals(fittedModel=m3,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m3)

tab_model(m3)
plot_model(m3,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m3)

## visualize model

# Flow with CI and psite_genus
mydf <- ggpredict(m3, terms= c("CI")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlab("Streamflow (m3/sec)")+ylab("# of pseudocysts per fish")


### Impact of point-source pollution - NOTATH MYX.SP only CI----

notath_count_myxsp <- subset(notath_count,psite_spp.x == "MYX.SP")


### Which link function is the best?
m1 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = notath_count_myxsp,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = notath_count_myxsp,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = notath_count_myxsp,
              family = nbinom2(link="log")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = notath_count_myxsp,
              family = nbinom2(link="sqrt")) # with log-link function, nbinom1 does not converge

AIC(m1,m2,m3,m4) # best is nbinom2(log)

summary(m1)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m3)

## visualize model

# Flow with CI and psite_genus
mydf <- ggpredict(m1, terms= c("CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlab("Streamflow (m3/sec)")+ylab("# of pseudocysts per fish")+ylim(0,50)


### Impact of point-source pollution - PIMVIG MYX.GO only CI----

pimvig_count_myxgo <- subset(pimvig_count,psite_spp.x == "MYX.GO")


### Which link function is the best?
m1 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = pimvig_count_myxgo,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = pimvig_count_myxgo,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = pimvig_count_myxgo,
              family = nbinom2(link="log")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = pimvig_count_myxgo,
              family = nbinom2(link="sqrt")) # with log-link function, nbinom1 does not converge

AIC(m1,m2,m3,m4) # best is nbinom1(log)

summary(m1)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m3)

## visualize model

# Flow with CI and psite_genus
mydf <- ggpredict(m1, terms= c("CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlab("Streamflow (m3/sec)")+ylab("# of pseudocysts per fish")


### Impact of point-source pollution - PIMVIG MYX.THEL only CI----

pimvig_count_myxthel <- subset(pimvig_count,psite_spp.x == "MYX.THEL")


### Which link function is the best?
m1 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = pimvig_count_myxthel,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = pimvig_count_myxthel,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = pimvig_count_myxthel,
              family = nbinom2(link="log")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = pimvig_count_myxthel,
              family = nbinom2(link="sqrt")) # with log-link function, nbinom1 does not converge

AIC(m1,m2,m3,m4) # best is nbinom1(log)

summary(m1)

#Evaluate residuals
s=simulateResiduals(fittedModel=m4,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m4)

tab_model(m4)
plot_model(m4,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m4)

## visualize model

# Flow with CI and psite_genus
mydf <- ggpredict(m4, terms= c("CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlab("Streamflow (m3/sec)")+ylab("# of pseudocysts per fish")


### Impact of point-source pollution - PIMVIG MYX.SBAD only CI----

pimvig_count_myxsbad <- subset(pimvig_count,psite_spp.x == "MYXO.SBAD")


### Which link function is the best?
m1 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = pimvig_count_myxsbad,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = pimvig_count_myxsbad,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = pimvig_count_myxsbad,
              family = nbinom2(link="log")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = pimvig_count_myxsbad,
              family = nbinom2(link="sqrt")) # with log-link function, nbinom1 does not converge

AIC(m1,m2,m3,m4) # best is nbinom1(log)

summary(m3)

#Evaluate residuals
s=simulateResiduals(fittedModel=m3,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m3)

tab_model(m3)
plot_model(m3,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m3)

## visualize model

# Flow with CI and psite_genus
mydf <- ggpredict(m3, terms= c("CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlab("Streamflow (m3/sec)")+ylab("# of pseudocysts per fish")+ylim(0,75)


### Impact of point-source pollution - CARVEL MYX.CM only CI----

caarvel_presence_myxcm <- subset(carvel_presence,psite_spp.x == "MYX.CM")


### Which link function is the best?
m1 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = caarvel_presence_myxcm,
              family = binomial(link="logit")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = caarvel_presence_myxcm,
              family = binomial(link="probit")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = caarvel_presence_myxcm,
              family = binomial(link="inverse")) # with log-link function, nbinom1 does not converge


AIC(m1,m2,m3) # best is binomial(probit), the logit did not converge

summary(m2)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m2)

tab_model(m2)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus
mydf <- ggpredict(m2, terms= c("CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlab("Streamflow (m3/sec)")+ylab("# of pseudocysts per fish")+ylim(0,1)


### Impact of point-source pollution - CARVEL MYX.STWY only CI----
gamaff_presence_myxstwy <- subset(gamaff_presence,psite_spp.x == "MYX.STWY")


### Which link function is the best?
m1 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = gamaff_presence_myxstwy,
              family = binomial(link="logit")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = gamaff_presence_myxstwy,
              family = binomial(link="probit")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ 
                CI+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected/season),
              data = gamaff_presence_myxstwy,
              family = binomial(link="inverse")) # with log-link function, nbinom1 does not converge


AIC(m1,m2,m3) # best is binomial(probit)

summary(m2)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m2)

tab_model(m2)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus
mydf <- ggpredict(m2, terms= c("CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlab("Streamflow (m3/sec)")+ylab("# of pseudocysts per fish")+ylim(0,1)


### Holding on to, The effect of time on abundance data----

higherthanfive <- subset(myxo_count, psite_spp.x == "MYX.F"|
                           psite_spp.x == "MYX.G"|
                           psite_spp.x == "MYX.TAIL"|
                           psite_spp.x == "MYX.SP"|
                           psite_spp.x == "MYX.GO"|
                           psite_spp.x == "MYX.THEL"|
                           psite_spp.x == "MYXO.SBAD")

higherthanfive$psite_spp.x <- as.factor(higherthanfive$psite_spp.x)

m1 <- glmmTMB(psite_count ~ scale(YearCollected)*CI+
                offset(logTL_mm)+(1|season)+
                (1+scale(YearCollected)|psite_spp.x)+
                (1+CI|psite_spp.x),
              data = higherthanfive,
              family = nbinom1(link="log")) 

m2 <- glmmTMB(psite_count ~ scale(YearCollected)*CI+
                offset(logTL_mm)+(1|season)+
                (YearCollected|psite_spp.x)+
                (CI|psite_spp.x),
              data = higherthanfive,
              family = nbinom1(link="sqrt")) 

m3 <- glmmTMB(psite_count ~ scale(YearCollected)*CI+
                offset(logTL_mm)+(1|season)+
                (1+scale(YearCollected)|psite_spp.x)+
                (1+CI|psite_spp.x),
              data = higherthanfive,
              family = nbinom2(link="log")) 

m4 <- glmmTMB(psite_count ~ scale(YearCollected)*CI+
                offset(logTL_mm)+(1|season)+
                (1+scale(YearCollected)|psite_spp.x)+
                (1+CI|psite_spp.x),
              data = higherthanfive,
              family = nbinom2(link="sqrt")) 

AIC(m1,m2,m3,m4)
summary(m2)


## Carpiodes velifer - MYX.G
unique(carvel_count$psite_spp.x)

carvel_count_myxg <- subset(carvel_count,psite_spp.x == "MYX.G")

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = carvel_count_myxg,
              family = nbinom1(link="log")) 

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = carvel_count_myxg,
              family = nbinom2(link="log")) 


m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = carvel_count_myxg,
              family = nbinom1(link="sqrt")) 

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = carvel_count_myxg,
              family = nbinom2(link="sqrt")) 

summary(m2)
AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

#Check over dispersion
performance::check_overdispersion(m4)

#Evaluate model
tab_model(m4)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m2, terms= c("YearCollected[n=100]","CI")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/Carvel_MyxG_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/Carvel_MyxG_Year.pdf", width=150, height=120, dpi=1000, units = "mm")


## Carpiodes velifer - MYX.F
carvel_count_myxf <- subset(carvel_count,psite_spp.x == "MYX.F")

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = carvel_count_myxf,
              family = nbinom1(link="log")) 

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = carvel_count_myxf,
              family = nbinom2(link="log")) 

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = carvel_count_myxf,
              family = nbinom1(link="sqrt")) 

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = carvel_count_myxf,
              family = nbinom2(link="sqrt")) 

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

#Check over dispersion
performance::check_overdispersion(m2)

#Evaluate model
tab_model(m2)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m2, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/Carvel_MyxF_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/Carvel_MyxF_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

### Ictalurus punctatus
## ICtalurus puncatus - MYX.TAIL
ictpun_count_myxtail <- subset(ictpun_count,psite_spp.x == "MYX.TAIL")

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = ictpun_count_myxtail,
              family = nbinom1(link="log")) 

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = ictpun_count_myxtail,
              family = nbinom2(link="log")) 

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = ictpun_count_myxtail,
              family = nbinom1(link="sqrt")) 

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = ictpun_count_myxtail,
              family = nbinom2(link="sqrt")) 

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

#Check over dispersion
performance::check_overdispersion(m2)

#Evaluate model
tab_model(m2)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m2, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/ICTPUN_MyxTail_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/ICTPUN_MyxTail_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

### Notropis atherinoides
## Notropis atherinoides - MYX.SP
notath_count_myxsp <- subset(notath_count,psite_spp.x == "MYX.SP")

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = notath_count_myxsp,
              family = nbinom1(link="log")) 

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = notath_count_myxsp,
              family = nbinom2(link="log")) 

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = notath_count_myxsp,
              family = nbinom1(link="sqrt")) 

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = notath_count_myxsp,
              family = nbinom2(link="sqrt")) 

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

#Check over dispersion
performance::check_overdispersion(m1)

#Evaluate model
tab_model(m2)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m2, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/NOTATH_MyxSP_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/NOTATH_MyxSP_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

### Pimephales vigilax
## Pimephales vigilax - MYX.GO
pimvig_count_myxgo <- subset(pimvig_count,psite_spp.x == "MYX.GO")

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = pimvig_count_myxgo,
              family = nbinom1(link="log")) 

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = pimvig_count_myxgo,
              family = nbinom2(link="log")) 

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = pimvig_count_myxgo,
              family = nbinom1(link="sqrt")) 

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = pimvig_count_myxgo,
              family = nbinom2(link="sqrt")) 

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

#Check over dispersion
performance::check_overdispersion(m1)

#Evaluate model
tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m1, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/PIMVIG_MyxGO_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/PIMVIG_MyxGO_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")


### Pimephales vigilax
## Pimephales vigilax - MYX.THEL
pimvig_count_myxthel <- subset(pimvig_count,psite_spp.x == "MYX.THEL")

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = pimvig_count_myxthel,
              family = nbinom1(link="log")) 

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = pimvig_count_myxthel,
              family = nbinom2(link="log")) 

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = pimvig_count_myxthel,
              family = nbinom1(link="sqrt")) 

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = pimvig_count_myxthel,
              family = nbinom2(link="sqrt")) 

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m4,n=250)
s$scaledResiduals
plot(s)

#Check over dispersion
performance::check_overdispersion(m1)

#Evaluate model
tab_model(m4)
plot_model(m4,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m4, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/PIMVIG_MyxTHEL_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/PIMVIG_MyxTHEL_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

### Pimephales vigilax
## Pimephales vigilax - MYXO.SBAD
pimvig_count_myxosbad <- subset(pimvig_count,psite_spp.x == "MYXO.SBAD")

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = pimvig_count_myxosbad,
              family = nbinom1(link="log")) 

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = pimvig_count_myxosbad,
              family = nbinom2(link="log")) 

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = pimvig_count_myxosbad,
              family = nbinom1(link="sqrt")) 

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|CI/site)+
                (1|season),
              data = pimvig_count_myxosbad,
              family = nbinom2(link="sqrt")) 

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

#Check over dispersion
performance::check_overdispersion(m2)

#Evaluate model
tab_model(m2)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m2, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Figures/Year/PIMVIG_MyxSBAD_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Figures/Year/PIMVIG_MyxSBAD_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

#### Things to hold on to----

# Some data only spans between 1973 and 1994. So we will perform the analyses 
# only for that time period and for the control because I cannot assume that the 
# downstram impact sites are going to be represented by this data.
mixN_year <- subset(mixN_unitall,  
                    YearCollected > 1972 &
                      YearCollected < 1995 )
mixN_control <- subset(mixN_year, CI == "control")
mixN <- subset(mixN_control, Unit_full == "mg/l")
