


### Load packages----

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

### Load parasite data and build data frames----

## Import parasite dataset

full_dataset <- read_csv("data/processed/Full_dataset_physical_2024.09.23.csv")

## Make a subset for myxozoans

#full_dataset_myxo <- subset(Full_dataset_elemN, Parasite_taxonomic_group == "Myxozoa")
full_dataset_myxo <- subset(full_dataset, Parasite_taxonomic_group == "Myxozoa")

## Create column for presence and absegroup_by()## Create column for presence and absence of myxozoans
full_dataset_myxo$psite_presence <- ifelse(full_dataset_myxo$psite_count > 0,                # condition
                                           1,    # what if condition is TRUE
                                           0)       # what if condition is FALSE

## Remove columns with myxozoans that aren't confirmed to be myxozoans and therefore are not identified
full_dataset_myxo <- subset(full_dataset_myxo, !is.na(Parasite_genus))
full_dataset_myxo <- subset(full_dataset_myxo, psite_spp.x!="MYX.GM") # According to Stephen, this is likely a contamination and also it is in the gall bladder so does not count for abundace
full_dataset_myxo <- subset(full_dataset_myxo, psite_spp.x!="MYX.GEOM") # According to Stephen, this is likely a contamination and also it is in the gall bladder so does not count for abundace

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


## Add sites for random structure of model

# Combine lat and long to create a unique identifier for site
full_dataset_myxo$lat_long <- paste(full_dataset_myxo$Latitude,
                                    full_dataset_myxo$Longitude,
                                    sep="_")

full_dataset_myxo$lat_long <- as.factor(full_dataset_myxo$lat_long)

full_dataset_myxo_control <- subset(full_dataset_myxo, CI=="control")
full_dataset_myxo_impact <- subset(full_dataset_myxo, CI=="impact")

# What are the unique combinations of lat and long for control and for impact
top.xy_control <- unique(full_dataset_myxo_control[c("lat_long")])
top.xy_control

ggplot(top.xy_control, aes(x= Longitude,y=Latitude))+
  geom_point()+apatheme

top.xy_impact <- unique(full_dataset_myxo_impact[c("lat_long")])
top.xy_impact

ggplot(top.xy_impact, aes(x= Longitude,y=Latitude))+
  geom_point()+apatheme


# Set longitude as factor and use it to establish site categories (four sites in total)
full_dataset_myxo <- full_dataset_myxo %>% 
  mutate(site = case_when(  lat_long == '30.76805_-89.83083' ~ "CA",
                            lat_long == '30.76639_-89.83195' ~ "CA",
                            lat_long == '30.7675_-89.83028' ~ "CA",
                            lat_long == '30.76528_-89.83222' ~ "CA",
                            lat_long == '30.76222_-89.83111' ~ "CA",
                            lat_long == '30.77611_-89.82777' ~ "CB",
                            lat_long == '30.77861_-89.82972' ~ "CB",
                            lat_long == '30.77611_-89.82722' ~ "CB",
                            lat_long == '30.77583_-89.82722' ~ "CB",
                            lat_long == '30.77472_-89.82806' ~ "CB",
                            lat_long == '30.7825_-89.82806'  ~ "CB",
                            lat_long == '30.7825_-89.82027' ~ "CC",
                            lat_long == '30.78472_-89.81944' ~ "CC",
                            lat_long == '30.705_-89.84611' ~ "IA",
                            lat_long == '30.70222_-89.84417' ~ "IA",
                            lat_long == '30.70389_-89.84472' ~ "IA",
                            lat_long == '30.73694_-89.82694' ~ "IB",
                            lat_long == '30.74195_-89.82528' ~ "IB",
                            lat_long == '30.74472_-89.82528' ~ "IB",
                            lat_long == '30.74055_-89.82694' ~ "IB",
                            lat_long == '30.74_-89.82777' ~ "IB",
                            lat_long == '30.74083_-89.82611' ~ "IB",
                            lat_long == '30.75361_-89.82639' ~ "IC",
                            lat_long == '30.75667_-89.82611' ~ "IC",
                            lat_long == '30.75417_-89.82694' ~ "IC",
                            lat_long == '30.75444_-89.82722' ~ "IC"))

full_dataset_myxo$site <- as.factor(full_dataset_myxo$site)

# check with 30.76528_-89.83222 & 30.76222_-89.83111
check <- subset(full_dataset_myxo, lat_long =='30.76528_-89.83222'|
                lat_long=='30.76222_-89.83111')

check <- subset(full_dataset_myxo, lat_long=='30.76222_-89.83111')


#Revise variable types
full_dataset_myxo$IndividualFishID <- as.factor(full_dataset_myxo$IndividualFishID)
full_dataset_myxo$CI <- as.factor(full_dataset_myxo$CI)
full_dataset_myxo$Fish_sp.x <- as.factor(full_dataset_myxo$Fish_sp.x)
full_dataset_myxo$Parasite_genus <- as.factor(full_dataset_myxo$Parasite_genus)
full_dataset_myxo$before_after <- as.factor(full_dataset_myxo$before_after)
full_dataset_myxo$season <- as.factor(full_dataset_myxo$season)

## Revise reference levels

#full_dataset_myxo$before_after <- relevel(full_dataset_myxo$before_after, ref = "before")
#full_dataset_myxo$CI <- relevel(full_dataset_myxo$before_after, ref = "control")


## Create subset per fish species

pimvig_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Pimephales vigilax")
gamaff_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Gambusia affinis")
ictpun_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Ictalurus punctatus")
notath_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Notropis atherinoides")
carvel_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Carpiodes velifer")
hybnuc_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Hybognathus nuchalis")







## Restrict to year and control

# Restrict to the years that we are going to work with
full_dataset_years <- subset(full_dataset, 
                             YearCollected > 1972 &
                               YearCollected < 1995 )

# Restrict to the years that we are going to work with
full_dataset_yearscontrol <- subset(full_dataset_years, 
                                    CI == "control" )

# Merge with parasite data

#Full_dataset_elements <- merge(full_dataset, mod, by.x = "YearCollected", by.y = "YearCollected", all.x = TRUE)
Full_dataset_elements <- merge(full_dataset_yearscontrol, mod, by.x = "YearCollected", by.y = "YearCollected", all.x = TRUE)

## Make a subset for myxozoans

#full_dataset_myxo <- subset(Full_dataset_elemN, Parasite_taxonomic_group == "Myxozoa")
full_dataset_myxo <- subset(Full_dataset_elements, Parasite_taxonomic_group == "Myxozoa")

## Create column for presence and absence of myxozoans
full_dataset_myxo$psite_presence <- ifelse(full_dataset_myxo$psite_count > 0,                # condition
                                           1,    # what if condition is TRUE
                                           0)       # what if condition is FALSE

## Remove columns with myxozoans that aren't confirmed to be myxozoans and therefore are not identified
full_dataset_myxo <- subset(full_dataset_myxo, !is.na(Parasite_genus))
full_dataset_myxo <- subset(full_dataset_myxo, psite_spp.x!="MYX.GM") # According to Stephen, this is likely a contamination and also it is in the gall bladder so does not count for abundace

#Revise variable types
full_dataset_myxo$IndividualFishID <- as.factor(full_dataset_myxo$IndividualFishID)
full_dataset_myxo$CI <- as.factor(full_dataset_myxo$CI)
full_dataset_myxo$Fish_sp.x <- as.factor(full_dataset_myxo$Fish_sp.x)
full_dataset_myxo$Parasite_genus <- as.factor(full_dataset_myxo$Parasite_genus)
full_dataset_myxo$before_after <- as.factor(full_dataset_myxo$before_after)
full_dataset_myxo$MonthCollected <- as.factor(full_dataset_myxo$MonthCollected)
full_dataset_myxo$season <- as.factor(full_dataset_myxo$season)

## Revise reference levels

#full_dataset_myxo$before_after <- relevel(full_dataset_myxo$before_after, ref = "before")
#full_dataset_myxo$CI <- relevel(full_dataset_myxo$before_after, ref = "control")


## Create subset per fish species

pimvig_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Pimephales vigilax")
gamaff_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Gambusia affinis")
ictpun_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Ictalurus punctatus")
notath_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Notropis atherinoides")
carvel_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Carpiodes velifer")
hybnuc_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Hybognathus nuchalis")

## Create subset per parasite genus with river physical data

myxobolus_physical <- subset(full_dataset_myxo, Parasite_genus == "Myxobolus")
chloromyxum_physical <- subset(full_dataset_myxo, Parasite_genus == "Chloromyxum")
henneguya_physical <- subset(full_dataset_myxo, Parasite_genus == "Henneguya")
myxidium_physical<- subset(full_dataset_myxo, Parasite_genus == "Myxidium")
unicauda_physical <- subset(full_dataset_myxo, Parasite_genus == "Unicauda")
thelohanellus_physical <- subset(full_dataset_myxo, Parasite_genus == "Thelohanellus")

## Crease subset for myxozoans that can be counted

myxo_count <- subset(full_dataset_myxo, Parasite_genus == "Myxobolus"|Parasite_genus == "Henneguya"|Parasite_genus == "Unicauda"|Parasite_genus == "Thelohanellus")


# Scale variables

myxo_count$logTL_mm <- log(myxo_count$TotalLength_mm)
myxo_count$scaledmean_streamflow <- scale(myxo_count$mean_streamflow)
myxo_count$scaledmean_nitrogen <- scale(myxo_count$mean_nitrogen)
myxo_count$scaledmean_temperature <- scale(myxo_count$mean_temperature)
myxo_count$scaledYear <- scale(myxo_count$YearCollected)

## Create subset of myxo_count per fish species

pimvig_count <- subset(myxo_count, Fish_sp.x == "Pimephales vigilax")
gamaff_count <- subset(myxo_count, Fish_sp.x == "Gambusia affinis")
ictpun_count <- subset(myxo_count, Fish_sp.x == "Ictalurus punctatus")
notath_count <- subset(myxo_count, Fish_sp.x == "Notropis atherinoides")
carvel_count <- subset(myxo_count, Fish_sp.x == "Carpiodes velifer")
hybnuc_count <- subset(myxo_count, Fish_sp.x == "Hybognathus nuchalis")

# Set a theme for all your plots
apatheme= theme_bw(base_size = 14,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
