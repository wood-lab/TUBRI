### This script is meant to analyse the impact of multiple stressors on myxozoa abundance
### Written by Daki Diaz-Morales
### Revised August 17, 2024



### Preparation of data frames----

##Install packages

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


## Import data

full_dataset <- read_csv("data/processed/Full_dataset_with_psite_life_history_info_2024.08.27.csv")

## Merge with stream flow data


# Create a yearflow column for later merging

full_dataset$MonthCollected <- as.numeric(full_dataset$MonthCollected)

full_dataset$season_flow <- ifelse(full_dataset$MonthCollected > 5 &
                                     full_dataset$MonthCollected < 12,                # condition
                                             'low',    # what if condition is TRUE
                                             'high')       # what if condition is FALSE

full_dataset$yearflow <- paste(as.character(full_dataset$YearCollected),as.character(full_dataset$season_flow),sep="_")

# Read physical data
physicalUSGS <- read_excel("data/geospatial/Physicaldata_USGS.xlsx")

# Make a subset only for stream flow in m3/sec
streamflow <- subset(physicalUSGS, Result_Characteristic=="Stream flow")
streamflow_ms <- subset(streamflow, Unit=="m3/sec")

# Create a yearflow column for later merging

streamflow_ms$MonthCollected <- as.numeric(streamflow_ms$MonthCollected)

streamflow_ms$season_flow <- ifelse(streamflow_ms$MonthCollected > 5 &
                                      streamflow_ms$MonthCollected < 12,                # condition
                                   'low',    # what if condition is TRUE
                                   'high')       # what if condition is FALSE

streamflow_ms$yearflow <- paste(as.character(streamflow_ms$YearCollected),as.character(streamflow_ms$season_flow),sep="_")

# Get rid of columns we do not need
streamflow_clean <- cbind.data.frame(streamflow_ms$yearflow,streamflow_ms$Result)

# Name the columns

colnames(streamflow_clean)[1]<-"yearflow"
colnames(streamflow_clean)[2]<-"Result"

# Summarize by taking the mean stream flow per 'yearflow'

streamflow_mean <- streamflow_clean %>% group_by(yearflow) %>% 
  summarise(meanFlow=mean(Result),
            .groups = 'drop')

Full_dataset_physical <-merge(full_dataset, streamflow_mean, by.x = "yearflow", by.y = "yearflow", all.x = TRUE)


## Make a subset for myxozoans

full_dataset_myxo <- subset(full_dataset, Parasite_taxonomic_group == "Myxozoa")

## Create column for presence and absence of myxozoans
full_dataset_myxo$psite_presence <- ifelse(full_dataset_myxo$psite_count > 0,                # condition
                                               1,    # what if condition is TRUE
                                               0)       # what if condition is FALSE

## Remove columns with myxozoans that aren't confirmed to be myxozoans and therefore are not identified
full_dataset_myxo <- subset(full_dataset_myxo, !is.na(Parasite_genus))

#Revise variable types
full_dataset_myxo$IndividualFishID <- as.factor(full_dataset_myxo$IndividualFishID)
full_dataset_myxo$CI <- as.factor(full_dataset_myxo$CI)
full_dataset_myxo$Fish_sp.x <- as.factor(full_dataset_myxo$Fish_sp.x)
full_dataset_myxo$Parasite_genus <- as.factor(full_dataset_myxo$Parasite_genus)


## Create subset per fish species

pimvig_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Pimephales vigilax")
gamaff_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Gambusia affinis")
ictpun_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Ictalurus punctatus")
notath_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Notropis atherinoides")
carvel_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Carpiodes velifer")
hybnuc_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Hybognathus nuchalis")



## Make a subset for data before the clean water act

full_dataset_myxo_bca <- subset(full_dataset_myxo, YearCollected < 1973)


### Before clean water act----

# Include only impacted sites (between latitudes 30.60-30.76)
full_dataset_myxo_bcalat <- subset(full_dataset_myxo_bca, CI == "impact")

# Make a subset for each parasite genera
full_dataset_myxo_bcamyxobolus <- subset(full_dataset_myxo_bcalat, Parasite_genus == "Myxobolus")
full_dataset_myxo_bcamyxidium <- subset(full_dataset_myxo_bcalat, Parasite_genus == "Myxidium")
full_dataset_myxo_bcaChloromyxum <- subset(full_dataset_myxo_bcalat, Parasite_genus == "Chloromyxum")
full_dataset_myxo_bcaHenneguya <- subset(full_dataset_myxo_bcalat, Parasite_genus == "Henneguya")
full_dataset_myxo_bcaThelohanellus <- subset(full_dataset_myxo_bcalat, Parasite_genus == "Thelohanellus")
full_dataset_myxo_bcaUnicauda <- subset(full_dataset_myxo_bcalat, Parasite_genus == "Unicauda")

# Make subsets only for CARVEL
full_dataset_myxo_bcamyxoboluscarvel <- subset(full_dataset_myxo_bcamyxobolus, Fish_sp.x == "Carpiodes velifer")
full_dataset_myxo_bcaChloromyxumcarvel <- subset(full_dataset_myxo_bcaChloromyxum, Fish_sp.x == "Carpiodes velifer")

## What is the impact of distance from pulp mill on the abundance of myxozoans?

#Set a theme for all your plots
apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# Error plot for prevalence for all myxozoan genera along latitude

lat_myxo <-ggerrorplot(full_dataset_myxo_bca, x = "Latitude", y = "psite_presence",
                            ggtheme = theme_bw(), color="Parasite_genus",rawdata=TRUE,
                            position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Prevalence of infection")+
  ylim(0,1)+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)

# Error plot for prevalence for each myxozoan genera along latitude

lat_myxobolus <-ggerrorplot(full_dataset_myxo_bcamyxobolus, x = "Latitude", y = "psite_presence",
                                ggtheme = theme_bw(), rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Prevalence of infection")+
  ylim(0,1)+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)
  
lat_myxidium <-ggerrorplot(full_dataset_myxo_bcamyxidium, x = "Latitude", y = "psite_presence",
                            ggtheme = theme_bw(), rawdata=TRUE,
                            position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Prevalence of infection")+
  ylim(0,1)+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)

lat_chloro <-ggerrorplot(full_dataset_myxo_bcaChloromyxumcarvel, x = "Latitude", y = "psite_presence",
                           ggtheme = theme_bw(), rawdata=TRUE,
                           position=position_dodge(0.5),width=0.00, size=0.3)+
  apatheme+
  facet_wrap("Fish_sp.x")+
  xlab("Latitude")+ylab("Prevalence of infection")+
  ylim(0,1)+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)

lat_henne <-ggerrorplot(full_dataset_myxo_bcaHenneguya, x = "Latitude", y = "psite_presence",
                         ggtheme = theme_bw(), rawdata=TRUE,
                         position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Prevalence of infection")+
  ylim(0,1)+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)

lat_thelo <-ggerrorplot(full_dataset_myxo_bcaThelohanellus, x = "Latitude", y = "psite_presence",
                        ggtheme = theme_bw(), rawdata=TRUE,
                        position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Prevalence of infection")+
  ylim(0,1)+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)

lat_uni <-ggerrorplot(full_dataset_myxo_bcaUnicauda, x = "Latitude", y = "psite_presence",
                        ggtheme = theme_bw(), rawdata=TRUE,
                        position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Prevalence of infection")+
  ylim(0,1)+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)


# Error plot for abundance

lat_myxobolus <-ggerrorplot(full_dataset_myxo_bcamyxobolus, x = "Latitude", y = "psite_count",
                            ggtheme = theme_bw(), rawdata=TRUE,
                            position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Parasite abundance")+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)

lat_myxidium <-ggerrorplot(full_dataset_myxo_bcamyxidium, x = "Latitude", y = "psite_count",
                           ggtheme = theme_bw(), rawdata=TRUE,
                           position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Prevalence of infection")+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)

lat_chloro <-ggerrorplot(full_dataset_myxo_bcaChloromyxum, x = "Latitude", y = "psite_count",
                         ggtheme = theme_bw(), rawdata=TRUE,
                         position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Parasite abundance")

lat_henne <-ggerrorplot(full_dataset_myxo_bcaHenneguya, x = "Latitude", y = "psite_count",
                        ggtheme = theme_bw(), rawdata=TRUE,
                        position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Prevalence of infection")+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)

lat_thelo <-ggerrorplot(full_dataset_myxo_bcaThelohanellus, x = "Latitude", y = "psite_count",
                        ggtheme = theme_bw(), rawdata=TRUE,
                        position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Prevalence of infection")+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)

lat_uni <-ggerrorplot(full_dataset_myxo_bcaUnicauda, x = "Latitude", y = "psite_count",
                      ggtheme = theme_bw(), rawdata=TRUE,
                      position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Prevalence of infection")+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)


# GLM for the probability of finding myxozoans along the river

glm_presence<-glmmTMB(psite_presence ~ Latitude+
                                      scaled_TL_mm+
                                      (1|MonthCollected),
                      data = full_dataset_myxo_bcamyxoboluscarvel,
                      family=binomial())

glm_presence<-glmmTMB(psite_presence ~ Latitude+
                        (1|MonthCollected),
                      data = full_dataset_myxo_bcaChloromyxumcarvel,
                      family=binomial()) #removed the term length because it gave issues

summary(glm_presence)

# Evaluate residuals
# Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_presence,n=250)
s$scaledResiduals
plot(s)


# With the plot()function
mydf <- ggpredict(glm_presence, c("Latitude")) 

apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

plot(mydf,show_data=TRUE)+
  labs(x = 'Latitude', y = 'Probability of myxozoan presence (%)',title=NULL)+
  apatheme+
  scale_x_reverse()

# GLM for the abundance of myxozoans along the river
glm_count<-glmmTMB(psite_count ~ Latitude+
                        scaled_TL_mm+
                        (1|MonthCollected),
                      data = full_dataset_myxo_bcamyxoboluscarvel,
                      family=nbinom1())

glm_count<-glmmTMB(psite_count ~ Latitude+
                     scaled_TL_mm+
                     (1|MonthCollected),
                   data = full_dataset_myxo_bcaChloromyxumcarvel,
                   family=nbinom1()) # not converging even with nbinom2 and removing fish length

summary(glm_count)

#Evaluate residuals
#Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_count,n=250)
s$scaledResiduals
plot(s)


#With the plot()function
mydf <- ggpredict(glm_count, c("Latitude")) 

apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

plot(mydf,rawdata=TRUE)+
  labs(x = 'Latitude', y = 'Myxozoan abundance',title=NULL)+
  apatheme+
  scale_x_reverse()

### Myxozoan prevalence across time----

#Set a theme for all your plots
apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# Plot myxozoan prevalence against time 

pimvig_prevalence <-ggerrorplot(pimvig_myxo, x = "YearCollected", y = "psite_presence",
                      ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                      position=position_dodge(0.5),width=0.00, size=0.3)+
                      facet_wrap("Parasite_genus")+
                      apatheme+ggtitle("Pimephales vigilax")+
                      xlab("Year")+ylab("Prevalence of infection")+
                                          ylim(0,1)

pimvig_prevalence

# Plot myxozoan prevalence against time - GAMAFF

gamaff_prevalence <-ggerrorplot(gamaff_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
                                facet_wrap("Parasite_genus")+
                                apatheme+ggtitle("Gambusia affinis")+
                                xlab("Year")+ylab("Prevalence of infection")+
                                ylim(0,1)

gamaff_prevalence

# Plot myxozoan prevalence against time - ICTPUN

ictpun_prevalence <-ggerrorplot(ictpun_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Ictalurus punctatus")+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

ictpun_prevalence


# Plot myxozoan prevalence against time - NOTATH

notath_prevalence <-ggerrorplot(notath_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Notropis atherinoides")+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

notath_prevalence

# Plot myxozoan prevalence against time - CARVEL

carvel_prevalence <-ggerrorplot(carvel_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Carpiodes velifer")+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

carvel_prevalence


# Plot myxozoan prevalence against time - HYBNUC

hybnuc_prevalence <-ggerrorplot(hybnuc_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Hybognathus nuchalis")+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

hybnuc_prevalence


# GLM for the probability of finding myxozoans
glm_presence <- glmer(psite_presence ~ YearCollected+CI*Parasite_genus+
                        (1|Fish_sp.x/CatalogNumber),
                  data = full_dataset_myxo,family=binomial())

glm_presence <- glmer(psite_presence ~ poly(Latitude,3)+YearCollected+
                        CI*Parasite_genus+scaled_TL_mm+
                      (1|Fish_sp.x/IndividualFishID)+
                      (1|MonthCollected),
                    data = full_dataset_myxo,family=binomial())


glm_presence <- glmmTMB(psite_presence ~ YearCollected+CI*Parasite_genus+scaled_TL_mm+
                        (1|Fish_sp.x/IndividualFishID)+
                        (1|MonthCollected),
                      data = full_dataset_myxo,family=binomial())



summary(glm_presence)

#Evaluate residuals
#Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_presence,n=250)
s$scaledResiduals
plot(s)


#With the plot()function
mydf <- ggpredict(glm_presence, c("Latitude[all]","Parasite_genus")) 

apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

plot(mydf,rawdata=TRUE)+
  labs(x = 'Latitude', y = 'Probability of myxozoan presence (%)',title=NULL)+
  apatheme+
  geom_vline(xintercept=30.76, linetype="dashed", color = "black", size=0.5)+
  scale_x_reverse()

### Myxozoan abundace across time----

## Make a subset for data before the clean water act

full_dataset_myxo_bca <- subset(full_dataset_myxo, YearCollected < 1973)


glm_presence <- glmmTMB(psite_presence ~ poly(Latitude,2)*Parasite_genus+YearCollected+scaled_TL_mm+
                          (1|Fish_sp.x/IndividualFishID)+
                          (1|MonthCollected),
                        data = full_dataset_myxo_bcalat,family=binomial())


#Set a theme for all your plots
apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# Plot myxozoan prevalence against time 

pimvig_prevalence <-ggerrorplot(pimvig_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Pimephales vigilax")+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

pimvig_prevalence

# Plot myxozoan prevalence against time - GAMAFF

gamaff_prevalence <-ggerrorplot(gamaff_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Gambusia affinis")+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

gamaff_prevalence

# Plot myxozoan prevalence against time - ICTPUN

ictpun_prevalence <-ggerrorplot(ictpun_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Ictalurus punctatus")+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

ictpun_prevalence


# Plot myxozoan prevalence against time - NOTATH

notath_prevalence <-ggerrorplot(notath_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Notropis atherinoides")+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

notath_prevalence

# Plot myxozoan prevalence against time - CARVEL

carvel_prevalence <-ggerrorplot(carvel_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Carpiodes velifer")+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

carvel_prevalence


# Plot myxozoan prevalence against time - HYBNUC

hybnuc_prevalence <-ggerrorplot(hybnuc_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Hybognathus nuchalis")+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)


hybnuc_prevalence


# GLM for the probability of finding myxozoans
glm_presence <- glmer(psite_presence ~ YearCollected+CI*Parasite_genus+
                        (1|Fish_sp.x/CatalogNumber),
                      data = full_dataset_myxo,family=binomial())

glm_presence <- glmer(psite_presence ~ poly(Latitude,3)+YearCollected+
                        CI*Parasite_genus+scaled_TL_mm+
                        (1|Fish_sp.x/IndividualFishID)+
                        (1|MonthCollected),
                      data = full_dataset_myxo,family=binomial())


glm_presence <- glmmTMB(psite_presence ~ poly(Latitude,3)+YearCollected+CI*Parasite_genus+scaled_TL_mm+
                          (1|Fish_sp.x/IndividualFishID)+
                          (1|MonthCollected),
                        data = full_dataset_myxo,family=binomial())



summary(glm_presence)

#Evaluate residuals
#Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_presence,n=250)
s$scaledResiduals
plot(s)


#With the plot()function
mydf <- ggpredict(glm_presence, c("Latitude [all]","CI","Parasite_genus")) 

apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

plot(mydf,rawdata=TRUE)+
  labs(x = 'Latitude', y = 'Probability of myxozoan presence (%)',title=NULL)+
  apatheme+
  geom_vline(xintercept=30.76, linetype="dashed", color = "black", size=0.5)

### Streamflow----

## Import data

physicalUSGS <- read_excel("data/geospatial/Physicaldata_USGS.xlsx")

View(physicalUSGS)


streamflow <- subset(physicalUSGS, Result_Characteristic=="Stream flow")
streamflow_ms <- subset(streamflow, Unit=="m3/sec")

# Plot stream flow against time 

streamflow_plot <-ggerrorplot(streamflow_ms, x = "Year", y = "Result",
                                ggtheme = theme_bw(),rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  ggtitle("Streamflow")+
  xlab("Year")+ylab("Streamflow (m3/sec)")

streamflow_plot

ggplot(streamflow_ms, aes(x= as.factor(Year),
                  y=Result))+
  geom_point()+apatheme+
  ggtitle("Streamflow per year")+
  xlab("Year")+ylab("Streamflow (m3/sec)")

ggplot(streamflow_ms, aes(x= as.factor(Month),
                           y=Result))+
  geom_point()+apatheme+
  ggtitle("Streamflow per month")+
  xlab("Month")+ylab("Streamflow (m3/sec)")


### Water Temperature USGS----

## Import data

physicalUSGS <- read_excel("data/geospatial/Physicaldata_USGS.xlsx")

View(physicalUSGS)


wtemp <- subset(physicalUSGS, Result_Characteristic=="Temperature, water")

# Plot stream flow against time 

wtemp_plot <-ggerrorplot(wtemp, x = "YearCollected", y = "Result",
                              ggtheme = theme_bw(),rawdata=TRUE,
                              position=position_dodge(0.5),width=0.00, size=0.3)+
  ggtitle("Water temperature")+
  xlab("Year")+ylab("Temperature (C)")

wtemp_plot

# Water temperature per month

apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())


ggplot(wtemp, aes(x= as.factor(Month),
                            y=Result))+
  geom_point()+apatheme+
  ggtitle("Water temperature per month")+
  xlab("Month")+ylab("Temperature (degC)")+
  geom_hline(yintercept=20, linetype="dashed", color = "black", size=0.5)




