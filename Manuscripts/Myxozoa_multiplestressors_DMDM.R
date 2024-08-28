### This script is meant to analyse the impact of multiple stressors on myxozoa abundance
### Written by Daki Diaz-Morales
### Revised August 17, 2024



### Preparation of data frames----

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


## Import data

full_dataset <- read_csv("data/processed/Full_dataset_with_psite_life_history_info_2024.08.27.csv")

## Merge with physical data

# Stream flow
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

# Merge with parasite data

Full_dataset_streamflow <- merge(full_dataset, streamflow_mean, by.x = "yearflow", by.y = "yearflow", all.x = TRUE)

# Adding water temperature data
# Create a yearseason column in our data for later merging

Full_dataset_streamflow$MonthCollected <- as.numeric(Full_dataset_streamflow$MonthCollected)

Full_dataset_streamflow$season_temp <- ifelse(Full_dataset_streamflow$MonthCollected > 4 &
                                                Full_dataset_streamflow$MonthCollected < 10,                # condition
                                   'summer',    # what if condition is TRUE
                                   'winter')       # what if condition is FALSE

Full_dataset_streamflow$yearseason <- paste(as.character(Full_dataset_streamflow$YearCollected),as.character(Full_dataset_streamflow$season_temp),sep="_")

# Make a subset of only water temperature
wtemp <- subset(physicalUSGS, Result_Characteristic=="Temperature, water")

# Create a yearseason column in USGS data for later merging

wtemp$MonthCollected <- as.numeric(wtemp$MonthCollected)

wtemp$season_temp <- ifelse(wtemp$MonthCollected > 4 &
                              wtemp$MonthCollected < 10,                # condition
                                    'summer',    # what if condition is TRUE
                                    'winter')       # what if condition is FALSE

wtemp$yearseason <- paste(as.character(wtemp$YearCollected),as.character(wtemp$season_temp),sep="_")

# Get rid of columns we do not need
wtemp_clean <- cbind.data.frame(wtemp$yearseason,wtemp$Result)

# Name the columns

colnames(wtemp_clean)[1]<-"yearseason"
colnames(wtemp_clean)[2]<-"Result"

# Summarize by taking the mean stream flow per 'yearseason'

wtemp_mean <- wtemp_clean %>% group_by(yearseason) %>% 
  summarise(meanTemp=mean(Result),
            .groups = 'drop')

# Merge with parasite data

Full_dataset_physical <- merge(Full_dataset_streamflow, wtemp_mean, by.x = "yearseason", by.y = "yearseason", all.x = TRUE)

## Make a subset for myxozoans

full_dataset_myxo <- subset(Full_dataset_physical, Parasite_taxonomic_group == "Myxozoa")

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

## Create subset per paraasite genus with river physical data

myxobolus_physical <- subset(full_dataset_myxo, Parasite_genus == "Myxobolus")
chloromyxum_physical <- subset(full_dataset_myxo, Parasite_genus == "Chloromyxum")
henneguya_physical <- subset(full_dataset_myxo, Parasite_genus == "Henneguya")
myxidium_physical<- subset(full_dataset_myxo, Parasite_genus == "Myxidium")
unicauda_physical <- subset(full_dataset_myxo, Parasite_genus == "Unicauda")
thelohanellus_physical <- subset(full_dataset_myxo, Parasite_genus == "Thelohanellus")


### Before clean water act, impact of distance from pulp mill on prevalence----

## What is the impact of distance from pulp mill on the abundance of myxozoans?


## Make a subset for data before the clean water act
full_dataset_myxo_bca <- subset(full_dataset_myxo, YearCollected < 1973)


## General plots
# Error plot for prevalence for all myxozoan genera along latitude

lat_myxo <-ggerrorplot(full_dataset_myxo_bca, x = "Latitude", y = "psite_presence",
                       ggtheme = theme_bw(), color="Parasite_genus",rawdata=TRUE,
                       position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+
  xlab("Latitude")+ylab("Prevalence of infection")+
  ylim(0,1)+
  geom_vline(xintercept=30.76222, linetype="dashed", color = "black", size=0.5)

## Test the effect of distance from the pulp mill (with latitude as proxy for distance) on the prevalence of myxozoanss

# GLM for the probability of finding myxozoans along the river

glm_presence<-glmmTMB(psite_presence ~ Latitude*CI*Parasite_genus+
                        Latitude*CI*Fish_sp.x+
                        scaled_TL_mm+
                        (1|MonthCollected),
                      data = full_dataset_myxo_bca,
                      family=binomial())

      # This model failed to converge. Let's make it a bit simpler by running the model per fish species

# GLM for the probability of finding myxozoans in CARVEL along the river
glm_presence<-glmmTMB(psite_presence ~ Latitude*CI*Parasite_genus+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = carvel_myxo,
                      family=binomial())


summary(glm_presence) 
      # Model converged, diagnostics are great, no significances at all


# GLM for the probability of finding myxozoans in PIMVIG along the river
glm_presence<-glmmTMB(psite_presence ~ Latitude*CI*Parasite_genus+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = pimvig_myxo,
                      family=binomial())


summary(glm_presence) 
      # Model converged, diagnostics are great, only fish size was significant (p < 0.01, estimate = -4.182e-02)

# GLM for the probability of finding myxozoans in GAMAFF along the river
glm_presence<-glmmTMB(psite_presence ~ Latitude*CI*Parasite_genus+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = gamaff_myxo,
                      family=binomial())


summary(glm_presence) 
      # Model did not converge, however not much to see here since prevalence of infection was really low overall


# GLM for the probability of finding myxozoans in ICTPUN along the river, in this case no interaction with parasite genus because there is only one genus for this fish
glm_presence<-glmmTMB(psite_presence ~ Latitude*CI+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = ictpun_myxo,
                      family=binomial())

summary(glm_presence) 
      # Model did not converge,  however not much to see here since prevalence of infection was really low overall


# GLM for the probability of finding myxozoans in NOTATH along the river, , in this case no interaction with parasite genus because there is only one genus for this fish
glm_presence<-glmmTMB(psite_presence ~ Latitude*CI+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = notath_myxo,
                      family=binomial())


summary(glm_presence) 
      # Model converged, diagnostics are great, no significances


# GLM for the probability of finding myxozoans in HYBNUC along the river, , in this case no interaction with parasite genus because there is only one genus for this fish
glm_presence<-glmmTMB(psite_presence ~ Latitude*CI*Parasite_genus+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = hybnuc_myxo,
                      family=binomial())


summary(glm_presence) 
      # Model did not converge, however not much to see here since prevalence of infection was really low


# Evaluate residuals
# Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_presence,n=250)
s$scaledResiduals
plot(s)


# With the plot()function
mydf <- ggpredict(glm_presence, c("Latitude[all]","Parasite_genus")) 

apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

plot(mydf,show_data=TRUE)+
  labs(x = 'Latitude', y = 'Probability of myxozoan presence (%)',title=NULL)+
  apatheme

# In conclusion, although there was a trend for Myxobolous and Chloromyxum in CARVEl to decrease in prevalence with distance downstream of the pulpmill, there was not a significant effect of distance from pulp mill on the probability of finding myxozoans in fish. 

### Impact of stream flow on parasite presence----
#Set a theme for all your plots
apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# Plot myxozoan prevalence against streamflow 

myxobolus_flow_plot <-ggerrorplot(myxobolus_physical, x = "meanFlow", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+ggtitle("Myxobolus")+
  xlab("Stream flow (m3/sec)")+ylab("Prevalence of infection")

myxobolus_flow_plot

chloromyxum_flow_plot <-ggerrorplot(chloromyxum_physical, x = "meanFlow", y = "psite_presence",
                                  ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                  position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+ggtitle("Chloromyxum")+
  xlab("Stream flow (m3/sec)")+ylab("Prevalence of infection")

chloromyxum_flow_plot


henneguya_flow_plot <-ggerrorplot(henneguya_physical, x = "meanFlow", y = "psite_presence",
                                    ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                    position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+ggtitle("Henneguya")+
  xlab("Stream flow (m3/sec)")+ylab("Prevalence of infection")

henneguya_flow_plot


myxidium_flow_plot <-ggerrorplot(myxidium_physical, x = "meanFlow", y = "psite_presence",
                                  ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                  position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+ggtitle("Myxidium")+
  xlab("Stream flow (m3/sec)")+ylab("Prevalence of infection")

myxidium_flow_plot


thelohanellus_flow_plot <-ggerrorplot(thelohanellus_physical, x = "meanFlow", y = "psite_presence",
                                 ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                 position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+ggtitle("Thelohanellus")+
  xlab("Stream flow (m3/sec)")+ylab("Prevalence of infection")

thelohanellus_flow_plot

Unicauda_flow_plot <-ggerrorplot(unicauda_physical, x = "meanFlow", y = "psite_presence",
                                      ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                      position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Fish_sp.x")+
  apatheme+ggtitle("Unicauda")+
  xlab("Stream flow (m3/sec)")+ylab("Prevalence of infection")

Unicauda_flow_plot


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


glm_presence <- glmmTMB(psite_presence ~ YearCollected+
                          CI*Fish_sp.x+
                          meanFlow*Fish_sp.x+scaled_TL_mm+
                        (1|MonthCollected),
                      data = myxobolus_physical,family=binomial())



summary(glm_presence)

#Evaluate residuals
#Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_presence,n=250)
s$scaledResiduals
plot(s)


#With the plot()function
plot_model(glm_presence)+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)


mydf <- ggpredict(glm_presence, c("meanFlow[all]","CI","Fish_sp.x")) 
mydf <- ggpredict(glm_presence, c("YearCollected[all]","CI","Fish_sp.x")) 

apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

plot(mydf,rawdata=TRUE)+
  labs(x = 'Stream flow (m3/sec)', y = 'Probability of myxozoan presence (%)',title=NULL)+
  apatheme

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

### Streamflow visualization----

# Plot stream flow as means throught time

streamflow_plot <-ggerrorplot(streamflow_ms, x = "Year", y = "Result",
                                ggtheme = theme_bw(),rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  ggtitle("Streamflow")+
  xlab("Year")+ylab("Streamflow (m3/sec)")

streamflow_plot

# Plot stream flow as raw data per year

ggplot(streamflow_ms, aes(x= as.factor(Year),
                  y=Result))+
  geom_point()+apatheme+
  ggtitle("Streamflow per year")+
  xlab("Year")+ylab("Streamflow (m3/sec)")

# Plot stream flow as raw data per month

ggplot(streamflow_ms, aes(x= as.factor(Month),
                           y=Result))+
  geom_point()+apatheme+
  ggtitle("Streamflow per month")+
  xlab("Month")+ylab("Streamflow (m3/sec)")


### Water Temperature USGS visualization----
# Plot stream flow as means throughout time 
Full_dataset_physical
wtemp_plot <-ggerrorplot(wtemp, x = "YearCollected", y = "Result",
                              ggtheme = theme_bw(),rawdata=TRUE,
                              position=position_dodge(0.5),width=0.00, size=0.3)+
  ggtitle("Water temperature")+
  xlab("Year")+ylab("Temperature (C)")

wtemp_plot

# Water temperature per month

ggplot(wtemp, aes(x= as.factor(MonthCollected),
                            y=Result))+
  geom_point()+apatheme+
  ggtitle("Water temperature per month")+
  xlab("Month")+ylab("Temperature (degC)")+
  geom_hline(yintercept=20, linetype="dashed", color = "black", size=0.5)


# USGS data merged with our data, to double check
# Water temperature per month

ggplot(Full_dataset_physical, aes(x= as.factor(season_temp),
                  y=meanTemp))+
  geom_point()+apatheme+
  ggtitle("Water temperature per month")+
  xlab("Month")+ylab("Temperature (degC)")+
  geom_hline(yintercept=20, linetype="dashed", color = "black", size=0.5)


