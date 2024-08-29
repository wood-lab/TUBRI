### This script is meant to analyse the impact of multiple stressors on myxozoa abundance
### Written by Daki Diaz-Morales
### Revised August 17, 2024



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

###Prepare data frames-----
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

#Full_dataset_streamflow <- Full_dataset_streamflow %>% 
 # mutate(season_temp = case_when(MonthCollected >= 6 &
  #                                 MonthCollected <= 8 ~ "summer", # both tests: group A
   #                              MonthCollected >= 9 &
    #                               MonthCollected <= 11 ~ "fall",
     #                            MonthCollected >= 1 &
      #                             MonthCollected <= 2 
       #                          ~ "winter",
        #                         MonthCollected >= 3 &
         #                          MonthCollected <= 5 ~ "fall",
          #                       MonthCollected == 12 ~ "winter"
  #))


Full_dataset_streamflow$season_temp <- ifelse(Full_dataset_streamflow$MonthCollected > 4 &
                                                Full_dataset_streamflow$MonthCollected < 10,                # condition
                                   'summer',    # what if condition is TRUE
                                   'winter')       # what if condition is FALSE

Full_dataset_streamflow$yearseason <- paste(as.character(Full_dataset_streamflow$YearCollected),as.character(Full_dataset_streamflow$season_temp),sep="_")

# Make a subset of only water temperature
wtemp <- subset(physicalUSGS, Result_Characteristic=="Temperature, water")

# Create a yearseason column in USGS data for later merging

wtemp$MonthCollected <- as.numeric(wtemp$MonthCollected)


#wtemp <- wtemp %>% 
 # mutate(season_temp = case_when(MonthCollected >= 6 &
  #                           MonthCollected <= 8 ~ "summer", # both tests: group A
   #                        MonthCollected >= 9 &
    #                         MonthCollected <= 11 ~ "fall",
     #                      MonthCollected >= 1 &
      #                       MonthCollected <= 2 
       #                       ~ "winter",
        #                   MonthCollected >= 3 &
         #                    MonthCollected <= 5 ~ "fall",
          #                 MonthCollected == 12 ~ "winter"
  #))

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

# Summarize by taking the mean temperature 'yearseason'

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
full_dataset_myxo$before_after <- as.factor(full_dataset_myxo$before_after)
full_dataset_myxo$season_temp <- as.factor(full_dataset_myxo$season_temp)
full_dataset_myxo$MonthCollected <- as.factor(full_dataset_myxo$MonthCollected)


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

## Create subset of myxo_count per fish species

pimvig_count <- subset(myxo_count, Fish_sp.x == "Pimephales vigilax")
gamaff_count <- subset(myxo_count, Fish_sp.x == "Gambusia affinis")
ictpun_count <- subset(myxo_count, Fish_sp.x == "Ictalurus punctatus")
notath_count <- subset(myxo_count, Fish_sp.x == "Notropis atherinoides")
carvel_count <- subset(myxo_count, Fish_sp.x == "Carpiodes velifer")
hybnuc_count <- subset(myxo_count, Fish_sp.x == "Hybognathus nuchalis")

## Do a quick check for collinearity
# Get only columns for which we want to check collinearity -scaled
colfull_dataset_myxo <- cbind.data.frame(scale(full_dataset_myxo$meanFlow),
                                         scale(full_dataset_myxo$meanTemp),
                                         scale(full_dataset_myxo$YearCollected),
                                         full_dataset_myxo$scaled_TL_mm)


cor(colfull_dataset_myxo) 

# Get only columns for which we want to check collinearity - unscaled
colfull_dataset_myxo <- cbind.data.frame(full_dataset_myxo$meanFlow,
                                         full_dataset_myxo$meanTemp,
                                         full_dataset_myxo$YearCollected,
                                         full_dataset_myxo$scaled_TL_mm)


cor(colfull_dataset_myxo) 


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
                      family=binomial()) # This model failed to converge. Let's make it a bit simpler by running the model per fish species

# GLM for the probability of finding myxozoans in CARVEL along the river
glm_presence<-glmmTMB(psite_presence ~ Latitude*CI*Parasite_genus+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = carvel_myxo,
                      family=binomial())


summary(glm_presence) # Model converged, diagnostics are great, no significances at all


# GLM for the probability of finding myxozoans in PIMVIG along the river
glm_presence<-glmmTMB(psite_presence ~ Latitude*CI*Parasite_genus+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = pimvig_myxo,
                      family=binomial())


summary(glm_presence) # Model converged, diagnostics are great, only fish size was significant (p < 0.01, estimate = -4.182e-02)

# GLM for the probability of finding myxozoans in GAMAFF along the river
glm_presence<-glmmTMB(psite_presence ~ Latitude*CI*Parasite_genus+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = gamaff_myxo,
                      family=binomial())


summary(glm_presence) # Model did not converge, however not much to see here since prevalence of infection was really low overall


# GLM for the probability of finding myxozoans in ICTPUN along the river, in this case no interaction with parasite genus because there is only one genus for this fish
glm_presence<-glmmTMB(psite_presence ~ Latitude*CI+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = ictpun_myxo,
                      family=binomial())

summary(glm_presence)  # Model did not converge,  however not much to see here since prevalence of infection was really low overall


# GLM for the probability of finding myxozoans in NOTATH along the river, , in this case no interaction with parasite genus because there is only one genus for this fish
glm_presence<-glmmTMB(psite_presence ~ Latitude*CI+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = notath_myxo,
                      family=binomial())


summary(glm_presence) # Model converged, diagnostics are great, no significances


# GLM for the probability of finding myxozoans in HYBNUC along the river, in this case no interaction with parasite genus because there is only one genus for this fish
glm_presence<-glmmTMB(psite_presence ~ Latitude*CI*Parasite_genus+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = hybnuc_myxo,
                      family=binomial())


summary(glm_presence) # Model did not converge, however not much to see here since prevalence of infection was really low


# Evaluate residuals
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

### Myxozoan prevalence across time----

#Set a theme for all your plots
apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# Plot myxozoan prevalence against time - PIMVIG

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


# GLM for the probability of finding myxozoans taking all fish in consideration

glm_presence <- glmmTMB(psite_presence ~ meanTemp*Parasite_genus*Fish_sp.x+
                          CI*Parasite_genus*Fish_sp.x+
                          meanFlow*Parasite_genus*Fish_sp.x+
                        meanTemp*CI*meanFlow+
                        scaled_TL_mm+
                          (1|IndividualFishID)+
                        (1|MonthCollected),
                      data = full_dataset_myxo,family=binomial())

glm_presence <- glmer(psite_presence ~ meanTemp*Parasite_genus*Fish_sp.x+
                          CI*Parasite_genus*Fish_sp.x+
                          meanFlow*Parasite_genus*Fish_sp.x+
                          meanTemp*CI*meanFlow+
                          scaled_TL_mm+
                        (1|IndividualFishID)+
                          (1|MonthCollected),
                        data = full_dataset_myxo,family=binomial())


# The last two models did not converge. Let's make it simpler by rerunning the models per fish species


## GLM for the probability of finding myxozoans for CARVEL
carvel_myxo$Parasite_genus <- relevel(carvel_myxo$Parasite_genus, ref = "Myxobolus")
carvel_myxo$before_after <- relevel(carvel_myxo$before_after, ref = "before")

glm_presence <- glmer(psite_presence ~ scale(YearCollected)*
                        scale(meanTemp)*CI*scale(meanFlow)*Parasite_genus+
                        scaled_TL_mm+
                        (1|IndividualFishID)+
                        (1|MonthCollected),
                      data = carvel_myxo,family=binomial()) #did not converge


glm_presence <- glmmTMB(psite_presence ~ scale(YearCollected)*
                        scale(meanTemp)*CI*scale(meanFlow)*Parasite_genus+
                        scaled_TL_mm+
                      (1|IndividualFishID)+
                      (1|MonthCollected),
                      data = carvel_myxo,family=binomial()) #did not converge


glm_presence <- glmer(psite_presence ~ scale(YearCollected)+
                        scale(meanTemp)*CI*scale(meanFlow)+
                        CI*before_after+
                        scaled_TL_mm+
                        scale(meanTemp)*Parasite_genus+
                        scale(meanFlow)*Parasite_genus+
                        CI*Parasite_genus+
                        (1|IndividualFishID)+
                        (1|MonthCollected),
                      data = carvel_myxo,family=binomial()) #singular fit, diagnostics are OK


summary(glm_presence)


#Evaluate residuals
#Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_presence,n=250)
s$scaledResiduals
plot(s)


#With the plot()function
plot_model(glm_presence,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

# Flow with CI and psite_genus

mydf <- ggpredict(glm_presence, c("meanFlow[n=200]","CI","Parasite_genus")) 

    plot(mydf,rawdata=TRUE,jitter=0.05)+
      labs(x = 'Stream flow (m3/sec)', y = 'Probability of myxozoan presence (%)',title=NULL)+
      apatheme

# Temperature with CI and psite_genus
mydf <- ggpredict(glm_presence, c("meanTemp[n=200]","CI","Parasite_genus")) 

    plot(mydf,show_data = TRUE, show_ci = TRUE, jitter=0.05)+
      labs(x = 'Temperature (C)', y = 'Probability of myxozoan presence (%)',title=NULL)+
      apatheme

# only temperature
mydf <- ggpredict(glm_presence, c("meanTemp[n=200]")) 

        plot(mydf,show_data = TRUE, show_ci = TRUE, jitter=0.05)+
          labs(x = 'Temperature (C)', y = 'Probability of myxozoan presence (%)',title=NULL)+
          apatheme

# only flow
mydf <- ggpredict(glm_presence, c("meanFlow[n=200]")) 

          plot(mydf,rawdata=TRUE,jitter=0.05)+
        labs(x = 'Stream flow (m3/sec)', y = 'Probability of myxozoan presence (%)',title=NULL)+
        apatheme
          
          
## GLM for the probability of finding myxozoans for HYBNUC
glm_presence <- glmer(psite_presence ~ scale(YearCollected)*
                                  scale(meanTemp)*CI*scale(meanFlow)*Parasite_genus+
                                  scaled_TL_mm+
                                  (1|IndividualFishID)+
                                  (1|MonthCollected),
                                data = hybnuc_myxo,family=binomial()) #did not converge
          
          
glm_presence <- glmmTMB(psite_presence ~ scale(YearCollected)*
                          scale(meanTemp)*CI*scale(meanFlow)*Parasite_genus+
                          scaled_TL_mm+
                          (1|IndividualFishID)+
                          (1|MonthCollected),
                        data = hybnuc_myxo,family=binomial()) #did not converge
          
          
glm_presence <- glmer(psite_presence ~ scale(YearCollected)+
                        scale(meanTemp)*CI*scale(meanFlow)+
                        CI*before_after+
                        scaled_TL_mm+
                        scale(meanTemp)*Parasite_genus+
                        scale(meanFlow)*Parasite_genus+
                        CI*Parasite_genus+
                        (1|IndividualFishID)+
                        (1|MonthCollected),
                      data = hybnuc_myxo,family=binomial()) #singular fit, diagnostics are OK


summary(glm_presence)
          
          
#Evaluate residuals
  #Not suitable for GLM-quasipoisson but for the other models
  s=simulateResiduals(fittedModel=glm_presence,n=250)
  s$scaledResiduals
  plot(s)
          
          
# With the plot()function
plot_model(glm_presence,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
  
# only temperature
mydf <- ggpredict(glm_presence, c("CI","before_after")) 

plot(mydf,show_data = TRUE, show_ci = TRUE, jitter=0.05)+
  labs(x = 'Temperature (C)', y = 'Probability of myxozoan presence (%)',title=NULL)+
  apatheme

# only flow
mydf <- ggpredict(glm_presence, c("meanFlow[n=200]")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Stream flow (m3/sec)', y = 'Probability of myxozoan presence (%)',title=NULL)+
  apatheme

### Myxozoan abundace across time----

#Set a theme for all your plots
apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


# Plot myxozoan abundance against time 

pimvig_abundance <-ggerrorplot(pimvig_count, x = "YearCollected", y = "psite_count",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Pimephales vigilax")+
  xlab("Year")+ylab("Abundance")

pimvig_abundance

# Plot myxozoan abundance against time - GAMAFF

gamaff_abundance <-ggerrorplot(gamaff_count, x = "YearCollected", y = "psite_count",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Gambusia affinis")+
  xlab("Year")+ylab("Abundance")

gamaff_abundance

# Plot myxozoan abundance against time - ICTPUN

ictpun_abundance <-ggerrorplot(ictpun_count, x = "YearCollected", y = "psite_count",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Ictalurus punctatus")+
  xlab("Year")+ylab("Abundance")

ictpun_abundance


# NOTATH does not play a role here because it is only infected with myxidium and that can only be analyzed in terms of presence and absence

# Plot myxozoan prevalence against time - CARVEL

carvel_abundance <-ggerrorplot(carvel_count, x = "YearCollected", y = "psite_count",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Carpiodes velifer")+
  xlab("Year")+ylab("Abundance")

carvel_abundance


# Plot myxozoan prevalence against time - HYBNUC

hybnuc_abundance <-ggerrorplot(hybnuc_count, x = "YearCollected", y = "psite_count",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Hybognathus nuchalis")+
  xlab("Year")+ylab("Abundance")


hybnuc_abundance


# GLM for the abundance of  myxozoans
glm_abundance <- glmmTMB(psite_count ~ 
                           meanTemp*CI*meanFlow+
                           CI*before_after+
                           meanTemp*Parasite_genus+
                           meanFlow*Parasite_genus+
                           CI*Parasite_genus+
                           StandardLength_mm+
                           (1|MonthCollected),
                         data = myxo_count,
                         family = nbinom1) #not converging

glm_abundance <- glmer.nb(psite_count ~ YearCollected+
                        meanTemp*CI*meanFlow+
                        CI*before_after+
                        StandardLength_mm+
                        (1|MonthCollected),
                      data = carvel_count) #terrible fit

glm_abundance <- glmmTMB(psite_count ~ 
                            meanTemp*CI*meanFlow+
                            CI*before_after+
                            StandardLength_mm+
                            (1|MonthCollected),
                          data = carvel_count,
                         family = nbinom2) #terrible fit


summary(glm_abundance)

#Evaluate residuals
#Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_abundance,n=250)
s$scaledResiduals
plot(s)


#With the plot()function
mydf <- ggpredict(glm_abundance, c("meanFlow [all]","CI","before_after")) 
mydf <- ggpredict(glm_abundance, c("meanTemp [all]","CI","before_after")) 
mydf <- ggpredict(glm_abundance, c("before_after","CI")) 

apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

plot(mydf,show_data=TRUE, jitter = 0.3,dodge = TRUE)+
  labs(x = 'Treatment', y = 'Abundance of myxozoans',title=NULL)+
  apatheme

### Streamflow visualization----

# Plot stream flow as means throught time

streamflow_plot <-ggerrorplot(streamflow_ms, x = "Year", y = "Result",
                                ggtheme = theme_bw(),rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  ggtitle("Streamflow")+
  xlab("Year")+ylab("Streamflow (m3/sec)")

streamflow_plot

# Plot stream flow as raw data per year

ggplot(streamflow_ms, aes(x= as.factor(YearCollected),
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

# Plot stream flow as raw data per month

ggplot(Full_dataset_physical, aes(x= meanFlow,
                          y=psite_count))+
  facet_wrap("before_after")+
  geom_point()+apatheme+
  ggtitle("Streamflow per month")+
  xlab("Month")+ylab("Streamflow (m3/sec)")


### Water Temperature USGS visualization----
# Plot stream flow as means throughout time 

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
  geom_hline(yintercept=15, linetype="dashed", color = "black", size=0.5)+
  geom_hline(yintercept=25, linetype="dashed", color = "black", size=0.5)+
  geom_hline(yintercept=35, linetype="dashed", color = "black", size=0.5)

hist(wtemp$Result)
# USGS data merged with our data, to double check
# Water temperature per month

ggplot(Full_dataset_physical, aes(x= as.factor(season_temp),
                  y=meanTemp))+
  geom_point()+apatheme+
  ggtitle("Water temperature per month")+
  xlab("Month")+ylab("Temperature (degC)")+
  geom_hline(yintercept=20, linetype="dashed", color = "black", size=0.5)


# Plot Temperature against flow

ggplot(Full_dataset_physical, aes(x= meanTemp,
                                  y=meanFlow))+
  geom_point()+apatheme+
  ggtitle("Streamflow as as function of temperature")+
  xlab("Temperature")+ylab("Streamflow")


ggplot(carvel_myxo, aes(x= YearCollected,
                                  y=meanTemp,color=season_temp))+
  geom_point()+apatheme+
  facet_wrap("Parasite_genus")+
  ggtitle("Streamflow as as function of temperature")+
  xlab("Temperature")+ylab("Streamflow")

# Plot Temperature against year

ggplot(Full_dataset_physical, aes(x= YearCollected,
                                  y=meanTemp))+
  geom_point()+apatheme+
  ggtitle("Streamflow as as function of temperature")+
  xlab("Temperature")+ylab("Streamflow")
