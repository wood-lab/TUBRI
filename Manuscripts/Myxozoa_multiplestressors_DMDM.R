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

# Adding stream flow to data
# Read physical data
physicalUSGS <- read_excel("data/geospatial/Physicaldata_USGS.xlsx")

# Make a subset only for stream flow in m3/sec
streamflow <- subset(physicalUSGS, Result_Characteristic=="Stream flow")
streamflow_ms <- subset(streamflow, Unit=="m3/sec")

# Summarize by taking the mean stream flow per 'year'

streamflow_mean <- streamflow_ms %>% group_by(YearCollected) %>% 
  summarise(meanFlow=mean(Result),
            .groups = 'drop')

# Merge with parasite data

Full_dataset_streamflow <- merge(full_dataset, streamflow_mean, by.x = "YearCollected", by.y = "YearCollected", all.x = TRUE)

# Adding water temperature data
# Make a subset of only water temperature
wtemp <- subset(physicalUSGS, Result_Characteristic=="Temperature, water")

# Summarize by taking the mean temperature 'yearseason'

wtemp_summ <- wtemp %>% group_by(YearCollected) %>% 
  summarise(meanTemp=mean(Result),medianTemp=median(Result),maxTemp=max(Result),
            .groups = 'drop')

# Merge with parasite data

Full_dataset_physical <- merge(Full_dataset_streamflow, wtemp_summ, by.x = "YearCollected", by.y = "YearCollected", all.x = TRUE)

## Create a season column

Full_dataset_physical$MonthCollected <- as.numeric(Full_dataset_physical$MonthCollected)

Full_dataset_physical <- Full_dataset_physical %>% 
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

Full_dataset_physical$season <- as.factor(Full_dataset_physical$season)


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
full_dataset_myxo$MonthCollected <- as.factor(full_dataset_myxo$MonthCollected)

## Revise reference levels

full_dataset_myxo$before_after <- relevel(full_dataset_myxo$before_after, ref = "before")


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


# Set a theme for all your plots
apatheme= theme_bw(base_size = 14,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

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

### Myxozoan abundace across time. CARVEL----

# Plot myxozoan prevalence against time - CARVEL

carvel_abundance <-ggerrorplot(carvel_count, x = "YearCollected", y = "psite_count",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Carpiodes velifer")+
  xlab("Year")+ylab("Abundance")+
  geom_vline(xintercept="1972", linetype="dashed", color = "black", size=0.5)


carvel_abundance

## GLM for CARVEL

# Both CatalogNumber and season as as random effects
glm_count <- glmmTMB(psite_count ~ scale(YearCollected)*
                       scale(medianTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+ 
                       (1|season),
                     data = carvel_count,family=nbinom1()) #horrible fit


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)*
                       scale(medianTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+ 
                       (1|season),
                     data = carvel_count,family=nbinom2()) #no convergence

# Both CatalogNumber and season as as random effects
glm_count1 <- glmmTMB(psite_count ~ 
                       scale(meanTemp)*scale(meanFlow)*
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+ 
                       (1|season),
                     data = carvel_count,family=nbinom1()) #horrible fit

glm_count2 <- glmmTMB(psite_count ~ 
                       scale(meanTemp)*scale(meanFlow)*
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+ 
                       (1|season),
                     data = carvel_count,family=nbinom2()) #horrible fit

glm_count3 <- glmmTMB(psite_count ~ 
                        scale(meanTemp)*scale(meanFlow)*
                        CI*before_after+
                        offset(scaled_TL_mm)+
                        (1|CatalogNumber)+ 
                        (1|season),
                      ziformula=~1,
                      data = carvel_count,family=nbinom1()) #horrible fit

glm_count4 <- glmmTMB(psite_count ~ 
                        scale(meanTemp)*scale(meanFlow)*
                        CI*before_after+
                        offset(scaled_TL_mm)+
                        (1|CatalogNumber)+ 
                        (1|season),
                      ziformula=~1,
                      data = carvel_count,family=nbinom2()) #horrible fit



AIC(glm_count1,glm_count2,glm_count3)

# only CatalogNumber as random effect
glm_count <- glmmTMB(psite_count ~ scale(YearCollected)*
                       scale(medianTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber),
                     data = carvel_count,family=nbinom2()) #no estimates possible, with neither nbinom1 or nbinom2

# only season as random effect
glm_count <- glmmTMB(psite_count ~ scale(YearCollected)*
                       scale(medianTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|season),
                     data = carvel_count,family=nbinom1()) #no estimates possible, with neither nbinom1 or nbinom2

summary(glm_count)

## For the random structure, I tried with season and catalognumber as REs and also without one or the other, with both nbinom1 and nibinom2 and the one with both REs performs better. 

# Both CatalogNumber and season as as random effects
glm_count <- glmmTMB(psite_count ~ scale(YearCollected)*scale(medianTemp)+
                       scale(YearCollected)*scale(meanFlow)+
                       scale(medianTemp)*scale(meanFlow)+
                       CI*scale(medianTemp)+
                       CI*scale(meanFlow)+
                       CI*before_after+
                       scaled_TL_mm+
                       (1|CatalogNumber)+ 
                       (1|season),
                     data = carvel_count,family=nbinom1()) # good convergence, residuals bad, 

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)*scale(medianTemp)+
                       scale(YearCollected)*scale(meanFlow)+
                       scale(medianTemp)*scale(meanFlow)+
                       CI*scale(medianTemp)+
                       CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+ 
                       (1|season),
                     data = carvel_count,family=nbinom2()) # overdispersed


glm_count <- glmer.nb(psite_count ~ scale(YearCollected)*scale(medianTemp)+
                       scale(YearCollected)*scale(meanFlow)+
                       scale(medianTemp)*scale(meanFlow)+
                       CI*scale(medianTemp)+
                       CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+ 
                       (1|season),
                     data = carvel_count) # did not converge



glm_count <- glmmTMB(psite_count ~ scale(YearCollected)*scale(medianTemp)+
                       scale(YearCollected)*scale(meanFlow)+
                       scale(medianTemp)*scale(meanFlow)+
                       CI*scale(medianTemp)+
                       CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+ 
                       (1|season),
                     ziformula=~1,
                     data = carvel_count,family=nbinom2()) # nbinom2 very bad fit, nbinom1 no convergence


summary(glm_count)

#Evaluate residuals
s=simulateResiduals(fittedModel=glm_count1,n=250)
s$scaledResiduals
plot(s)


#With the plot()function
plot_model(glm_count1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

plot_model(glm_count, type = "pred")

# Flow with CI and psite_genus

mydf <- ggpredict(glm_count, c("meanFlow[n=200]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  labs(x = 'Stream flow (m3/sec)', y = 'Abundance of myxozoans',title=NULL)+
  apatheme

# Temperature with CI and psite_genus

mydf <- ggpredict(glm_count, c("meanTemp[n=200]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  labs(x = 'Temperature (°C)', y = 'Abundance of myxozoans',title=NULL)+
  apatheme


# Temperature with CI and psite_genus

mydf <- ggpredict(glm_count, c("YearCollected[n=200]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Year', y = 'Abundance of myxozoans',title=NULL)+
  apatheme


# Flow with CI and psite_genus

mydf <- ggpredict(glm_count, c("before_after","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Clean water act', y = 'Abundance of myxozoans',title=NULL)+
  apatheme


### Myxozoan abundace across time. HYBNUC----

# Plot myxozoan prevalence against time - HYBNUC

hybnuc_abundance <-ggerrorplot(hybnuc_count, x = "YearCollected", y = "psite_count",
                               ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                               position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Hybognathus nuchalis")+
  xlab("Year")+ylab("Abundance")


hybnuc_abundance

## GLM for CARVEL

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)*Parasite_genus*
                       scale(meanTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       scaled_TL_mm+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = hybnuc_count,family=nbinom1()) #five way interaction is too much with both nbinom1, nbinom2


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)*Parasite_genus+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = hybnuc_count,family=nbinom1()) #bad convergence, four way interaction too much with both nbinom1, nbinom2

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       scale(meanTemp)*Parasite_genus+
                       scale(meanFlow)*Parasite_genus+
                       CI*Parasite_genus+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = hybnuc_count,family=nbinom2()) #bad convergence we simplify by creating a subset by parasite genus

# Create a subset for unicauda. For Myxobolus abundances are below three cysts per fish.

hybnuc_count_uni <- subset(hybnuc_count,Parasite_genus == "Unicauda")

# Full interaction model
glm_count <- glmmTMB(psite_count ~ scale(YearCollected)*
                       scale(meanTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       scaled_TL_mm+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = hybnuc_count_uni,family=nbinom1()) #four way interaction is too much with both nbinom1, nbinom2


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = hybnuc_count_uni,family=poisson()) #bad convergence

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = hybnuc_count_uni,family=nbinom1()) #bad convergence

glm_count <- glmmTMB(psite_count ~ 
                       scale(meanTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = hybnuc_count_uni,family=nbinom1()) #bad convergence

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = hybnuc_count_uni,family=nbinom1()) #bad convergence

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*scale(meanFlow)+
                       CI+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = hybnuc_count_uni,family=nbinom1()) #bad convergence




#Evaluate residuals
#Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_count,n=250)
s$scaledResiduals
plot(s)


#With the plot()function
plot_model(glm_count,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

# Flow with CI and psite_genus

mydf <- ggpredict(glm_count, c("meanFlow[n=200]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Stream flow (m3/sec)', y = 'Abundance of myxozoans',title=NULL)+
  apatheme

# Temperature with CI and psite_genus

mydf <- ggpredict(glm_count, c("meanTemp[n=200]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Temperature', y = 'Abundance of myxozoans',title=NULL)+
  apatheme


# Temperature with CI and psite_genus

mydf <- ggpredict(glm_count, c("YearCollected[n=200]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Temperature', y = 'Abundance of myxozoans',title=NULL)+
  apatheme


# Flow with CI and psite_genus

mydf <- ggpredict(glm_count, c("CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Clean water act', y = 'Abundance of myxozoans',title=NULL)+
  apatheme


### Myxozoan abundace across time. ICTPUN----


# Plot myxozoan abundance against time - ICTPUN

ictpun_abundance <-ggerrorplot(ictpun_count, x = "YearCollected", y = "psite_count",
                               ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                               position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Ictalurus punctatus")+
  xlab("Year")+ylab("Abundance")

ictpun_abundance



## GLM for CARVEL

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)*
                       scale(meanTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       scaled_TL_mm+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = ictpun_count,family=nbinom2()) #four way interaction is too much


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = ictpun_count,family=nbinom1()) #no convergence

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = ictpun_count,family=nbinom2()) #no convergence

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)+CI+scale(meanFlow)+
                       before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = ictpun_count,family=nbinom1()) #no convergence

summary(glm_count)

# nothing is converging, with neither nbinom1 and nbinom2, but makes sense because only 1963 there was an observation of 25 cysts per fish and the rest was near zero, so there is no change to analyze

### Myxozoan abundace across time. GAMAFF----

# Plot myxozoan abundance against time - GAMAFF

gamaff_abundance <-ggerrorplot(gamaff_count, x = "YearCollected", y = "psite_count",
                               ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                               position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Gambusia affinis")+
  xlab("Year")+ylab("Abundance")

gamaff_abundance


## GLM for CARVEL

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)*
                       scale(meanTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       scaled_TL_mm+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = gamaff_count,family=nbinom2()) #four way interaction is too much


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = gamaff_count,family=nbinom1()) # did not converge with either nbinom1 or nbinom2


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = gamaff_count,family=nbinom2()) # did not converge with either nbinom1 or nbinom2

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                        scale(meanTemp)+scale(meanFlow)+
                        CI*before_after+
                        offset(scaled_TL_mm)+
                        (1|CatalogNumber)+
                        (1|MonthCollected),
                      data = gamaff_count, family = nbinom2) # did not converge with either nbinom1 or nbinom2 


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)+scale(meanFlow)+
                       CI+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = gamaff_count, family = nbinom1) # did not converge with either nbinom1 or nbinom2 


glm_count <- glmmTMB(psite_count ~ 
                       scale(meanTemp)+scale(meanFlow)+
                       CI+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = gamaff_count, family = nbinom2) # converged and diagnostics are good

glm_count <- glmmTMB(psite_count ~ 
                       scale(meanTemp)+scale(meanFlow)+
                       CI+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = gamaff_count, family = nbinom1) # did not converge

summary(glm_count)


#Evaluate residuals
#Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_count,n=250)
s$scaledResiduals
plot(s)


#With the plot()function
plot_model(glm_count,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

# Flow with CI and psite_genus

mydf <- ggpredict(glm_count, c("meanFlow[n=200]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Stream flow (m3/sec)', y = 'Abundance of myxozoans',title=NULL)+
  apatheme

# Temperature with CI and psite_genus

mydf <- ggpredict(glm_count, c("meanTemp[n=200]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Temperature', y = 'Abundance of myxozoans',title=NULL)+
  apatheme


# Flow with CI 

mydf <- ggpredict(glm_count, c("CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Clean water act', y = 'Abundance of myxozoans',title=NULL)+
  apatheme



### Myxozoan abundace across time. PIMVIG----


# Plot myxozoan abundance against time 

pimvig_abundance <-ggerrorplot(pimvig_count, x = "YearCollected", y = "psite_count",
                               ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                               position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+ggtitle("Pimephales vigilax")+
  xlab("Year")+ylab("Abundance")

pimvig_abundance

## GLM for PIMVIG

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)*
                       scale(meanTemp)*CI*scale(meanFlow)*Parasite_genus+
                       CI*before_after+
                       scaled_TL_mm+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count,family=nbinom2()) #5 way interaction is too much


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       Parasite_genus*scale(meanTemp)+
                       Parasite_genus*scale(meanFlow)+
                       Parasite_genus*CI+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count,family=nbinom2()) #converged, diagnostics are good, but fit terrible

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       Parasite_genus*scale(meanTemp)+
                       Parasite_genus*scale(meanFlow)+
                       Parasite_genus*CI+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count,family=nbinom1()) #did not converge

glm_count <- glmer.nb(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       Parasite_genus*scale(meanTemp)+
                       Parasite_genus*scale(meanFlow)+
                       Parasite_genus*CI+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count) #did not converge


# To make the model a bit simpler, let's split the analysis in the two parasite genera

pimvig_count_myxobolus <- subset(pimvig_count, Parasite_genus == "Myxobolus")


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count_myxobolus,family=nbinom2()) #did not converge

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count_myxobolus,family=nbinom1()) #did not converge

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*scale(meanFlow)+
                       CI*scale(meanTemp)+
                       CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count_myxobolus,family=nbinom1()) # did not converge converged


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count_myxobolus,family=nbinom2()) # converged with nibnon2 (not with nbinom1) but bad diagnostics


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     zi=~1,
                     data = pimvig_count_myxobolus,family=nbinom2()) # converged, good diagnostics, horrible fit



glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)+scale(meanFlow)+
                       CI+before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     zi=~1,
                     data = pimvig_count_myxobolus,family=nbinom2()) # converged, good diagnostics, horrible fit


glm_count <- glmmTMB(psite_count ~ 
                       scale(meanTemp)+scale(meanFlow)+CI+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     zi = ~1,
                     data = pimvig_count_myxobolus,family=nbinom2()) # converged, good diagnostics, horrible fit


glm_count <- glmmTMB(psite_count ~ 
                       meanTemp+meanFlow+CI+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count_myxobolus,family=nbinom1()) # did not converge


summary(glm_count)


#Evaluate residuals
#Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_count,n=250)
s$scaledResiduals
plot(s)


#With the plot()function
plot_model(glm_count,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

# Flow with CI and psite_genus

mydf <- ggpredict(glm_count, c("meanFlow[n=200]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Stream flow (m3/sec)', y = 'Abundance of myxozoans',title=NULL)+
  apatheme

# Temperature with CI and psite_genus

mydf <- ggpredict(glm_count, c("meanTemp[all]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Temperature', y = 'Abundance of myxozoans',title=NULL)+
  apatheme


# Temperature with CI

mydf <- ggpredict(glm_count, c("YearCollected[n=200]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Temperature', y = 'Abundance of myxozoans',title=NULL)+
  apatheme


# now for thelohanellus

pimvig_count_thelohanellus <- subset(pimvig_count, Parasite_genus == "Thelohanellus")


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count_thelohanellus,family=nbinom1()) #did not converge

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)*CI*scale(meanFlow)+
                       scale(meanTemp)*scale(meanFlow)+
                       CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count_thelohanellus,family=nbinom2()) #did not converge

glm_count <- glmmTMB(psite_count ~ 
                       scale(meanTemp)*CI*scale(meanFlow)+
                       scale(meanTemp)*scale(meanFlow)+
                       CI*scale(meanFlow)+
                       CI*before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count_thelohanellus,family=nbinom1()) #did not converge

glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)+scale(meanFlow)+
                       CI+
                       before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count_thelohanellus,family=nbinom1()) #converged, bad diagnostics


glm_count <- glmmTMB(psite_count ~ scale(YearCollected)+
                       scale(meanTemp)+scale(meanFlow)+
                       CI+
                       before_after+
                       offset(scaled_TL_mm)+
                       (1|CatalogNumber)+
                       (1|MonthCollected),
                     data = pimvig_count_thelohanellus,family=nbinom2()) #converged, good diagnostics

summary(glm_count)


#Evaluate residuals
#Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_count,n=250)
s$scaledResiduals
plot(s)

#With the plot()function
plot_model(glm_count,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

# Flow with CI and psite_genus

mydf <- ggpredict(glm_count, c("meanFlow[n=200]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Stream flow (m3/sec)', y = 'Abundance of myxozoans',title=NULL)+
  apatheme

# Temperature with CI and psite_genus

mydf <- ggpredict(glm_count, c("meanTemp[all]","CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Temperature', y = 'Abundance of myxozoans',title=NULL)+
  apatheme

# Temperature with CI and psite_genus

mydf <- ggpredict(glm_count, c("CI")) 

plot(mydf,rawdata=TRUE,jitter=0.05)+
  labs(x = 'Temperature', y = 'Abundance of myxozoans',title=NULL)+
  apatheme


### Streamflow visualization----

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

ggplot(streamflow_mean, aes(x= as.factor(YearCollected),
                          y=meanFlow))+
  geom_bar(stat = "identity",color="#9ecae1",fill="#9ecae1")+apatheme+
  ggtitle("Streamflow per year")+
  xlab("Year")+ylab("Streamflow (m3/sec)")+
  geom_hline(yintercept=214, linetype="dashed", color = "black", size=0.5)

# what is the mean flow of the river throughout the years?
mean(streamflow_mean$meanFlow)

### Water Temperature USGS visualization----

# wtemp change throughout time

m1 <- lm(meanTemp~YearCollected, data=wtemp_summ)
m2 <- lm(medianTemp~YearCollected, data=wtemp_summ)

#Evaluate residuals
#Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

summary(m1)
summary(m2)

#With the plot()function
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
plot_model(m2,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

tab_model(m1)

# Flow with CI and psite_genus

mydf <- ggpredict(m1, c("YearCollected[n=100]")) 

plot(mydf,rawdata=TRUE,color="red")+
  labs(x = 'Year', y = 'Mean Annual Temperature (C)',title=NULL)+
  apatheme


# Plot water temperature as means throughout time 

wtemp_plot <-ggerrorplot(wtemp, x = "YearCollected", y = "Result",
                              ggtheme = theme_bw(),rawdata=TRUE,
                              position=position_dodge(0.5),width=0.00, size=0.3)+
  ggtitle("Water temperature")+
  xlab("Year")+ylab("Temperature (C)")

wtemp_plot


# Plot water temperature calculated mean against years 

meanwtemp_plot <-ggerrorplot(wtemp_summ, x = "YearCollected", y = "meanTemp",
                         ggtheme = theme_bw(),rawdata=TRUE,
                         position=position_dodge(0.5),width=0.00, size=0.3)+
  ggtitle("Water temperature")+
  xlab("Year")+ylab("mean Temperature (C)")

meanwtemp_plot

# Plot water temoerature median across years 

medianwtemp_plot <-ggerrorplot(wtemp_summ, x = "YearCollected", y = "medianTemp",
                             ggtheme = theme_bw(),rawdata=TRUE,
                             position=position_dodge(0.5),width=0.00, size=0.3)+
  ggtitle("Water temperature")+
  xlab("Year")+ylab("median Temperature (C)")

medianwtemp_plot

# Plot water temoerature max across years 

maxwtemp_plot <-ggerrorplot(wtemp_summ, x = "YearCollected", y = "maxTemp",
                               ggtheme = theme_bw(),rawdata=TRUE,
                               position=position_dodge(0.5),width=0.00, size=0.3)+
  ggtitle("Water temperature")+
  xlab("Year")+ylab("max Temperature (C)")

maxwtemp_plot


# Water temperature per month

ggplot(wtemp, aes(x= as.factor(MonthCollected),
                            y=Result))+
  geom_point()+apatheme+
  ggtitle("Water temperature per month")+
  xlab("Month")+ylab("Temperature (degC)")+
  geom_hline(yintercept=15, linetype="dashed", color = "black", size=0.5)+
  geom_hline(yintercept=25, linetype="dashed", color = "black", size=0.5)+
  geom_hline(yintercept=35, linetype="dashed", color = "black", size=0.5)
