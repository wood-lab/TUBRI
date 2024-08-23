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

full_dataset <- read_csv("data/processed/Full_dataset_with_psite_life_history_info_2024.08.17.csv")

## Make a subset for myxozoans

full_dataset_myxo <- subset(full_dataset, Parasite_taxonomic_group == "Myxozoa")

## Make Parasite_genus a factor

full_dataset_myxo$Parasite_genus <- as.factor(full_dataset_myxo$Parasite_genus)

## Create column for presence and absence of myxozoans
full_dataset_myxo$psite_presence <- ifelse(full_dataset_myxo$psite_count > 0,                # condition
                                               1,    # what if condition is TRUE
                                               0)       # what if condition is FALSE

## Remove columns with myxozoans that aren't confirmed to be myxozoans and therefore are not identified
full_dataset_myxo <- subset(full_dataset_myxo, !is.na(Parasite_genus))


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

##What is the impact of distance from pulp mill on the abundance of myxozoans?

#Set a theme for all your plots
apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

# Plot myxozoan abundance against distance from pulp mill (latitude)

ggplot(full_dataset_myxo_bca,aes(Latitude,psite_count, color = Parasite_genus,group = Parasite_genus))+
  geom_point(size=4)+
  facet_wrap("Fish_sp.x")+
  geom_vline(xintercept=30.76, linetype="dashed", color = "black", size=0.5)+
  apatheme


# Plot myxozoan abundance against time 

ggplot(full_dataset_myxo,aes(YearCollected,psite_count,shape = CI, color = Parasite_genus,group = Parasite_genus))+
  geom_point(size=2)+
  facet_wrap("Fish_sp.x")+
  geom_vline(xintercept=1972, linetype="dashed", color = "black", size=0.5)+
  apatheme


# GLM for the probability of finding myxozoans along the river
glm_presence<-glm(psite_presence ~ poly(Latitude,2),
             data = full_dataset_myxo_bca,family=binomial())

glm_presence<-glmmTMB(psite_presence ~ poly(Latitude,2)*Parasite_genus+
                        +(1|CatalogNumber),
                  data = full_dataset_myxo_bca,family=binomial())


summary(glm_presence)

#Evaluate residuals
#Not suitable for GLM-quasipoisson but for the other models
s=simulateResiduals(fittedModel=glm_presence,n=250)
s$scaledResiduals
plot(s)


#With the plot()function
mydf <- ggpredict(glm_presence, c("Latitude")) 

apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

plot(mydf,rawdata=TRUE)+
  labs(x = 'Latitude', y = 'Probability of myxozoan presence (%)',title=NULL)+
  apatheme+
  geom_vline(xintercept=30.76, linetype="dashed", color = "black", size=0.5)
  



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
                      apatheme+
                      xlab("Year")+ylab("Prevalence of infection")+
                      ylim(0,1)

pimvig_prevalence

# Plot myxozoan prevalence against time - GAMAFF

gamaff_prevalence <-ggerrorplot(gamaff_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

gamaff_prevalence

# Plot myxozoan prevalence against time - ICTPUN

ictpun_prevalence <-ggerrorplot(ictpun_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

ictpun_prevalence


# Plot myxozoan prevalence against time - NOTATH

notath_prevalence <-ggerrorplot(notath_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

notath_prevalence

# Plot myxozoan prevalence against time - CARVEL

carvel_prevalence <-ggerrorplot(carvel_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

carvel_prevalence


# Plot myxozoan prevalence against time - HYBNUC

hybnuc_prevalence <-ggerrorplot(hybnuc_myxo, x = "YearCollected", y = "psite_presence",
                                ggtheme = theme_bw(), color="CI",rawdata=TRUE,
                                position=position_dodge(0.5),width=0.00, size=0.3)+
  facet_wrap("Parasite_genus")+
  apatheme+
  xlab("Year")+ylab("Prevalence of infection")+
  ylim(0,1)

hybnuc_prevalence





