## This script is meant to do the formal analysis on myxozoans
## The first set of analysis is looking at time. The second is looking at the impact of temperature, nitrogen, and element concentrations
## Created by Daki Diaz Morales (diazdakeishla@gmail.com)

## Required libraries----
library(readr)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(glmmTMB)
library(DHARMa)
library(sjPlot)
library(ggeffects)
library(ggplot2)
library(cowplot)

#### Data import and daughter dataframe creation----

# Set a theme for all your plots
apatheme= theme_bw(base_size = 14,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

## For the time analysis
full_dataset_myxo <- read_csv("Manuscripts/Myxozoans/data/myxozoans_time.csv")

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

## Summary stats
prevalence_myxo <- subset(full_dataset_myxo, !is.na(psite_presence))
abundance_myxo <- subset(full_dataset_myxo, !is.na(psite_count))

## Which parasites have a prevalence of more than 5%
summary_prevalence <- prevalence_myxo %>% group_by(Fish_sp.x,psite_spp.x) %>% 
  summarise(Prevalence=mean(psite_presence),
            .groups = 'drop') 

prevalence_higherthan5percent <- subset(summary_prevalence, Prevalence > 0.05)

# The parasites with more than 5% of prevalence are: MYX.CM, MYX.F, MYX.G, MYX.SWTY, MYX.TAIL, MYX.SP, MYX.GO, MYX.THEL, MYXO.SBAD

## Create subset per myxozoan species for the time analysis (myxozoans with a prevalence higher than 5%)
myx_cm <- subset(full_dataset_myxo, psite_spp.x == "MYX.CM")
myx_f <- subset(full_dataset_myxo, psite_spp.x == "MYX.F")
myx_g <- subset(full_dataset_myxo, psite_spp.x == "MYX.G")
myx_stwy <- subset(full_dataset_myxo, psite_spp.x == "MYX.SWTY")
myx_tail <- subset(full_dataset_myxo, psite_spp.x == "MYX.TAIL")
myx_sp <- subset(full_dataset_myxo, psite_spp.x == "MYX.SP")
myx_go <- subset(full_dataset_myxo, psite_spp.x == "MYX.GO")
myx_thel <- subset(full_dataset_myxo, psite_spp.x == "MYX.THEL")
myxo_sbad <- subset(full_dataset_myxo, psite_spp.x == "MYXO.SBAD")

## For the multiple stressor analysis
myxo_multiple_stressors <- read_csv("Manuscripts/Myxozoans/data/myxozoans_multiplestressors.csv")

## Create subset for myx.g for the non-point source pollution analysis
# We only do it for myx.g because it is the only that changed singificantly through time
# Restrict to the years that we are going to work with
myx_g_ms <- subset(myxo_multiple_stressors, psite_spp.x == "MYX.G")

myx_g_y <- subset(myx_g_ms, 
                  YearCollected > 1972 &
                    YearCollected < 1995 )

# Restrict to control because we have environmental data for that area of the river
myx_g_yc <- subset(myx_g_y, 
                   CI == "control" )

#### FIGURE 1-----

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


#### The effect of time on abundance data----

## Carpiodes velifer - MYX.G
# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_g,
              family = nbinom1(link="log")) # AIC = 1530.354

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_g,
              family = nbinom2(link="log")) # AIC = 1535.634

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_g,
              family = nbinom1(link="sqrt")) # AIC = 1634.865

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_g,
              family = nbinom2(link="sqrt")) # AIC = 1535.841

summary(m1)
AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=myx_g_time,n=250)
plot(s)

# m2 has the best diagnostics and an acceptable AIC. The AIC is 5 units larger than m1, but it provides better diagnostics.

myx_g_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_g,
              family = nbinom2(link="log")) # AIC = 1534.881

summary <- summary(myx_g_time)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_g/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_g/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

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

#Evaluate model
tab_model(myx_g_time)
plot_model(myx_g_time,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## Visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(myx_g_time, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/mm fish length)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Myxozoans/Figures/Time/Carvel_MyxG_Year.png", width=150, height=150, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/Carvel_MyxG_Year.pdf", width=150, height=150, dpi=1000, units = "mm")

## Carpiodes velifer - MYX.F

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_f,
              family = nbinom1(link="log")) # 266.5972

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_f,
              family = nbinom2(link="log")) # 266.6983

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_f,
              family = nbinom1(link="sqrt")) # AIC = 279.3361

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_f,
              family = nbinom2(link="sqrt")) # did not converge

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=myx_f_time,n=250)
s$scaledResiduals
plot(s)

# m1 has the best AIC and the best residuals
myx_f_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_f,
              family = nbinom1(link="log")) # 266.5972

summary <- summary(myx_f_time)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_f/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_f/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

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

#Evaluate model
tab_model(myx_f_time)
plot_model(myx_f_time,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(myx_f_time, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/mm fish length)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Myxozoans/Figures/Time/Carvel_MyxF_Abundance_Year.png", width=150, height=180, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/Carvel_MyxF_Abundance_Year.pdf", width=150, height=180, dpi=1000, units = "mm")

### Ictalurus punctatus
## Ictalurus puncatus - MYX.TAIL

# GLMM - which family is the best?

m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_tail,
              family = nbinom1(link="log")) # 91.75001

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_tail,
              family = nbinom2(link="log")) # AIC = 87.96106, but very bad fit!!

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_tail,
              family = nbinom1(link="sqrt")) # AIC = 107.96720

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_tail,
              family = nbinom2(link="sqrt")) # Did not converge

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=myx_tail_time,n=250)
s$scaledResiduals
plot(s)

# m1 has the best AIC and diagnostics and fit

myx_tail_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                           offset(logTL_mm)+
                           (1|site)+
                           (1|season),
                         data = myx_tail,
                         family = nbinom1(link="log")) # 91.75001

summary <- summary(myx_tail_time)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_tail/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_tail/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

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

#Evaluate model
tab_model(myx_tail_time)
plot_model(myx_tail_time,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(myx_tail_time, terms= c("YearCollected")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/mm fish length)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Myxozoans/Figures/Time/ICTPUN_MyxTail_Abundance_Year.png", width=150, height=180, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/ICTPUN_MyxTail_Abundance_Year.pdf", width=150, height=180, dpi=1000, units = "mm")

### Notropis atherinoides
## Notropis atherinoides - MYX.SP
# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_sp,
              family = nbinom1(link="log")) # AIC = 479.7415

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_sp,
              family = nbinom2(link="log")) # AIC = 477.0371

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_sp,
              family = nbinom1(link="sqrt")) # AIC = 461.1927

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_sp,
              family = nbinom2(link="sqrt")) # AIC = NA

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=myx_sp_time,n=250)
s$scaledResiduals
plot(s)

# m3 has the best AIC

myx_sp_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                         offset(sqrt(TotalLength_mm))+
                         (1|site)+
                         (1|season),
                       data = myx_sp,
                       family = nbinom1(link="sqrt")) # AIC = 461.1927

summary <- summary(myx_sp_time)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_sp/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_sp/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

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

#Evaluate model
tab_model(myx_sp_time)
plot_model(myx_sp_time,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(myx_sp_time, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Myxozoans/Figures/Time/NOTATH_MyxSP_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/NOTATH_MyxSP_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

### Pimephales vigilax
## Pimephales vigilax - MYX.GO

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_go,
              family = nbinom1(link="log")) # AIC = 136.2273

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_go,
              family = nbinom2(link="log")) # 136.6936

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_go,
              family = nbinom1(link="sqrt")) # AIC = 157.5028

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_go,
              family = nbinom2(link="sqrt")) # did not converge

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

#Check over dispersion
performance::check_overdispersion(m1)

# m2 has the best AIC and the best diagnostics, although there is a significant pattern in the residuals

myx_go_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                         offset(logTL_mm)+
                         (1|site)+
                         (1|season),
                       data = myx_go,
                       family = nbinom2(link="log")) # 136.6936

summary <- summary(myx_go_time)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_go/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_go/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

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
tab_model(myx_go_time)
plot_model(myx_go_time,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus
mydf <- ggpredict(myx_go_time, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/mm fish length)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Myxozoans/Figures/Time/PIMVIG_MyxGO_Abundance_Year.png", width=150, height=180, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/PIMVIG_MyxGO_Abundance_Year.pdf", width=150, height=180, dpi=1000, units = "mm")

### Pimephales vigilax
## Pimephales vigilax - MYX.THEL

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_thel,
              family = nbinom1(link="log")) # AIC = 277.1146

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_thel,
              family = nbinom2(link="log")) # AIC = 272.2518

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_thel,
              family = nbinom1(link="sqrt")) # AIC = 293.5116

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_thel,
              family = nbinom2(link="sqrt")) # did not converge

AIC(m1,m2,m3,m4)

# Evaluate residuals
s=simulateResiduals(fittedModel=myx_thel_time,n=250)
s$scaledResiduals
plot(s)

# m2 has the best diagnostics. 

myx_thel_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                           offset(logTL_mm)+
                           (1|site)+
                           (1|season),
                         data = myx_thel,
                         family = nbinom2(link="log")) # AIC = 272.2518

summary <- summary(myx_thel_time)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_thel/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_thel/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

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
tab_model(myx_thel_time)
plot_model(myx_thel_time,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(myx_thel_time, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/mm fish length)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Myxozoans/Figures/Time/PIMVIG_MyxTHEL_Abundance_Year.png", width=150, height=180, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/PIMVIG_MyxTHEL_Abundance_Year.pdf", width=150, height=180, dpi=1000, units = "mm")

### Pimephales vigilax
## Pimephales vigilax - MYXO.SBAD

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myxo_sbad,
              family = nbinom1(link="log")) # AIC = 189.0973

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myxo_sbad,
              family = nbinom2(link="log")) # AIC = 185.9324

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myxo_sbad,
              family = nbinom1(link="sqrt")) # AIC = 196.6316

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myxo_sbad,
              family = nbinom2(link="sqrt")) # did not converge

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=myxo_sbad_time,n=250)
s$scaledResiduals
plot(s)

# m2 has the best AIC and the best diagnostics

myxo_sbad_time <- glmmTMB(psite_count ~ scale(YearCollected)+
                            offset(logTL_mm)+
                            (1|site)+
                            (1|season),
                          data = myxo_sbad,
                          family = nbinom2(link="log")) # AIC = 185.9324

summary <- summary(myxo_sbad_time)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_sbad/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_sbad/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

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

#Evaluate model
tab_model(myxo_sbad_time)
plot_model(myxo_sbad_time,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(myxo_sbad_time, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Parasite abundance (# pseudocysts/mm fish length)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Myxozoans/Figures/Time/PIMVIG_MyxSBAD_Abundance_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/PIMVIG_MyxSBAD_Abundance_Year.pdf", width=150, height=120, dpi=1000, units = "mm")


## Combine all estimates and confidence intervals to make summary plot

## combine  data frames

results_combined <- bind_rows(results_myxg,results_myxf, results_myxtail,results_myxsp,results_myxthel,results_myxsbad,results_myxgo)

# What is the mean estimate across all fishes?
estimates_abundance_mean <- results_combined %>%  
  summarise(Est_mean=mean(Estimate),
            .groups = 'drop') 

# Data for horizontal lines: the mean estimate observed across all fish species. This helps you to see which fish is significantly deviating from the rest of the fish
# I took the values "-0.936" and "510" from the code lines 168-170.

# Set a theme for all your plots
apatheme2= theme_bw(base_size = 14,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Plot the estimates with confidence intervals
gg <- ggplot(results_combined, aes(x = Parasite, y = Estimate, ymin = Lower_CI, ymax = Upper_CI, color = Parasite)) +
  geom_pointrange(size=0.6) +
  labs(x = "Parasite - Host - Organ", y = "Effect of time") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.0, color = "darkgrey") +
  apatheme2 + coord_flip() +
  scale_x_discrete(labels = c(
    "MYX.G" = expression(italic("Myxobolus")*" spp. 2 "*italic("- C. velifer ")*"gills"),
    "MYX.F" = expression(italic("Myxobolus")*" spp. 1 "*italic("- C. velifer ")*"fins"),
    "MYX.THEL" = expression(italic("Thelohanellus")*" sp."*italic("- P. vigilax ")*"gills"),
    "MYX.SBAD" = expression(italic("Myxobolus")*" spp. 4"*italic("- P. vigilax ")*"fins"),
    "MYX.GO" = expression(italic("Myxobolus")*" spp. 5 "*italic("- P. vigilax ")*"gills"),
    "MYX.TAIL" = expression(italic("Henneguya")*" sp."*italic("- I. punctatus ")*"gills"),
    "MYX.SP" = expression(italic("Myxobolus")*" spp. 3" *italic("- N. atherinoides ")*"fins")
  ),
  limits = c("MYX.GO","MYX.SBAD","MYX.THEL","MYX.SP","MYX.TAIL", "MYX.G","MYX.F")) +
  scale_color_manual(values = c(
    "MYX.G" = "red",
    "MYX.F" = "black",
    "MYX.THEL" = "black",
    "MYX.SBAD" = "black",
    "MYX.GO" = "black",
    "MYX.TAIL" = "black",
    "MYX.SP" = "black"
  )) +
  guides(color = "none")  # removes the legend

ggsave(file="Manuscripts/Myxozoans/Figures/Time/estimates_all.png", width=150, height=150, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/estimates_all.pdf", width=150, height=150, dpi=1000, units = "mm")

### Figure 2----
# --- Myx.G plot---
p1 <- plot(mydf, show_data = TRUE, show_residuals = FALSE, jitter = 0.01, color = c("#5aae61")) +
  labs(
    x = "Year",
    y = expression(atop(
      italic("Myxobolus")*" spp. 2 "*italic("- C. velifer ")*"gills",
      "Parasite abundance (# pseudocysts/fish)"
    )),
    title = NULL
  ) +
  apatheme
# --- Combine them with bold labels ---
plot_grid(gg, p1, labels = c("A", "B"), label_fontface = "bold")

ggsave(file="Manuscripts/Myxozoans/Figures/Figure2.png", width=300, height=150, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Figure2.pdf", width=300, height=150, dpi=1000, units = "mm")

