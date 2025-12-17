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
library(performance)
library(patchwork)
library(forcats)

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
myx_stwy <- subset(full_dataset_myxo, psite_spp.x == "MYX.STWY")
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

#### FIGURE 1----
## Making map representing the full myxozoan sampling points (based on "full_dataset_myxo")

# Just make sure everything looks okay and is making sense
devtools::install_github("stadiamaps/ggmap") 
library("ggmap")
register_stadiamaps("45588e55-a703-4f6e-9ef0-ee6ed672b73a", write = TRUE)
library(ggspatial)

# Add the pulp mill outflow and USGS gauge (30.765582, -89.832)
extra_points <- data.frame(
  Longitude = c(-89.8209072,-89.832),  # Replace with your desired longitudes
  Latitude = c(30.79324276, 30.765582),    # Replace with your desired latitudes
  CI = c("Gauge", "Outflow")            # Optional: for mapping or color/shape grouping
)

all_points <- rbind(
  full_dataset_myxo[, c("Longitude", "Latitude", "CI")],
  extra_points[, c("Longitude", "Latitude", "CI")]
)

bounds <- c(
  left = min(all_points$Longitude) - 0.01,
  bottom = min(all_points$Latitude) - 0.01,
  right = max(all_points$Longitude) + 0.01,
  top = max(all_points$Latitude) + 0.01
)

## Create map
map_check <- get_stadiamap(bounds, zoom = 13, maptype = "stamen_toner_lite") %>%
  ggmap() +
  geom_point(
    data = all_points,
    aes(x = Longitude, y = Latitude, shape = CI, fill = CI),
    size = 4, color = "black"
  ) +
  annotation_north_arrow(
    location = "br",     # "tl", "tr", "bl", "br" for corner position
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1, "cm"),
    width = unit(1, "cm")
  ) +
  scale_shape_manual(  name = "Legend:",values = c("control" = 21, "impact" = 21, "Gauge" = 22, "Outflow" = 24),
                       labels = c("control" = "Control", "impact" = "Impact", "Gauge" = "Gauge","Outflow" = "Outflow")
  ) +
  scale_fill_manual(name = "Legend:",values = c("control" = "#4393c3", "impact" = "#b2182b", "Gauge" = "black", "Outflow" = "black"),
                    labels = c("control" = "Control", "impact" = "Impact", "Gauge" = "Gauge","Outflow" = "Outflow")
  ) +
  xlab("") + ylab("") +
  theme(
    legend.position = "right",
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 8),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1
  ))

map_check

ggsave(file="/Users/dakeishladiaz-morales/TUBRI/Manuscripts/Myxozoans/Figures/Figure1_map.png", width=240, height=180, dpi=1000, units = "mm")
ggsave(file="/Users/dakeishladiaz-morales/TUBRI/Manuscripts/Myxozoans/Figures/Figure1_map.pdf", width=240, height=180, dpi=1000, units = "mm")
ggsave(file="/Users/dakeishladiaz-morales/TUBRI/Manuscripts/Myxozoans/Figures/Figure1_map.tif", width=240, height=180, dpi=1000, units = "mm")

#### FIGURE 2-----

## Make plot showing the prevalence of infection of fish per parasite
# What is the number of fish we dissected per fish species
fish_dissected <- full_dataset_myxo %>% group_by(Fish_sp.x) %>% 
  summarise(number_dissected=length(unique(IndividualFishID)),
            .groups = 'drop') 

# Overview of abundance
summary_abundance <- full_dataset_myxo %>%
  group_by(psite_spp.x, Fish_sp.x) %>%            
  summarize(min = min(psite_count),
            max = max(psite_count),
            mean = mean(psite_count),.groups = "drop") # This part ensures that we have one fish per parasite genus

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
                     labels = c("Chloromyxum" = expression(italic("Chloromyxum")), 
                                "Henneguya"= expression(italic("Henneguya")), 
                                "Myxidium" = expression(italic("Myxidium")), 
                                "Myxobolus" = expression(italic("Myxobolus")),
                                "Thelohanellus" = expression(italic("Thelohanellus")),
                                "Unicauda" = expression(italic("Unicauda")))) +
  scale_fill_manual(name = "Parasite genus:",
                    values = c("Chloromyxum" = "#543005", 
                               "Henneguya" = "#8c510a", 
                               "Myxidium" = "#dfc27d", 
                               "Myxobolus" = "#c7eae5",
                               "Thelohanellus" = "#35978f",
                               "Unicauda" = "#01665e"), 
                    labels = c("Chloromyxum" = expression(italic("Chloromyxum")), 
                               "Henneguya"= expression(italic("Henneguya")), 
                               "Myxidium" = expression(italic("Myxidium")), 
                               "Myxobolus" = expression(italic("Myxobolus")),
                               "Thelohanellus" = expression(italic("Thelohanellus")),
                               "Unicauda" = expression(italic("Unicauda")))) +
  
  xlab("Fish species")+ylab("Prevalence of infection (%)")+
  scale_x_discrete(limits = c("Carpiodes velifer", 
                              "Pimephales vigilax", 
                              "Notropis atherinoides",
                              "Gambusia affinis",
                              "Hybognathus nuchalis",
                              "Ictalurus punctatus"),
                   labels = c("Carpiodes velifer"=expression(italic("Carpiodes velifer")*"\n"*"n=180"),
                              "Pimephales vigilax"=expression(italic("Pimephales vigilax")*"\n"*
                                                               "n = 193"),
                              "Notropis atherinoides"=expression(italic("Notropis atherinoides")*"\n"*
                                                               "n = 206"),
                              "Gambusia affinis"=expression(italic("Gambusia affinis")*"\n"*
                                                               "n = 208"),
                              "Hybognathus nuchalis"=expression(italic("Hybognathus nuchalis")*"\n"*
                                                               "n = 221"),
                              "Ictalurus punctatus"=expression(italic("Ictalurus punctatus")*"\n"*
                                                               "n = 88"))) + ylim(0,100)



ggsave(file="/Users/dakeishladiaz-morales/TUBRI/Manuscripts/Myxozoans/Figures/Figure2.png", width=240, height=180, dpi=1000, units = "mm")
ggsave(file="/Users/dakeishladiaz-morales/TUBRI/Manuscripts/Myxozoans/Figures/Figure2.pdf", width=240, height=180, dpi=1000, units = "mm")
ggsave(file="/Users/dakeishladiaz-morales/TUBRI/Manuscripts/Myxozoans/Figures/Figure2.tif", width=240, height=180, dpi=1000, units = "mm")


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
ggsave(file="Manuscripts/Myxozoans/Figures/Time/Carvel_MyxG_Year.tif", width=150, height=150, dpi=1000, units = "mm")

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
  labs(x = "Parasite + Host + Organ", y = "Effect of time") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.0, color = "darkgrey") +
  apatheme2 + coord_flip() +
  scale_x_discrete(labels = c(
    "MYX.G" = expression(italic("Myxobolus")*"+"*italic("C. velifer")*"+gills"),
    "MYX.F" = expression(italic("Myxobolus")*"+"*italic("C. velifer")*"+fins"),
    "MYX.THEL" = expression(italic("Thelohanellus")*"+"*italic("P. vigilax")*"+skin"),
    "MYX.SBAD" = expression(italic("Myxobolus")*"+"*italic("P. vigilax")*"+gills"),
    "MYX.GO" = expression(italic("Myxobolus")*"+"*italic("P. vigilax")*"+conn. tissue"),
    "MYX.TAIL" = expression(italic("Henneguya")*"+"*italic("I. punctatus")*"+fins"),
    "MYX.SP" = expression(italic("Myxobolus")*"+" *italic("N. atherinoides")*"+gills")
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
ggsave(file="Manuscripts/Myxozoans/Figures/Time/estimates_all.pdf", width=150, height=150, dpi=1000, units = "mm")

#### FIGURE 3----
apatheme5= theme_bw(base_size = 14,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



# --- Myx.G plot---
p1 <- plot(mydf, show_data = TRUE, show_residuals = FALSE, jitter = 0.01, color = c("#5aae61")) +
  labs(
    x = "Year",
    y = expression(atop(
      italic("Myxobolus")*" + "*italic("C. velifer")*" + gills",
      "Parasite abundance"
    )),
    title = NULL
  ) +
  apatheme5

p2 <- ggplot(myxgo_predict,
             aes(x = x, y = predicted)) +
  
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 6,
    shape = 21,
    fill = c("#4393c3", "#b2182b")
  ) +
  
  xlab("Site") +
  ylab(expression(atop(
    italic("Myxobolus")*" + "*italic("P. vigilax")*" + connective tissue",
    "Parasite abundance"
  ))) +
  
  theme_minimal() + apatheme5

# --- Combine them with bold labels ---
plot_grid(gg, p1, gg2, p2, labels = c("A", "B","C","D"), label_fontface = "bold")

ggsave(file="Manuscripts/Myxozoans/Figures/Figure3.png", width=300, height=300, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Figure3.pdf", width=300, height=300, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Figure3.tif", width=300, height=300, dpi=1000, units = "mm")


#### The effect of time on prevalence data----

## Carpiodes velifer - MYX.CM

# GLMM - which family is the best?
m1 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_cm,
              family = binomial(link = "logit")) 

m2 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_cm,
              family = binomial(link = "probit")) 

AIC(m1,m2)
summary(m1)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
plot(s)

# This one coverged and it's good
myxo_cm_time <- glmmTMB(psite_presence ~ scale(YearCollected)+
                          offset(logTL_mm)+
                          (1|site)+
                          (1|season),
                        data = myx_cm,
                        family = binomial(link = "logit")) 

summary <- summary(myxo_cm_time)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_cm/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_cm/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

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

#Evaluate model
tab_model(myxo_cm_time)
plot_model(myxo_cm_time,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model
mydf <- ggpredict(myxo_cm_time, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Probability of myxozoan infection (%)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Myxozoans/Figures/Time/myx_cm/Carvel_MyxCM_Presence_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/myx_cm/Carvel_MyxCM_Presence_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

## Gambusia affinis - MYX.STWY

# GLMM - which family is the best?
m1 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_stwy,
              family = binomial(link = "logit")) #116.772

m2 <- glmmTMB(psite_presence ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_stwy,
              family = binomial(link = "probit")) #118.394

AIC(m1,m2)
summary(m1)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

# The best is m1

myxo_stwy_time <- glmmTMB(psite_presence ~ scale(YearCollected)+
                            offset(logTL_mm)+
                            (1|site)+
                            (1|season),
                          data = myx_stwy,
                          family = binomial(link = "logit")) 

summary <- summary(myxo_stwy_time)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_stwy/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_stwy/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

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

#Evaluate model
tab_model(myxo_stwy_time)
plot_model(myxo_stwy_time,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(myxo_stwy_time, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Probability of myxozoan infection (%)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Myxozoans/Figures/Time/myx_stwy/GAMAFF_MyxSTWY_Presence_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/myx_stwy/GAMAFF_MyxSTWY_Presence_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

#### CI analysis ----
## Carpiodes velifer - MYX.CM

# GLMM - which family is the best?
m1 <- glmmTMB(psite_presence ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_cm,
              family = binomial(link = "logit")) 

m2 <- glmmTMB(psite_presence ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_cm,
              family = binomial(link = "probit")) 

AIC(m1,m2)
summary(m2)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
plot(s)

# This one coverged and it's good
myxo_cm_ci <- glmmTMB(psite_presence ~ CI+
                          offset(logTL_mm)+
                          (1|site)+
                          (1|season),
                        data = myx_cm,
                        family = binomial(link = "probit")) 

summary <- summary(myxo_cm_ci)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_cm/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_cm/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

## Extract estimates and confidence intervals for myx_cm

# Calculate confidence intervals and estimates
conf_int_myxscm <- confint(myxo_cm_ci)

# Combine estimates and confidence intervals into a data frame
results_myxscm <- data.frame(
  Estimate = conf_int_myxscm[2,3],
  Upper_CI = conf_int_myxscm[2,2],
  Lower_CI = conf_int_myxscm[2,1]
)

results_myxscm$Fish_sp <- "Carpiodes velifer"
results_myxscm$Parasite <- "MYX.CM"

#Evaluate model
tab_model(myxo_cm_ci)
plot_model(myxo_cm_ci,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model
mydf <- ggpredict(myxo_cm_ci, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Probability of myxozoan infection (%)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Myxozoans/Figures/Time/myx_cm/Carvel_MyxCM_Presence_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/myx_cm/Carvel_MyxCM_Presence_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

## Gambusia affinis - MYX.STWY

# GLMM - which family is the best?
m1 <- glmmTMB(psite_presence ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_stwy,
              family = binomial(link = "logit")) #116.772

m2 <- glmmTMB(psite_presence ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_stwy,
              family = binomial(link = "probit")) #118.394

AIC(m1,m2)
summary(m1)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

# The best is m1

myxo_stwy_CI <- glmmTMB(psite_presence ~ CI+
                            offset(logTL_mm)+
                            (1|site)+
                            (1|season),
                          data = myx_stwy,
                          family = binomial(link = "logit")) 

summary <- summary(myxo_stwy_CI)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_stwy/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_stwy/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

## Extract estimates and confidence intervals for myx_cm

# Calculate confidence intervals and estimates
conf_int_myxstwy <- confint(myxo_stwy_CI)

# Combine estimates and confidence intervals into a data frame
results_myxstwy <- data.frame(
  Estimate = conf_int_myxstwy[2,3],
  Upper_CI = conf_int_myxstwy[2,2],
  Lower_CI = conf_int_myxstwy[2,1]
)

results_myxstwy$Fish_sp <- "Gambusia affinis"
results_myxstwy$Parasite <- "MYX.STWY"

#Evaluate model
tab_model(myxo_stwy_CI)
plot_model(myxo_stwy_CI,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(myxo_stwy_time, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
  labs(x = 'Year', y = 'Probability of myxozoan infection (%)',title=NULL)+
  apatheme

ggsave(file="Manuscripts/Myxozoans/Figures/Time/myx_stwy/GAMAFF_MyxSTWY_Presence_Year.png", width=150, height=120, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/myx_stwy/GAMAFF_MyxSTWY_Presence_Year.pdf", width=150, height=120, dpi=1000, units = "mm")

## Carpiodes velifer - MYX.G
# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_g,
              family = nbinom1(link="log")) # AIC = 1530.354

m2 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_g,
              family = nbinom2(link="log")) # AIC = 1535.634

m3 <- glmmTMB(psite_count ~CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_g,
              family = nbinom1(link="sqrt")) # AIC = 1634.865

m4 <- glmmTMB(psite_count ~ CI+
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

myx_g_ci <- glmmTMB(psite_count ~CI+
                        offset(sqrt(TotalLength_mm))+
                        (1|site)+
                        (1|season),
                      data = myx_g,
                      family = nbinom1(link="sqrt"))

summary <- summary(myx_g_ci)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_g/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_g/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

## Extract estimates and confidence intervals for myx_g

# Calculate confidence intervals and estimates
conf_int_myxg <- confint(myx_g_ci)

# Combine estimates and confidence intervals into a data frame
results_myxg <- data.frame(
  Estimate = conf_int_myxg[2,3],
  Upper_CI = conf_int_myxg[2,2],
  Lower_CI = conf_int_myxg[2,1]
)

results_myxg$Fish_sp <- "Carpiodes velifer"
results_myxg$Parasite <- "MYX.G"

#Evaluate model
tab_model(myx_g_ci)
plot_model(myx_g_ci,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## Carpiodes velifer - MYX.F

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_f,
              family = nbinom1(link="log")) # 266.5972

m2 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_f,
              family = nbinom2(link="log")) # 266.6983

m3 <- glmmTMB(psite_count ~ CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_f,
              family = nbinom1(link="sqrt")) # AIC = 279.3361

m4 <- glmmTMB(psite_count ~ CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_f,
              family = nbinom2(link="sqrt")) # did not converge

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=myx_f_ci,n=250)
s$scaledResiduals
plot(s)

# m1 has the best AIC and the best residuals
myx_f_ci <- glmmTMB(psite_count ~ CI+
                        offset(logTL_mm)+
                        (1|site)+
                        (1|season),
                      data = myx_f,
                      family = nbinom1(link="log"))

summary <- summary(myx_f_ci)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_f/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_f/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

## Extract estimates and confidence intervals for myx_f

# Calculate confidence intervals and estimates
conf_int_myxf <- confint(myx_f_ci)

# Combine estimates and confidence intervals into a data frame
results_myxf <- data.frame(
  Estimate = conf_int_myxf[2,3],
  Upper_CI = conf_int_myxf[2,2],
  Lower_CI = conf_int_myxf[2,1]
)

results_myxf$Fish_sp <- "Carpiodes velifer"
results_myxf$Parasite <- "MYX.F"

#Evaluate model
tab_model(myx_f_ci)
plot_model(myx_f_ci,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

### Ictalurus punctatus
## Ictalurus puncatus - MYX.TAIL

# GLMM - which family is the best?

m1 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_tail,
              family = nbinom1(link="log")) # 91.75001

m2 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_tail,
              family = nbinom2(link="log")) # AIC = 87.96106, but very bad fit!!

m3 <- glmmTMB(psite_count ~ CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_tail,
              family = nbinom1(link="sqrt")) # AIC = 107.96720

m4 <- glmmTMB(psite_count ~ CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_tail,
              family = nbinom2(link="sqrt")) # Did not converge

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=myx_tail_ci,n=250)
s$scaledResiduals
plot(s)

# m1 has the best AIC and diagnostics and fit

myx_tail_ci <- glmmTMB(psite_count ~ CI+
                           offset(logTL_mm)+
                           (1|site)+
                           (1|season),
                         data = myx_tail,
                         family = nbinom1(link="log"))

summary <- summary(myx_tail_ci)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_tail/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_tail/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

## Extract estimates and confidence intervals for myx_f

# Calculate confidence intervals and estimates
conf_int_myxtail <- confint(myx_tail_ci)

# Combine estimates and confidence intervals into a data frame
results_myxtail <- data.frame(
  Estimate = conf_int_myxtail[2,3],
  Upper_CI = conf_int_myxtail[2,2],
  Lower_CI = conf_int_myxtail[2,1]
)

results_myxtail$Fish_sp <- "Ictalurus puncatus"
results_myxtail$Parasite <- "MYX.TAIL"

#Evaluate model
tab_model(myx_tail_ci)
plot_model(myx_tail_ci,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

### Notropis atherinoides
## Notropis atherinoides - MYX.SP
# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_sp,
              family = nbinom1(link="log")) # AIC = 479.7415

m2 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_sp,
              family = nbinom2(link="log")) # AIC = 477.0371

m3 <- glmmTMB(psite_count ~CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_sp,
              family = nbinom1(link="sqrt")) # AIC = 461.1927

m4 <- glmmTMB(psite_count ~ CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_sp,
              family = nbinom2(link="sqrt")) # AIC = NA

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel= m2,n=250)
s$scaledResiduals
plot(s)

# m3 has the best AIC

myx_sp_ci <- glmmTMB(psite_count ~ CI+
                      offset(logTL_mm)+
                      (1|site)+
                      (1|season),
                    data = myx_sp,
                    family = nbinom2(link="log")) 

summary <- summary(myx_sp_ci)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_sp/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_sp/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

## Extract estimates and confidence intervals for myx_f

# Calculate confidence intervals and estimates
conf_int_myxsp <- confint(myx_sp_ci)

# Combine estimates and confidence intervals into a data frame
results_myxsp <- data.frame(
  Estimate = conf_int_myxsp[2,3],
  Upper_CI = conf_int_myxsp[2,2],
  Lower_CI = conf_int_myxsp[2,1]
)

results_myxsp$Fish_sp <- "Notropis atherinoides"
results_myxsp$Parasite <- "MYX.SP"

#Evaluate model
tab_model(myx_sp_ci)
plot_model(myx_sp_ci,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

### Pimephales vigilax
## Pimephales vigilax - MYX.GO

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_go,
              family = nbinom1(link="log")) # AIC = 136.2273

m2 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_go,
              family = nbinom2(link="log")) # 136.6936

m3 <- glmmTMB(psite_count ~ CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_go,
              family = nbinom1(link="sqrt")) # AIC = 157.5028

m4 <- glmmTMB(psite_count ~ CI+
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

myx_go_ci <- glmmTMB(psite_count ~ CI+
                         offset(logTL_mm)+
                         (1|site)+
                         (1|season),
                       data = myx_go,
                       family = nbinom2(link="log")) 

summary <- summary(myx_go_ci)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_go/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_go/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

## Extract estimates and confidence intervals for myx_f

# Calculate confidence intervals and estimates
conf_int_myxgo <- confint(myx_go_ci)

# Combine estimates and confidence intervals into a data frame
results_myxgo <- data.frame(
  Estimate = conf_int_myxgo[2,3],
  Upper_CI = conf_int_myxgo[2,2],
  Lower_CI = conf_int_myxgo[2,1]
)

results_myxgo$Fish_sp <- "Pimephales vigilax"
results_myxgo$Parasite <- "MYX.GO"

#Evaluate model
tab_model(myx_go_ci)
plot_model(myx_go_ci,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Plot
myxgo_predict <-ggpredict(myx_go_ci,c("CI"))

myxgo_plot_ci <- ggplot(myxgo_predict,
                     aes(x = x, y = predicted)) +
  
  geom_errorbar(
    aes(ymin = conf.low, ymax = conf.high),
    width = 0.15,
    linewidth = 0.8
  ) +
  
  geom_point(
    size = 6,
    shape = 21,
    fill = c("#4393c3", "#b2182b")
  ) +
  
  xlab("Site") +
  ylab("Parasite abundance (# pseudocysts/mm fish total length)") +
  
  theme_minimal() + apatheme4

myxgo_plot

ggsave(file="Manuscripts/Myxozoans/Figures/Time/PIMVIG_MyxGO_Abundance_Year.png", width=150, height=180, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/PIMVIG_MyxGO_Abundance_Year.pdf", width=150, height=180, dpi=1000, units = "mm")

### Pimephales vigilax
## Pimephales vigilax - MYX.THEL

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_thel,
              family = nbinom1(link="log")) # AIC = 277.1146

m2 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_thel,
              family = nbinom2(link="log")) # AIC = 272.2518

m3 <- glmmTMB(psite_count ~ CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_thel,
              family = nbinom1(link="sqrt")) # AIC = 293.5116

m4 <- glmmTMB(psite_count ~ CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_thel,
              family = nbinom2(link="sqrt")) # did not converge

AIC(m1,m2,m3,m4)

# Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

# m2 has the best diagnostics. 

myx_thel_ci <- glmmTMB(psite_count ~ CI+
                           offset(logTL_mm)+
                           (1|site)+
                           (1|season),
                         data = myx_thel,
                         family = nbinom2(link="log"))

summary <- summary(myx_thel_ci)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_thel/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_thel/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

## Extract estimates and confidence intervals for myx_f

# Calculate confidence intervals and estimates
conf_int_myxthel <- confint(myx_thel_ci)

# Combine estimates and confidence intervals into a data frame
results_myxthel <- data.frame(
  Estimate = conf_int_myxthel[2,3],
  Upper_CI = conf_int_myxthel[2,2],
  Lower_CI = conf_int_myxthel[2,1]
)

results_myxthel$Fish_sp <- "Pimephales vigilax"
results_myxthel$Parasite <- "MYX.THEL"

#Evaluate model
tab_model(myx_thel_ci)
plot_model(myx_thel_ci,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

### Pimephales vigilax
## Pimephales vigilax - MYXO.SBAD

# GLMM - which family is the best?
m1 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myxo_sbad,
              family = nbinom1(link="log")) # AIC = 189.0973

m2 <- glmmTMB(psite_count ~ CI+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myxo_sbad,
              family = nbinom2(link="log")) # AIC = 185.9324

m3 <- glmmTMB(psite_count ~ CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myxo_sbad,
              family = nbinom1(link="sqrt")) # AIC = 196.6316

m4 <- glmmTMB(psite_count ~ CI+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myxo_sbad,
              family = nbinom2(link="sqrt")) # did not converge

AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

# m2 has the best AIC and the best diagnostics

myxo_sbad_ci <- glmmTMB(psite_count ~ CI+
                            offset(logTL_mm)+
                            (1|site)+
                            (1|season),
                          data = myxo_sbad,
                          family = nbinom2(link="log"))

summary <- summary(myxo_sbad_ci)

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Figures/Time/myx_sbad/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Figures/Time/myx_sbad/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

## Extract estimates and confidence intervals for myx_sbad

# Calculate confidence intervals and estimates
conf_int_myxsbad <- confint(myxo_sbad_ci)

# Combine estimates and confidence intervals into a data frame
results_myxsbad <- data.frame(
  Estimate = conf_int_myxsbad[2,3],
  Upper_CI = conf_int_myxsbad[2,2],
  Lower_CI = conf_int_myxsbad[2,1]
)

results_myxsbad$Fish_sp <- "Pimephales vigilax"
results_myxsbad$Parasite <- "MYX.SBAD"

#Evaluate model
tab_model(myxo_sbad_ci)
plot_model(myxo_sbad_ci,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)


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
gg2 <- ggplot(results_combined, aes(x = Parasite, y = Estimate, ymin = Lower_CI, ymax = Upper_CI, color = Parasite)) +
  geom_pointrange(size=0.6) +
  labs(x = "Parasite + Host + Organ", y = "Effect of pulp mill") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.0, color = "darkgrey") +
  apatheme2 + coord_flip() +
  scale_x_discrete(labels = c(
    "MYX.G" = expression(italic("Myxobolus")*"+"*italic("C. velifer")*"+gills"),
    "MYX.F" = expression(italic("Myxobolus")*"+"*italic("C. velifer")*"+fins"),
    "MYX.THEL" = expression(italic("Thelohanellus")*"+"*italic("P. vigilax")*"+skin"),
    "MYX.SBAD" = expression(italic("Myxobolus")*"+"*italic("P. vigilax")*"+gills"),
    "MYX.GO" = expression(italic("Myxobolus")*"+"*italic("P. vigilax")*"+conn. tissue"),
    "MYX.TAIL" = expression(italic("Henneguya")*"+"*italic("I. punctatus")*"+fins"),
    "MYX.SP" = expression(italic("Myxobolus")*"+" *italic("N. atherinoides")*"+gills")
  ),
  limits = c("MYX.GO","MYX.SBAD","MYX.THEL","MYX.SP","MYX.TAIL", "MYX.G","MYX.F")) +
  scale_color_manual(values = c(
    "MYX.G" = "black",
    "MYX.F" = "black",
    "MYX.THEL" = "black",
    "MYX.SBAD" = "black",
    "MYX.GO" = "red",
    "MYX.TAIL" = "black",
    "MYX.SP" = "black"
  )) +
  guides(color = "none")  # removes the legend

ggsave(file="Manuscripts/Myxozoans/Figures/Time/estimates_all.png", width=150, height=150, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/estimates_all.pdf", width=150, height=150, dpi=1000, units = "mm")
ggsave(file="Manuscripts/Myxozoans/Figures/Time/estimates_all.pdf", width=150, height=150, dpi=1000, units = "mm")

#### Non-point source - Only for MYX.G----

# Some data only spans between 1973 and 1994. So we will perform the analyses 
# only for that time period and for the control because I cannot assume that the 
# downstram impact sites are going to be represented by this data.

m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_temperature)+
                Elements_PC1+
                Elements_PC2+
                scale(mean_nitrogen)+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_g_yc,
              family=nbinom1(link="sqrt")) ## YearCollected and Elements_PC2 are highly correlated based on the VIF values

vif_values <- check_collinearity(m1)


m1 <- glmmTMB(psite_count ~ 
                scale(mean_temperature)*scale(mean_nitrogen)+
                scale(mean_temperature)*Elements_PC1+
                scale(mean_temperature)*Elements_PC2+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_g_yc,
              family=nbinom1(link="log")) #518.0079

m2 <- glmmTMB(psite_count ~ 
                scale(mean_temperature)*scale(mean_nitrogen)+
                scale(mean_temperature)*Elements_PC1+
                scale(mean_temperature)*Elements_PC2+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_g_yc,
              family=nbinom2(link="sqrt")) #503.0135

m3 <- glmmTMB(psite_count ~ 
                scale(mean_temperature)*scale(mean_nitrogen)+
                scale(mean_temperature)*Elements_PC1+
                scale(mean_temperature)*Elements_PC2+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_g_yc,
              family=nbinom1(link="sqrt")) #509.2888

m4 <- glmmTMB(psite_count ~ 
                scale(mean_temperature)*scale(mean_nitrogen)+
                scale(mean_temperature)*Elements_PC1+
                scale(mean_temperature)*Elements_PC2+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_g_yc,
              family=nbinom2(link="log")) #did not converge

m5 <- glmmTMB(psite_count ~ poly(mean_temperature,2) +
                scale(mean_temperature)*scale(mean_nitrogen)+
                scale(mean_temperature)*Elements_PC1+
                scale(mean_temperature)*Elements_PC2+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = myx_g_yc,
              family=nbinom1(link="log")) #519.1605

m6 <- glmmTMB(psite_count ~ poly(mean_temperature,2) +
                scale(mean_temperature)*scale(mean_nitrogen)+
                scale(mean_temperature)*Elements_PC1+
                scale(mean_temperature)*Elements_PC2+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_g_yc,
              family=nbinom2(link="sqrt")) #504.8136, really bad fit

m7 <- glmmTMB(psite_count ~ poly(mean_temperature,2) +
                scale(mean_temperature)*scale(mean_nitrogen)+
                scale(mean_temperature)*Elements_PC1+
                scale(mean_temperature)*Elements_PC2+
                offset(sqrt(TotalLength_mm))+
                (1|site)+
                (1|season),
              data = myx_g_yc,
              family=nbinom1(link="sqrt")) #509.4216

vif_values <- check_collinearity(m2)

AIC(m1,m2,m3,m4,m5,m6,m7)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

summary(m2)

# The best model is 
ms_model <- glmmTMB(psite_count ~ 
                      scale(mean_temperature)*scale(mean_nitrogen)+
                      scale(mean_temperature)*Elements_PC1+
                      scale(mean_temperature)*Elements_PC2+
                      offset(sqrt(TotalLength_mm))+
                      (1|site)+
                      (1|season),
                    data = myx_g_yc,
                    family=nbinom2(link="sqrt")) #503.0135


tab_model(ms_model)
plot_model(ms_model,auto.label = FALSE,type = "est",ci.lvl = 0.95,colors = c("black"))+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

summary(ms_model)

summary <- summary(ms_model)

# Calculate confidence intervals and estimates
conf_int_ms_model <- confint(ms_model)

# Combine estimates and confidence intervals into a data frame
results <- data.frame(
  Estimate = conf_int_ms_model[2,3],
  Upper_CI = conf_int_ms_model[2,2],
  Lower_CI = conf_int_ms_model[2,1]
)

# Combine estimates and confidence intervals into a data frame
results2 <- data.frame(
  Estimate = conf_int_ms_model[5,3],
  Upper_CI = conf_int_ms_model[5,2],
  Lower_CI = conf_int_ms_model[5,1]
)

results_myxstwy$Fish_sp <- "Gambusia affinis"
results_myxstwy$Parasite <- "MYX.STWY"

# Save model output 
capture.output(summary, file = "Manuscripts/Myxozoans/Results/multiple_stressors/summary.txt")

#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Results/multiple_stressors/diagnostics.pdf",width=8,height=5,colormodel="rgb")
plot(s) 
dev.off()

# Save estimate plot
#Save diagnostics
dev.off()
pdf("Manuscripts/Myxozoans/Results/multiple_stressors/estimates.pdf",width=8,height=5)
plot_model(ms_model,title="",axis.title = ("Estimates"),type = "est",ci.lvl = 0.95,colors = c("black"))+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
dev.off()

ggsave(file="Manuscripts/Myxozoans/Results/multiple_stressors/estimates.png", width=125, height=150, dpi=1000, units = "mm")

#Evaluate residuals
s=simulateResiduals(fittedModel=ms_model,n=250)
s$scaledResiduals
plot(s)

# Plot predictions

mydf <- ggpredict(ms_model, terms= c("mean_temperature[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=TRUE,color=c("#5aae61","#762a83"),jitter=0.1)+
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
mydf2 <- ggpredict(ms_model, terms= c("Elements_PC2[n=50]")) 

plot(mydf2,show_data=FALSE,show_residuals=TRUE,color=c("#5aae61"))+
  labs(x = "Elements PC2", y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme2

# Plot predictions

mydf3 <- ggpredict(ms_model, terms= c("Elements_PC1[n=100]","mean_temperature[18.5,20.5]")) 

p2 <- plot(mydf3,show_data=TRUE,show_residuals=FALSE,color=c("#74add1","#d73027"),jitter=0.1,alpha=0.2)+
  labs(color="Temperature (°C)",x = "Elements PC1", y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme
#"#fee090",

# Plot predictions

mydf4 <- ggpredict(ms_model, terms= c("Elements_PC2[n=100]","mean_temperature[18.2,20.2]")) 

p3 <- plot(mydf4,show_data=FALSE,show_residuals=TRUE,color=c("#74add1","#d73027"),jitter = 0.2)+
  labs(color="Temperature (°C)",x = "Elements PC2", y = 'Parasite abundance (# pseudocysts/fish)',title=NULL)+
  apatheme

#### FIGURE 4----

# Calculate confidence intervals and estimates
conf_int_myxg_ms <- confint(ms_model)

# Combine estimates and confidence intervals into a data frame
results_myxg_ms <- data.frame(
  Estimate = conf_int_myxg_ms[2,3],
  Upper_CI = conf_int_myxg_ms[2,2],
  Lower_CI = conf_int_myxg_ms[2,1]
)

# Convert matrix to data frame and add term names
conf_int_myxg_ms_df <- as.data.frame(conf_int_myxg_ms)
conf_int_myxg_ms_df$term <- rownames(conf_int_myxg_ms)

# Rename columns for clarity
colnames(conf_int_myxg_ms_df) <- c("Lower_CI", "Upper_CI", "Estimate", "term")

# Clean and order the terms
conf_int_myxg_ms_clean <- conf_int_myxg_ms_df %>%
  filter(!grepl("Std.Dev.|Intercept", term)) %>%
  mutate(
    term_clean = term %>%
      gsub("scale\\(", "", .) %>%
      gsub("\\)", "", .) %>%
      gsub("_", " ", .),
    color_group = ifelse(term_clean %in% c(
      "mean temperature",
      "Elements PC2",
      "mean temperature:Elements PC1",
      "mean temperature:Elements PC2"
    ), "red", "black")
  )

# Set term_clean as a factor in reversed order
conf_int_myxg_ms_clean$term_clean <- factor(conf_int_myxg_ms_clean$term_clean, levels = rev(conf_int_myxg_ms_clean$term_clean))

# Plot
gg_fig4 <- ggplot(conf_int_myxg_ms_clean,
             aes(x = term_clean, y = Estimate,
                 ymin = Lower_CI, ymax = Upper_CI, color = color_group)) +
  geom_pointrange(size = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.0, color = "darkgrey") +
  scale_color_identity() +
  coord_flip() +
  labs(x = "", y = "z-tranformed estimates (± 95% CI)") +
  apatheme2 +
  theme(axis.text.y = element_text(size = 10, face = "bold"))

# Save
ggsave("Manuscripts/Myxozoans/Figures/Figure4.png", gg_fig4, width = 150, height = 150, units = "mm", dpi = 1000)
ggsave("Manuscripts/Myxozoans/Figures/Figure4.pdf", gg_fig4, width = 150, height = 150, units = "mm", dpi = 1000)
ggsave("Manuscripts/Myxozoans/Figures/Figure4.tif", gg_fig4, width = 150, height = 150, units = "mm", dpi = 1000)


#### Figure 5----
# Set a theme for all your plots
apatheme3= theme_bw(base_size = 10,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Prepare each plot 
pA <- plot(mydf, show_data = TRUE, show_residuals = TRUE, color = c("#5aae61", "#762a83"), jitter = 0.1) +
  labs(x = "Temperature (°C)", y = "Parasite abundance", title = NULL) +
  apatheme3 +
  ggtitle("A") +
  theme(plot.title = element_text(face = "bold", hjust = -0.2))

mydf2 <- ggpredict(ms_model, terms = c("Elements_PC2[n=50]"))

pB <- plot(mydf2, show_data = FALSE, show_residuals = TRUE, color = c("#5aae61")) +
  labs(x = "Elements PC2", y = "Parasite abundance", title = NULL) +
  apatheme3 +
  ggtitle("B") +
  theme(plot.title = element_text(face = "bold", hjust = -0.2))

mydf3 <- ggpredict(ms_model, terms = c("Elements_PC1[n=100]", "mean_temperature[18.5,20.5]"))

pC <- plot(mydf3, show_data = TRUE, show_residuals = FALSE, color = c("#74add1", "#d73027"), jitter = 0.1, alpha = 0.2) +
  labs(color = "Temp (°C)", x = "Elements PC1", y = "Parasite abundance", title = NULL) +
  apatheme3 +
  ggtitle("C") +
  theme(plot.title = element_text(face = "bold", hjust = -0.2))

mydf4 <- ggpredict(ms_model, terms = c("Elements_PC2[n=100]", "mean_temperature[18.2,20.2]"))

pD <- plot(mydf4, show_data = FALSE, show_residuals = TRUE, color = c("#74add1", "#d73027"), jitter = 0.2) +
  labs(color = "Temp (°C)", x = "Elements PC2", y = "Parasite abundance", title = NULL) +
  apatheme3 +
  ggtitle("D") +
  theme(plot.title = element_text(face = "bold", hjust = -0.2))

# Combine in 2x2 layout
final_plot <- (pA | pB) / (pC | pD)

# Save
ggsave("Manuscripts/Myxozoans/Figures/Figure5.png", final_plot, width = 200, height = 125, units = "mm", dpi = 300)
ggsave("Manuscripts/Myxozoans/Figures/Figure5.tif", final_plot, width = 200, height = 125, units = "mm", dpi = 300)

#### Supplementary material----

## Supplementary Figure 1 - Prevalence
# Set plot aesthetic
apatheme4= theme_bw(base_size = 14,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        strip.text.y = element_text(
          face = "italic",
          size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Calculate prevalence for each year, CI, parasite genus, and fish

summary_psitegenus_prevalence <- prevalence_myxo %>%
  group_by(Parasite_genus, YearCollected, Fish_sp.x, CI, IndividualFishID) %>%
  summarise(
    infected = as.integer(sum(psite_presence, na.rm = TRUE) > 0),
    .groups = "drop"
  )

nfish_table_b <- summary_psitegenus_prevalence %>%
  distinct(Fish_sp.x, CI, YearCollected, IndividualFishID) %>%
  group_by(Fish_sp.x, CI, YearCollected) %>%
  summarise(
    nfish = n_distinct(IndividualFishID),
    .groups = "drop"
  )

summary_prevalence_c <- summary_psitegenus_prevalence %>%
  group_by(Parasite_genus, Fish_sp.x, YearCollected, CI) %>%
  summarise(
    n_infected = sum(infected),
    .groups = "drop"
  ) %>%
  left_join(
    nfish_table_b,
    by = c("Fish_sp.x", "CI", "YearCollected")
  ) %>%
  mutate(
    prevalence = (n_infected / nfish) * 100
  )

# Abbreviate the name of the fish species so it fits in the plot
prev_year <- summary_prevalence_c %>%
  mutate(
    Fish_sp.x = sub("^([A-Za-z])[a-z]+\\s+", "\\1. ", Fish_sp.x)
  )

# Make the plot
plot_prevalence <- ggplot(prev_year,
       aes(x = YearCollected,
           y = prevalence,
           color = Parasite_genus,
           group = Parasite_genus)) +
  
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  
  facet_grid(
    rows = vars(Fish_sp.x),
    cols = vars(CI),
  ) + ylim(0,100)+apatheme4 + labs(x="Year",y="Prevalence (%)")+
  scale_color_manual(name = "Parasite genus:",
                     values = c("Chloromyxum" = "#543005", 
                                "Henneguya" = "#8c510a", 
                                "Myxidium" = "#c7eae5", 
                                "Myxobolus" = "#dfc27d",
                                "Thelohanellus" = "#35978f",
                                "Unicauda" = "#01665e"),
                     labels = c("Chloromyxum" = expression(italic("Chloromyxum")), 
                                "Henneguya"= expression(italic("Henneguya")), 
                                "Myxidium" = expression(italic("Myxidium")), 
                                "Myxobolus" = expression(italic("Myxobolus")),
                                "Thelohanellus" = expression(italic("Thelohanellus")),
                                "Unicauda" = expression(italic("Unicauda")))) +
  scale_fill_manual(name = "Parasite genus:",
                    values = c("Chloromyxum" = "#543005", 
                               "Henneguya" = "#8c510a", 
                               "Myxidium" = "#c7eae5", 
                               "Myxobolus" = "#dfc27d",
                               "Thelohanellus" = "#35978f",
                               "Unicauda" = "#01665e"), 
                    labels = c("Chloromyxum" = expression(italic("Chloromyxum")), 
                               "Henneguya"= expression(italic("Henneguya")), 
                               "Myxidium" = expression(italic("Myxidium")), 
                               "Myxobolus" = expression(italic("Myxobolus")),
                               "Thelohanellus" = expression(italic("Thelohanellus")),
                               "Unicauda" = expression(italic("Unicauda")))) 
  
ggsave("Manuscripts/Myxozoans/Figures/Supp_Figure1.png", plot_prevalence, width = 210, height = 297, units = "mm", dpi = 300)
ggsave("Manuscripts/Myxozoans/Figures/Supp_Figure1.tif", plot_prevalence, width = 210, height = 297, units = "mm", dpi = 300)

## Supplementary Figure 2 - Abundance
# Set plot aesthetic
apatheme4= theme_bw(base_size = 14,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        strip.text.y = element_text(
          face = "italic",
          size = 12),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# Calculate abundance and descriptive stats for each year, CI, parasite genus, and fish

  abundance <- subset(full_dataset_myxo, Parasite_genus != "Chloromyxum"&
                        Parasite_genus != "Myxidium") # Remove species that do not present pseudocysts
  
  # Calculate prevalence and abundance for parasite genus, CI, fish species,
  fish_genus_abundance_year <- abundance %>%
    group_by(Parasite_genus, Fish_sp.x,YearCollected, CI, IndividualFishID) %>%
    summarise(
      fish_parasite_count = sum(psite_count, na.rm = TRUE),
      .groups = "drop"
    )
  
  summary_abundance_c <- fish_genus_abundance_year %>%
    group_by(Parasite_genus, Fish_sp.x, YearCollected,CI) %>%
    summarise(
      min_abundance  = min(fish_parasite_count),
      max_abundance  = max(fish_parasite_count),
      mean_abundance = mean(fish_parasite_count),
      .groups = "drop"
    )
  
# Shorten fish species names so that they fit in the plot
  abundance_year <- summary_abundance_c %>%
    mutate(
      Fish_sp.x = sub("^([A-Za-z])[a-z]+\\s+", "\\1. ", Fish_sp.x)
    )
  
# Plot it  
 abundance_year <- ggplot(abundance_year,
         aes(x = YearCollected,
             y = mean_abundance,
             color = Parasite_genus,
             group = Parasite_genus)) +
    
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    
    facet_grid(
      rows = vars(Fish_sp.x),
      cols = vars(CI),
      scales = "free_y"
    ) + apatheme4 + labs(x="Year",y="Mean abundance (# pseudocysts/fish)")+
    scale_color_manual(name = "Parasite genus:",
                       values = c(
                                  "Henneguya" = "#8c510a", 
                                  "Myxobolus" = "#dfc27d",
                                  "Thelohanellus" = "#35978f",
                                  "Unicauda" = "#01665e"),
                       labels = c( 
                                  "Henneguya"= expression(italic("Henneguya")), 
                                  "Myxobolus" = expression(italic("Myxobolus")),
                                  "Thelohanellus" = expression(italic("Thelohanellus")),
                                  "Unicauda" = expression(italic("Unicauda")))) +
    scale_fill_manual(name = "Parasite genus:",
                      values = c( 
                                 "Henneguya" = "#8c510a", 
                                 "Myxobolus" = "#dfc27d",
                                 "Thelohanellus" = "#35978f",
                                 "Unicauda" = "#01665e"), 
                      labels = c( 
                                 "Henneguya"= expression(italic("Henneguya")), 
                                 "Myxobolus" = expression(italic("Myxobolus")),
                                 "Thelohanellus" = expression(italic("Thelohanellus")),
                                 "Unicauda" = expression(italic("Unicauda")))) 
  
  ggsave("Manuscripts/Myxozoans/Figures/Supp_Figure2.png", abundance_year, width = 210, height = 297, units = "mm", dpi = 300)
  ggsave("Manuscripts/Myxozoans/Figures/Supp_Figure2.tif", abundance_year, width = 210, height = 297, units = "mm", dpi = 300)
  
  
## Supp Figure 3
  
  m2 <- lm(meanTemp~YearCollected, data=wtemp_summ)
  
  summary(m2) 
  
  mydf <- ggpredict(m2, terms= c("YearCollected[n=100]")) 
  
  temp_plot <- plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.01,color=c("#5aae61"))+
    labs(x = 'Year', y = 'Mean temperature (ºC)',title=NULL)+
    apatheme
  
  ggsave("Manuscripts/Myxozoans/Figures/Supp_Figure3.png", temp_plot, width = 150, height = 100, units = "mm", dpi = 300)
  ggsave("Manuscripts/Myxozoans/Figures/Supp_Figure3.tif", temp_plot, width = 150, height = 100, units = "mm", dpi = 300)

## Supp Figure 4
  cols_to_keep <- c(38, 40:51)
  
  subset_df <- myx_g_yc %>%
    select(YearCollected, all_of(cols_to_keep))

  plot_df <- subset_df %>%
    pivot_longer(
      cols = -YearCollected,
      names_to = "variable",
      values_to = "value"
    )

  facet_labels <- c(
    As = "As (µg/L)",
    Ba = "Ba (µg/L)",
    Cd = "Cd (µg/L)",
    Cr = "Cr (µg/L)",
    Cu = "Cu (µg/L)",
    Fe = "Fe (µg/L)",
    Hg = "Hg (µg/L)",
    mean_nitrogen = "N (mg/L)",
    Mg = "Mg (µg/L)",
    Mn = "Mn (µg/L)",
    Ni = "Ni (µg/L)",
    Pb = "Pb (µg/L)",
    Zn = "Zn (µg/L)"
  )
  
  
 plot_supf4 <- ggplot(plot_df, aes(x = YearCollected, y = value)) +
    geom_point(
      alpha = 0.7,
      size = 2,
      position = position_jitter(width = 0.2, height = 0)
    ) +
    facet_wrap(~ variable,
               scales = "free_y",
               ncol = 3,   
               labeller = as_labeller(facet_labels)) +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.minor = element_blank()
    ) +
    labs(
      x = "Year",
      y = "Concentration"
    ) + apatheme4

  ggsave("Manuscripts/Myxozoans/Figures/Supp_Figure4.png", plot_supf4, width = 210, height = 297, units = "mm", dpi = 300)
  ggsave("Manuscripts/Myxozoans/Figures/Supp_Figure4.tif", plot_supf4, width = 210, height = 297, units = "mm", dpi = 300)
  

## Table 2
  
# Calculate prevalence and abundance for parasite genus, CI, fish species,
  fish_genus_abundance <- full_dataset_myxo %>%
    group_by(Parasite_genus, Fish_sp.x, CI, IndividualFishID) %>%
    summarise(
      fish_parasite_count = sum(psite_count, na.rm = TRUE),
      .groups = "drop"
    )
  
  summary_abundance_b <- fish_genus_abundance %>%
    group_by(Parasite_genus, Fish_sp.x, CI) %>%
    summarise(
      min_abundance  = min(fish_parasite_count),
      max_abundance  = max(fish_parasite_count),
      mean_abundance = mean(fish_parasite_count),
      .groups = "drop"
    )
  
  
  ## Prevalence
  summary_psitegenus <- prevalence_myxo %>%
    group_by(Parasite_genus, Fish_sp.x, CI, IndividualFishID) %>%
    summarise(
      infected = as.integer(sum(psite_presence, na.rm = TRUE) > 0),
      .groups = "drop"
    )
  
  nfish_table <- summary_psitegenus %>%
    distinct(Fish_sp.x, CI, IndividualFishID) %>%
    group_by(Fish_sp.x, CI) %>%
    summarise(
      nfish = n_distinct(IndividualFishID),
      .groups = "drop"
    )
  
  summary_prevalence_b <- summary_psitegenus %>%
    group_by(Parasite_genus, Fish_sp.x, CI) %>%
    summarise(
      n_infected = sum(infected),
      .groups = "drop"
    ) %>%
    left_join(nfish_table, by = c("Fish_sp.x", "CI")) %>%
    mutate(
      prevalence = (n_infected / nfish)*100
    )
  
  # Merge tables
 merged_descriptive <- summary_abundance_b %>%
    left_join(summary_prevalence_b,
              by = c("Parasite_genus", "Fish_sp.x", "CI"))
 
 # Save to CSV
 write_csv(
   merged_descriptive,
   "/Users/dakeishladiaz-morales/TUBRI/Manuscripts/Myxozoans/Figures/prevalence_abundance_by_genus_species_site.csv"
 )
 
 
 ## Supplementary Table 1
 
 nfish_table_ST1 <- full_dataset_myxo %>%
   distinct(Fish_sp.x, CatalogNumber,collection_date,CI, Latitude,Longitude,
            IndividualFishID) %>%
   group_by(Fish_sp.x, CatalogNumber,collection_date,CI, Latitude,Longitude) %>%
   summarise(
     nfish = n_distinct(IndividualFishID),
     .groups = "drop"
   )
 
  View(nfish_table_ST1)
  
  # Save to CSV
  write_csv(
    nfish_table_ST1,
    "/Users/dakeishladiaz-morales/TUBRI/Manuscripts/Myxozoans/data/fish_metadata.csv"
  )
  
#### Explore other parasites----

## Import parasite dataset
full_dataset <- read_csv("data/processed/Full_dataset_with_psite_life_history_info_2024.11.21.csv")

## Take only myxozoan parasites
protozoa <- subset(full_dataset, Parasite_taxonomic_group == "Protozoa")

