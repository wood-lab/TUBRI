#This script is  ade analyze the monogenean data from the CAREER REU program
#Revised on July 17, 2024
#Written by Imani Jones and Daki Diaz


# If you need to change your working directory, paste the file structure between the quotation marks below.
# To get the file structure, open the folder and then "Get info" about that folder. Copy the file structure
# and paste it below.

setwd("C:/Users/imani/OneDrive - Tuskegee University/Desktop/TUBRI_Monogenea_Project/TUBRI")

full_dataset <- read_csv("data/processed/Full_dataset_physical_2024.09.25.csv")

full_dataset$logTL_mm <- log(full_dataset$TotalLength_mm)


### Code for monogenean MONO.IP in ICTPUNct abundance across time and pollution impact----
# Once your working directory is set, you're ready to read in the data!  If there are sub-folders inside your 
# working directory, you'll need to specify them as I have below.
# The <- command tells R what a thing is called.  So you can read the line below as,
# Look at this csv file, and name it pim_vig_data.  When I call pim_vig_data, I want you to give me the csv file.

#pim_vig_data<-read.csv("data/processed/Ictalurus_punctatus_processed_human_readable.csv")
#pim_vig_data<-read.csv("Pimephales_vigilax_processed_machine_readable")

library(readr)

#ICTPUNCT <- read_csv("data/processed/Ictalurus_punctatus_processed_machine_readable_UPDATED_2024.08.25.csv")

ICTPUNCT <- subset(full_dataset, Fish_sp.x=="Ictalurus punctatus")

# To see your data as a spreadsheet in a new tab, use the View command.

 #View(ICTPUNCT)

# You can also see your data in the console below by running the name of the dataset.

ICTPUNCT

#We made a subset based on our data to include only MONO.IP

ICTPUNCT_MONOIP <- subset(ICTPUNCT, psite_spp.x == 'MONO.IP')

#install.packages("ggplot2")
library(ggplot2)

#Set a theme for all your plots
apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

#gg plot is plotting function in r, it allows multiple variables within the plot.
#first part of code is showing what data set I want to be called from. "aes"-- command in r, x is the variable in my data set, y is the monogenes affected in this data set. 
#Color is showing my third variable in data aset. "+"- indicates adding aesthic to graph. 
ggplot(ICTPUNCT_MONOIP, aes(x= YearCollected,
                                y=psite_count,
                                ,color=CI,
                                group=CI))+
  geom_point()+ 
  labs(title = "Counts of Monogeneans in ICTPUNCT", x="Year", y="Monogenean abundance (# monogeneans/fish)",color= 
         "Pollution impact:")+apatheme+geom_vline(xintercept=1972, linetype="dashed", color = "black", size=0.5)

# geom point adds points to graph.
#labs=labels. title= title of graph, all words on graphs are quotes. X axis= label for x axis.
#y=label for y axis. color= represents the legend or key of labels.

#install.packages("stats")
#install.packages("lme4")
#install.packages("DHARMa")
#install.packages("glmmTMB")

library(stats)
library(lme4)
library(DHARMa)
library(glmmTMB)


#Is linear, no quadratic term. The diagnostics are terrible.
#Best one
glmmonoip1 <- glmmTMB(psite_count ~ YearCollected*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom1,data = ICTPUNCT_MONOIP)

#Models have good gianostics
glmmonoip1 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom1,data = ICTPUNCT_MONOIP)
glmmonoip2 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom2,data = ICTPUNCT_MONOIP)


#Models have good gianostics
glmmonoip1 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+offset(logTL_mm)+(1|CatalogNumber),family=nbinom1,data = ICTPUNCT_MONOIP)
glmmonoip2 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+offset(logTL_mm)+(1|CatalogNumber),family=nbinom2,data = ICTPUNCT_MONOIP)

# generalized linear model for count data with poisson or negative binomial distribution
# generalized linear mixed model for countr data with varying slopes

# YearCollected*Pollution = Year + Pollution + Year:Pollution

#Model selection. The lower the AIC, the better the model is.
AIC(glmmonoip1,glmmonoip2)

#After evaluating AIC, this one is the best.
glmmonoip1 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom1,data = ICTPUNCT_MONOIP)

summary(glmmonoip1)
tab_model(glmmonoip1)

#Evaluate diagnostics
library(DHARMa)
library(MASS)

s=simulateResiduals(fittedModel=glmmonoip1,n=250)
s$scaledResiduals
plot(s)

## final plot with predictions and CI
library(ggeffects)
plot_model(glmmonoip1)

mydf <- ggpredict(glmmonoip1, c("YearCollected [n=100]", "CI"), jitter=TRUE) 

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'))


plot(mydf, rawdata = TRUE, alpha = 0.3, dot.alpha = 0.6,colors=c("#969696","#ce1256"))+
  labs(x = 'Year', y = 'Monogenean abundance (# monogeneans/fish)',title=NULL)+
  apatheme

dev.off()


ggsave(file="REU/Imani J/Model plots/ICTPUNCT_Ligictaluridus pricei _modelplot.png", width=150, height=90, dpi=1000, units = "mm")
ggsave(file="REU/Imani J/Model plots/ICTPUNCT_Ligictaluridus pricei _modelplot.pdf", width=150, height=90, dpi=1000, units = "mm")



### Code for monogenean MONO.DACT in PIMVIG abundance across time and pollution impact----
# Once your working directory is set, you're ready to read in the data!  If there are sub-folders inside your 
# working directory, you'll need to specify them as I have below.
# The <- command tells R what a thing is called.  So you can read the line below as,
# Look at this csv file, and name it pim_vig_data.  When I call pim_vig_data, I want you to give me the csv file.

library(readr)

#PIMVIG <- read_csv("data/processed/Pimephales_vigilax_processed_machine_readable_UPDATED_2024.08.16.csv")

PIMVIG <- subset(full_dataset, Fish_sp.x=="Pimephales vigilax")

# To see your data as a spreadsheet in a new tab, use the View command.

# View(ICTPUNCT)

# You can also see your data in the console below by running the name of the dataset.

PIMVIG

#We made a subset based on our data to include only monogeneans

PIMVIG_MONO <- subset(PIMVIG, psite_spp.x == 'MONO.DACT'|psite_spp.x == 'MONO.GYRO'|psite_spp.x == 'MONO.LG')

View(PIMVIG_MONO)

#install.packages("ggplot2")
library(ggplot2)

#Set a theme for all your plots
apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

#gg plot is plotting function in r, it allows multiple variables within the plot.
#first part of code is showing what data set I want to be called from. "aes"-- command in r, x is the variable in my data set, y is the monogenes affected in this data set. 
#Color is showing my third variable in data aset. "+"- indicates adding aesthic to graph. 
ggplot(PIMVIG_MONO, aes(x= YearCollected,
                            y=psite_count,
                            ,color=CI,
                            group=CI))+
  geom_point()+ 
  labs(title = "Counts of monogenean species in Pimephales vigilax", x="Year", y="Monogenean abundance (# monogeneans/fish)",color= 
         "Pollution impact:")+apatheme+geom_vline(xintercept=1972, linetype="dashed", color = "black", size=0.5)+
  facet_wrap(~ psite_spp)


# Make a subset only for Pimvig mono.dact

PIMVIG_MONOdact <- subset(PIMVIG, psite_spp.x == 'MONO.DACT')

View(PIMVIG_MONOdact)

library(stats)
library(lme4)
library(DHARMa)
library(glmmTMB)

#Model passed diagnostics test
glmmonodact <- glmmTMB(psite_count ~ YearCollected*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom2,data = PIMVIG_MONOdact)
glmmonodact <- glmmTMB(psite_count ~ YearCollected*CI+offset(logTL_mm)+(1|CatalogNumber),family=nbinom2,data = PIMVIG_MONOdact)

# Did not pass KS Test. The Kolmogorov-Smirnov Test is a type of non-parametric test of the equality of discontinuous and continuous a 1D probability distribution that is used to compare the sample with the reference probability test (known as one-sample K-S Test) or among two samples (known as two-sample K-S test). A K-S Test quantifies the distance between the cumulative distribution function of the given reference distribution and the empirical distributions of given two samples, or between the empirical distribution of given two samples.
glmmonodact1 <- glmmTMB(psite_count ~ YearCollected*CI+TotalLength_mm+(1|CatalogNumber),family=poisson,data = PIMVIG_MONOdact)

#Best one
glmmonodact1 <- glmmTMB(psite_count ~ YearCollected*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom1,data = PIMVIG_MONOdact)
glmmonodact1 <- glmmTMB(psite_count ~ YearCollected*CI+offset(logTL_mm)+(1|CatalogNumber),family=nbinom1,data = PIMVIG_MONOdact)

AIC(glmmonodact,glmmonodact1)

summary(glmmonodact1)

#Evaluate residuals
library(DHARMa)
library(MASS)

par(mfrow = c(2, 2))
s=simulateResiduals(fittedModel=glmmonodact1,n=250)
s$scaledResiduals
plot(s)

## final plot with predictions and CI
library(ggeffects)
plot_model(glmmonodact1)

mydf <- ggpredict(glmmonodact1, c("YearCollected [n=100]", "CI"), jitter=TRUE) 

apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

plot(mydf, rawdata = TRUE, alpha = 0.3, dot.alpha = 0.6,colors=c("#969696","#ce1256"))+
  labs(x = 'Year', y = 'Monogenean abundance (# monogeneans/fish)',title=NULL)+
  apatheme

dev.off()


ggsave(file="REU/Imani J/Model plots/PIMVIG_Dactylogyrus_modelplot.png", width=150, height=90, dpi=1000, units = "mm")
ggsave(file="REU/Imani J/Model plots/PIMVIG_Dactylogyrus_modelplot.pdf", width=150, height=90, dpi=1000, units = "mm")



# Pimvig mono.gyro

PIMVIG_MONOgyro <- subset(PIMVIG, psite_spp.x == 'MONO.GYRO'|psite_spp.x == 'MONO.LG')

library(stats)
library(lme4)
library(DHARMa)
library(glmmTMB)


#Model passed diagnostics test
glmmonogyro1 <- glmmTMB(psite_count ~ YearCollected*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom2,data = PIMVIG_MONOgyro)

# Did not pass KS Test. The Kolmogorov-Smirnov Test is a type of non-parametric test of the equality of discontinuous and continuous a 1D probability distribution that is used to compare the sample with the reference probability test (known as one-sample K-S Test) or among two samples (known as two-sample K-S test). A K-S Test quantifies the distance between the cumulative distribution function of the given reference distribution and the empirical distributions of given two samples, or between the empirical distribution of given two samples.
glmmonogyro1 <- glmmTMB(psite_count ~ YearCollected*CI+TotalLength_mm+(1|CatalogNumber),family=poisson,data = PIMVIG_MONOgyro)

#Best one
glmmonogyro1 <- glmmTMB(psite_count ~ YearCollected*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom1,data = PIMVIG_MONOgyro)
glmmonogyro1 <- glmmTMB(psite_count ~ YearCollected*CI+offset(logTL_mm)+(1|CatalogNumber),family=nbinom1,data = PIMVIG_MONOgyro)

summary(glmmonogyro1)

#Evaluate residuals
library(DHARMa)
library(MASS)

par(mfrow = c(2, 2))
s=simulateResiduals(fittedModel=glmmonogyro1,n=250)
s$scaledResiduals
plot(s)

## final plot with predictions and CI
library(ggeffects)
plot_model(glmmonogyro1)

mydf <- ggpredict(glmmonogyro1, c("YearCollected [n=100]", "CI"), jitter=TRUE) 

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'))


plot(mydf, rawdata = TRUE, alpha = 0.3, dot.alpha = 0.6,colors=c("#969696","#ce1256"))+
  labs(x = 'Year', y = 'Monogenean abundance (# monogeneans/fish)',title=NULL)+
  apatheme


ggsave(file="REU/Imani J/Model plots/PIMVIG_Gyrodactylus_modelplot.png", width=150, height=90, dpi=1000, units = "mm")
ggsave(file="REU/Imani J/Model plots/PIMVIG_Gyrodactylus_modelplot.pdf", width=150, height=90, dpi=1000, units = "mm")



### Analaysis for Notropis atherinoides----

# geom point adds points to graph.
#labs=labels. title= title of graph, all words on graphs are quotes. X axis= label for x axis.
#y=label for y axis. color= represents the legend or key of labels.

library(readr)

#NOTATH <- read_csv("data/processed/Notropis_atherinoides_processed_machine_readable.csv")

NOTATH <- subset(full_dataset, Fish_sp.x=="Notropis atherinoides")

NOTATH_MONOS <- subset(NOTATH, psite_spp.x == "MONO.UNK"| psite_spp.x == "MONO.ALL")

library (ggplot2)

  apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

  ggplot(NOTATH_MONOS, aes(x= YearCollected,
                          y=psite_count,
                          ,color=CI,
                          group=CI))+
    geom_point()+ 
    labs(title = "Counts of monogenean species in Notrophis Anthernoides", x="Year", y="Monogenean abundance (# monogeneans/fish)",color= 
           "Pollution impact:")+apatheme+geom_vline(xintercept=1972, linetype="dashed", color = "black", size=0.5)+
    facet_wrap(~ psite_spp)
  

# Models
  
  library(stats)
  library(lme4)
  library(DHARMa)
  library(glmmTMB)
  
  #Is linear, no quadratic term. Good diagnostics
  glm_notathmono1 <- glmmTMB(psite_count ~ YearCollected*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom1,data = NOTATH_MONOS)
  glm_notathmono2 <- glmmTMB(psite_count ~ YearCollected*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom2,data = NOTATH_MONOS)
  
  #Model has good gianostics
  glm_notathmono3 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom1,data = NOTATH_MONOS)
  
  #Model has bad gianostics
  glm_notathmono2 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom2,data = NOTATH_MONOS)
  
  #Model selection. The lower the AIC, the better the model is.
  AIC(glm_notathmono1,glm_notathmono2,glm_notathmono3)
  
  #After evaluating AIC, this one is the best.
  glm_notathmono1 <- glmmTMB(psite_count ~ YearCollected*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom1,data = NOTATH_MONOS)
  
  glm_notathmono1 <- glmmTMB(psite_count ~ YearCollected*CI+offset(logTL_mm)+(1|CatalogNumber),family=nbinom1,data = NOTATH_MONOS)
  
  summary(glm_notathmono1)
  
  #Evaluate residuals
  library(DHARMa)
  library(MASS)
  
  par(mfrow = c(2, 2))
  s=simulateResiduals(fittedModel=glm_notathmono1,n=250)
  s$scaledResiduals
  plot(s)
  
  ## final plot with predictions and CI
  library(ggeffects)
  plot_model(glm_notathmono2)
  
  mydf <- ggpredict(glm_notathmono1, c("YearCollected [n=100]", "CI"), jitter=TRUE) 
  
  apatheme=theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.line=element_line(),
          text=element_text(family='Times'))
  
  
  plot(mydf, rawdata = TRUE, alpha = 0.3, dot.alpha = 0.6,colors=c("#969696","#ce1256"))+
    labs(x = 'Year', y = 'Monogenean abundance (# monogeneans/fish)',title=NULL)+
    apatheme
  
  
  ggsave(file="REU/Imani J/Model plots/NOTATH_DACTYLOGYRUS_modelplot.png", width=150, height=90, dpi=1000, units = "mm")
  ggsave(file="REU/Imani J/Model plots/NOTATH_DACTYLOGYRUS_modelplot.pdf", width=150, height=90, dpi=1000, units = "mm")
  
  

### Analaysis for Hybognathus nuchalis----
  
  # geom point adds points to graph.
  #labs=labels. title= title of graph, all words on graphs are quotes. X axis= label for x axis.
  #y=label for y axis. color= represents the legend or key of labels.
  
  library(readr)
  
  HYBNUC <- read_csv("data/processed/Hybognathus_nuchalis_processed_machine_readable.csv")
  HYBNUC_MONOS <- subset(HYBNUC, psite_spp == "MONO.DACT")
  
  library (ggplot2)
  
  apatheme= theme_bw(base_size = 11,base_family = "sans")+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.line=element_line())
  
  ggplot(HYBNUC_MONOS, aes(x= YearCollected,
                           y=psite_count,
                           ,color=CI,
                           group=CI))+
    geom_point()+ 
    labs(title = "Counts of monogenean species in Hybognathus nuchalis", x="Year", y="Monogenean abundance (# monogeneans/fish)",color= 
           "Pollution impact:")+apatheme+geom_vline(xintercept=1972, linetype="dashed", color = "black", size=0.5)+
    facet_wrap(~ psite_spp)
  
  
  # Models
  
  library(stats)
  library(lme4)
  library(DHARMa)
  library(glmmTMB)
  
  HYBNUC_MONOS$CI <- as.factor(HYBNUC_MONOS$CI)
  
  #Is linear, no quadratic term. Good diagnostics
  glm_hybnucmono1 <- glmmTMB(psite_count ~ YearCollected*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom1,data = HYBNUC_MONOS)

  #Is linear, no quadratic term. Bad diagnostics
    glm_hybnucmono2 <- glmmTMB(psite_count ~ YearCollected*CI+TotalLength_mm+(1|CatalogNumber),family=nbinom2,data = HYBNUC_MONOS)
  
  summary(glm_hybnucmono1)
  
  #Evaluate residuals
  library(DHARMa)
  library(MASS)
  
  par(mfrow = c(2, 2))
  s=simulateResiduals(fittedModel=glm_hybnucmono1,n=250)
  s$scaledResiduals
  plot(s)
  
  ## final plot with predictions and CI
  library(ggeffects)
  plot_model(glm_hybnucmono1)
  
  mydf <- ggpredict(glm_hybnucmono1, c("YearCollected [n=100]", "CI"), jitter=TRUE) 
  
  apatheme=theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.line=element_line(),
          text=element_text(family='Times'))
  
  
  plot(mydf, rawdata = TRUE, alpha = 0.3, dot.alpha = 0.6,colors=c("#969696","#ce1256"))+
    labs(x = 'Year', y = 'Monogenean abundance (# monogeneans/fish)',title=NULL)+
    apatheme
  
  
  ggsave(file="REU/Imani J/Model plots/NOTATH_DACTYLOGYRUS_modelplot.png", width=150, height=90, dpi=1000, units = "mm")
  ggsave(file="REU/Imani J/Model plots/NOTATH_DACTYLOGYRUS_modelplot.pdf", width=150, height=90, dpi=1000, units = "mm")
  
  
  