#This script is  ade analyze the monogenean data from the CAREER REU program
#Revised on July 17, 2024
#Written by Imani Jones and Daki Diaz


# If you need to change your working directory, paste the file structure between the quotation marks below.
# To get the file structure, open the folder and then "Get info" about that folder. Copy the file structure
# and paste it below.

setwd("C:/Users/imani/OneDrive - Tuskegee University/Desktop/TUBRI_Monogenea_Project/TUBRI")



## Stats with Jolee----
library(readr)

ICTPUNCT <- read_csv("data/processed/Ictalurus_punctatus_processed_machine_readable_UPDATED_2024.07.10.csv")


monogenes.ict<-ICTPUNCT[, c ("MONO.IP","MONO.UNK","YearCollected")]
monogenes.ict
#install.packages("tidyr")
library(tidyr)
#Label of data set<-command in r to change wide format to long format. The parentheses (monogenes.ict) I am telling it to switch from wide to long.
#I want the code to pull out columns with the name "Mono". Next line of code is being called "Species", the values thate in the one column are not being called "count"
monogenes.ict_table<-pivot_longer(monogenes.ict,
                                  cols=starts_with("MONO"),
                                  names_to = "Species",
                                  values_to = "count")

install.packages("ggplot2")
library(ggplot2)
#Data frame is a command in r, I want year to equal a list within the data set of "Years Collected"
#I want my monogenes to be called from the data set called "Species"
#I want my monogenes to come from the data set monogenes.ict column count
Monogenes.ict_table<- data.frame(Year=c (monogenes.ict_table$YearCollected),
                                 Monogenes=c(monogenes.ict_table$Species),
                                 monogenes_affected= c(monogenes.ict_table$count)
                                             )
#gg plot is plotting function in r, it allows multiple variables within the plot.
#first part of code is showing what data set I want to be called from. "aes"-- command in r, x is the variable in my data set, y is the monogenes affected in this data set. 
#Color is showing my third variable in data aset. "+"- indicates adding aesthic to graph. 
ggplot(Monogenes.ict_table, aes(x= Year,
                                y=monogenes_affected
                                ,color=Monogenes,
                                group=Monogenes))+
  geom_point()+ 
  labs(title = "Counts of Monogenes over the Years", x="Year", y="Monogenean abundance (# monogeneans/fish)",color= 
        "species")
# geom point adds points to graph.
#labs=labels. title= title of graph, all words on graphs are quotes. X axis= label for x axis.
#y=label for y axis. color= represents the legend or key of labels.

### Code for monogenean MONO.IP in ICTPUNct abundance across time and pollution impact----
# Once your working directory is set, you're ready to read in the data!  If there are sub-folders inside your 
# working directory, you'll need to specify them as I have below.
# The <- command tells R what a thing is called.  So you can read the line below as,
# Look at this csv file, and name it pim_vig_data.  When I call pim_vig_data, I want you to give me the csv file.

#pim_vig_data<-read.csv("data/processed/Ictalurus_punctatus_processed_human_readable.csv")
#pim_vig_data<-read.csv("Pimephales_vigilax_processed_machine_readable")

library(readr)

ICTPUNCT <- read_csv("data/processed/Ictalurus_punctatus_processed_machine_readable_UPDATED_2024.07.10.csv")

# To see your data as a spreadsheet in a new tab, use the View command.

 #View(ICTPUNCT)

# You can also see your data in the console below by running the name of the dataset.

ICTPUNCT

#We made a subset based on our data to include only MONO.IP

ICTPUNCT_MONOIP <- subset(ICTPUNCT, psite_spp == 'MONO.IP')

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
  labs(title = "Counts of Monogeneans over the Years", x="Year", y="Monogenean abundance (# monogeneans/fish)",color= 
         "Pollution impact:")+apatheme+geom_vline(xintercept=1972, linetype="dashed", color = "black", size=0.5)

# geom point adds points to graph.
#labs=labels. title= title of graph, all words on graphs are quotes. X axis= label for x axis.
#y=label for y axis. color= represents the legend or key of labels.

library(stats)
library(lme4)
library(DHARMa)
library(glmmTMB)

#Is linear, no quadratic term. The diagnostics are terrible.
glmmonoip1 <- glmmTMB(psite_count ~ YearCollected*CI+(1|CatalogNumber),family=nbinom1,data = ICTPUNCT_MONOIP)

#Models have good gianostics
glmmonoip1 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+(1|CatalogNumber),family=nbinom1,data = ICTPUNCT_MONOIP)
glmmonoip2 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+(1|CatalogNumber),family=nbinom2,data = ICTPUNCT_MONOIP)

#Model selection. The lower the AIC, the better the model is.
AIC(glmmonoip1,glmmonoip2)

#After evaluating AIC, this one is the best.
glmmonoip1 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+(1|CatalogNumber),family=nbinom1,data = ICTPUNCT_MONOIP)


summary(glmmonodact1)

#Evaluate residuals
library(DHARMa)
library(MASS)

par(mfrow = c(2, 2))
s=simulateResiduals(fittedModel=glmmonoip1,n=250)
s$scaledResiduals
plot(s)

#Generate table

library(sjPlot)
library(sjmisc)
library(sjlabelled)
tab_model(glmmonoip1, show.df = FALSE, show.aic = TRUE) #p.style = "a" #for asterics


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


plot(mydf, rawdata = TRUE, alpha = 0.3, dot.alpha = 0.6)+
  labs(x = 'Year', y = 'Ligictaluridus pricei abundance')+
  apatheme



### Code for monogenean MONO.DACT in PIMVIG abundance across time and pollution impact----
# Once your working directory is set, you're ready to read in the data!  If there are sub-folders inside your 
# working directory, you'll need to specify them as I have below.
# The <- command tells R what a thing is called.  So you can read the line below as,
# Look at this csv file, and name it pim_vig_data.  When I call pim_vig_data, I want you to give me the csv file.

#pim_vig_data<-read.csv("data/processed/Ictalurus_punctatus_processed_human_readable.csv")
#pim_vig_data<-read.csv("Pimephales_vigilax_processed_machine_readable")

library(readr)

PIMVIG <- read_csv("data/processed/Pimephales_vigilax_processed_machine_readable_UPDATED_2024.07.10.csv")

# To see your data as a spreadsheet in a new tab, use the View command.

# View(ICTPUNCT)

# You can also see your data in the console below by running the name of the dataset.





PIMVIG

#We made a subset based on our data to include only MONO.IP

PIMVIG_MONO <- subset(PIMVIG, psite_spp == 'MONO.DACT'|psite_spp == 'MONO.GYRO'|psite_spp == 'MONO.LG')

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


# Pimvig mono.dact

PIMVIG_MONOdact <- subset(PIMVIG, psite_spp == 'MONO.DACT')


library(stats)
library(lme4)
library(DHARMa)
library(glmmTMB)

glmmonodact1 <- glm(psite_count ~ YearCollected*CI,family=nbinom1, data = PIMVIG_MONOdact)
glmmonodact1 <- glmmTMB(psite_count ~ YearCollected*CI+(1|CatalogNumber),family=nbinom1,data = PIMVIG_MONOdact)
glmmonodact1 <- glmmTMB(psite_count ~ YearCollected*CI+(1|CatalogNumber),family=poisson,data = PIMVIG_MONOdact)

#Best one
glmmonodact1 <- glmmTMB(psite_count ~ YearCollected*CI+(1|CatalogNumber),family=nbinom1,data = PIMVIG_MONOdact)


summary(glmmonodact1)

#Evaluate residuals
library(DHARMa)
library(MASS)

par(mfrow = c(2, 2))
s=simulateResiduals(fittedModel=glmmonodact1,n=250)
s$scaledResiduals
plot(s)

#Generate table

library(sjPlot)
library(sjmisc)
library(sjlabelled)
tab_model(glmmonodact1, show.df = FALSE, show.aic = TRUE) #p.style = "a" #for asterics


## final plot with predictions and CI
library(ggeffects)
plot_model(glmmonodact1)

mydf <- ggpredict(glmmonodact1, c("YearCollected [n=100]", "CI"), jitter=TRUE) 

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'))


plot(mydf, rawdata = TRUE, alpha = 0.3, dot.alpha = 0.6)+
  labs(x = 'Year', y = 'Monogenean abundance')+
  apatheme


# Pimvig mono.gyro

PIMVIG_MONOgyro <- subset(PIMVIG, psite_spp == 'MONO.GYRO')

library(stats)
library(lme4)
library(DHARMa)
library(glmmTMB)


glmmonodact1 <- glm(psite_count ~ YearCollected*CI,family=nbinom1, data = PIMVIG_MONOgyro)
glmmonodact1 <- glm(psite_count ~ YearCollected*CI,family=quasipoisson, data = PIMVIG_MONOgyro)



#Evaluate residuals
library(DHARMa)
library(MASS)

s=stdres(glmmonodact1,arg="response")
plot(s)

par(mfrow = c(2, 2))
s=simulateResiduals(fittedModel=glmmonogyro1,n=250)
s$scaledResiduals
plot(s)

diagnostics(glmmsurv, ask=FALSE)

#Generate table

library(sjPlot)
library(sjmisc)
library(sjlabelled)
tab_model(glmmonodact1, show.df = FALSE, show.aic = TRUE) #p.style = "a" #for asterics


## final plot with predictions and CI
library(ggeffects)
plot_model(glmmonodact1)

mydf <- ggpredict(glmmonodact1, c("YearCollected [n=100]", "CI"), jitter=TRUE) 

apatheme=theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        text=element_text(family='Times'))


plot(mydf, rawdata = TRUE, alpha = 0.3, dot.alpha = 0.6)+
  labs(x = 'Year', y = 'Dactylogyrus abundance')+
  apatheme



# Analaysis for Notropis atherinoides----

# geom point adds points to graph.
#labs=labels. title= title of graph, all words on graphs are quotes. X axis= label for x axis.
#y=label for y axis. color= represents the legend or key of labels.

library(readr)
NOTATH <- read_csv("data/processed/Notropis_atherinoides_processed_machine_readable.csv")
View(NOTATH)

library(readr)

NOTATH <- read_csv("data/processed/Notropis_atherinoides_processed_machine_readable.csv")
NOTATH_MONOS <- subset(NOTATH, psite_spp == "MONO.UNK"| psite_spp == "MONO.ALL")

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
  
  #Is linear, no quadratic term. The diagnostics are terrible.
  glm_notathmono1 <- glmmTMB(psite_count ~ YearCollected*CI+(1|CatalogNumber),family=nbinom1,data = NOTATH_MONOS)
  glm_notathmono2 <- glmmTMB(psite_count ~ YearCollected*CI+(1|CatalogNumber),family=nbinom2,data = NOTATH_MONOS)
  
  #Models have good gianostics
  glm_notathmono1 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+(1|CatalogNumber),family=nbinom1,data = NOTATH_MONOS)
  glm_notathmono2 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+(1|CatalogNumber),family=nbinom2,data = NOTATH_MONOS)
  
  #Model selection. The lower the AIC, the better the model is.
  AIC(glm_notathmono1,glm_notathmono2)
  
  #After evaluating AIC, this one is the best.
  glm_notathmono1 <- glmmTMB(psite_count ~ poly(YearCollected,2)*CI+(1|CatalogNumber),family=nbinom1,data = ICTPUNCT_MONOIP)
  
  
  summary(glmmonodact1)
  
  #Evaluate residuals
  library(DHARMa)
  library(MASS)
  
  par(mfrow = c(2, 2))
  s=simulateResiduals(fittedModel=glm_notathmono2,n=250)
  s$scaledResiduals
  plot(s)
  
  #Generate table
  
  library(sjPlot)
  library(sjmisc)
  library(sjlabelled)
  tab_model(glm_notathmono2, show.df = FALSE, show.aic = TRUE) #p.style = "a" #for asterics
  
  
  ## final plot with predictions and CI
  library(ggeffects)
  plot_model(glm_notathmono2)
  
  mydf <- ggpredict(glm_notathmono2, c("YearCollected [n=100]", "CI"), jitter=TRUE) 
  
  apatheme=theme_bw()+
    theme(panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          panel.border=element_blank(),
          axis.line=element_line(),
          text=element_text(family='Times'))
  
  
  plot(mydf, rawdata = TRUE, alpha = 0.3, dot.alpha = 0.6)+
    labs(x = 'Year', y = 'Monogenean abundance')+
    apatheme
  
  
