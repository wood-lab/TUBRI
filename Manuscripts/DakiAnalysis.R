


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

full_dataset <- read_csv("data/processed/Full_dataset_physical_2024.09.25.csv")

## Make a subset for myxozoans

#full_dataset_myxo <- subset(Full_dataset_elemN, Parasite_taxonomic_group == "Myxozoa")
full_dataset_myxo <- subset(full_dataset, Parasite_taxonomic_group == "Myxozoa")

## Create column for presence and absegroup_by()## Create column for presence and absence of myxozoans
full_dataset_myxo$psite_presence <- ifelse(full_dataset_myxo$psite_count > 0,                # condition
                                           1,    # what if condition is TRUE
                                           0)       # what if condition is FALSE
## Summary stats
prevalence_myxo <- subset(full_dataset_myxo, !is.na(psite_presence))
abundance_myxo <- subset(full_dataset_myxo, !is.na(psite_count))


summary_prevalence <- prevalence_myxo %>% group_by(Fish_sp.x,Parasite_genus) %>% 
  summarise(Prevalence=mean(psite_presence),
            .groups = 'drop') 

summary_abundance <- abundance_myxo %>% group_by(Fish_sp.x,Parasite_genus) %>% 
  summarise(Mean_abundance=mean(psite_count),
            Max_abundance=max(psite_count),
.groups = 'drop') 

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

## Create subset per fish species

pimvig_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Pimephales vigilax")
gamaff_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Gambusia affinis")
ictpun_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Ictalurus punctatus")
notath_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Notropis atherinoides")
carvel_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Carpiodes velifer")
hybnuc_myxo <- subset(full_dataset_myxo, Fish_sp.x == "Hybognathus nuchalis")


## Crease subset for myxozoans that can be counted

myxo_count <- subset(full_dataset_myxo, Parasite_genus == "Myxobolus"|Parasite_genus == "Henneguya"|Parasite_genus == "Unicauda"|Parasite_genus == "Thelohanellus")

# Scale variables

myxo_count$logTL_mm <- log(myxo_count$TotalLength_mm)
myxo_count$scaledmean_streamflow <- scale(myxo_count$mean_streamflow)
myxo_count$scaledmean_temperature <- scale(myxo_count$mean_temperature)
myxo_count$scaledYear <- scale(myxo_count$YearCollected)

## Create subset of myxo_count per fish species

pimvig_count <- subset(myxo_count, Fish_sp.x == "Pimephales vigilax")
gamaff_count <- subset(myxo_count, Fish_sp.x == "Gambusia affinis")
ictpun_count <- subset(myxo_count, Fish_sp.x == "Ictalurus punctatus")
notath_count <- subset(myxo_count, Fish_sp.x == "Notropis atherinoides")
carvel_count <- subset(myxo_count, Fish_sp.x == "Carpiodes velifer")
hybnuc_count <- subset(myxo_count, Fish_sp.x == "Hybognathus nuchalis")


## Restrict to year and control

# Restrict to the years that we are going to work with
myxo_count_y <- subset(myxo_count, 
                             YearCollected > 1972 &
                               YearCollected < 1995 )

# Restrict to the years that we are going to work with
myxo_count_yc <- subset(myxo_count_y, 
                                    CI == "control" )


### Add Nutrients
## Import data for nutrients

nutrients <- read_csv("data/Physicochemical/nutrients_resultphyschem.csv")

# Read USGS-gauge site metadata

siteinfo <-  read_csv("data/Physicochemical/station_info.csv")

# Add info on sampling sites (latitude, longitude, CI, etc.)

nutrients_withmetadata <- merge(nutrients, siteinfo, by.x = "MonitoringLocationIdentifier", by.y = "MonitoringLocationIdentifier", all.x = TRUE)

# Nitrogen species: nitrogen available to plants occurs in the form of 
# ammonia (NH3), ammonium (NH4+), nitrate (NO3−), nitrite (NO2-), and 
# organic nitrogen. Therefore, we used the mixed nitrogen measurement for our analyses.
# However it only spans between 1973 and 1994. So we will perform the analyses 
# only for that time period and for the control because I cannot assume that the 
# downstram impact sites are going to be represented by this data.


# Take only mixed forms of nitrogen in mg-N/L

mixN_unitall <- subset(nutrients_withmetadata, Measure == "Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3)")
mixN_year <- subset(mixN_unitall,  
                    YearCollected > 1972 &
                      YearCollected < 1995 )
mixN_control <- subset(mixN_year, CI == "control")
mixN <- subset(mixN_control, Unit_full == "mg/l")


#### The script below assigns to each fish individual the mean nutrients (mean_nitrogen) that the fish experienced over one year before its collection

# Add a new column to store the mean nutrient in the form of mixed nitrogen forms
myxo_count_yc$mean_nitrogen <- NA

library(lubridate)

# Combine year, month, and day columns into a date column
myxo_count_yc$collection_date <- make_date(myxo_count_yc$YearCollected, myxo_count_yc$MonthCollected, 
                                           myxo_count_yc$DayCollected)


mixN$measurement_date <- make_date(mixN$YearCollected, mixN$MonthCollected, 
                                   mixN$DayCollected)


# Loop through each row of 'myxo_count_yc'
for (i in 1:nrow(myxo_count_yc)) {
  # Extract the collection date for the individual
  collection_date <- as.Date(myxo_count_yc$collection_date[i])
  
  # Define the start and end date for one year before the collection date
  start_date <- collection_date - 365
  end_date <- collection_date
  
  # Filter 'mixN' to get the nitrogen data within that one-year range
  nutrient_data <- mixN[mixN$measurement_date >= start_date & mixN$measurement_date < end_date, "Result"]
  
  # Calculate the mean nitrogen and store it in 'myxo_count_yc'
  if (length(nutrient_data) > 0) {
    myxo_count_yc$mean_nitrogen[i] <- mean(nutrient_data, na.rm = TRUE)
  } else {
    myxo_count_yc$mean_nitrogen[i] <- NA  # If no data is found, store NA
  }
}

# View updated 'myxo_count_yc'
print(myxo_count_yc)

# rescale nutrients
myxo_count_yc$scaledmean_nitrogen <- scale(myxo_count_yc$mean_nitrogen)

### Elements
## Import data for metals

inorganics <- read_csv("data/Physicochemical/inorganics_resultphyschem.csv")

# Add info on sampling sites (latitude, longitude, CI, etc.)

inorganics_withmetadata <- merge(inorganics, siteinfo, by.x = "MonitoringLocationIdentifier", by.y = "MonitoringLocationIdentifier", all.x = TRUE)

# Summarize data table to year, day, and month collected
inorganics_means_ymd <- inorganics_withmetadata %>% group_by(YearCollected,MonthCollected,DayCollected,CI,Unit,Measure) %>% 
  summarise(Result=mean(Result),
            .groups = 'drop') 

# make the table wide
elements_na_yc_wide <- spread(inorganics_means_ymd, Measure,Result)

# take only the years and control sites
# Make a subset for element concentration in water

elements_water_w <- subset(elements_na_yc_wide, Unit=="ug/l"|
                             Unit=="mg/l")

# Restrict to the years that we are going to work with
elements_water_years_w <- subset(elements_water_w, 
                                 YearCollected > 1972 &
                                   YearCollected < 1995 )

# Restrict to the years that we are going to work with
elements_water_yearscontrolw <- subset(elements_water_years_w, 
                                       CI == "control" )

# Combine year, month, and day columns into a date column
elements_water_yearscontrolw$measurement_date <- make_date(elements_water_yearscontrolw$YearCollected, 
                                                           elements_water_yearscontrolw$MonthCollected, 
                                                           elements_water_yearscontrolw$DayCollected)

# Initialize new columns in myxo_count to store mean values for each variable

myxo_count_yc$As <- NA
myxo_count_yc$Ba <- NA
myxo_count_yc$Cd <- NA
myxo_count_yc$Cr <- NA
myxo_count_yc$Cu <- NA
myxo_count_yc$Pb <- NA
myxo_count_yc$Fe <- NA
myxo_count_yc$Mg <- NA
myxo_count_yc$Mn <- NA
myxo_count_yc$Hg <- NA
myxo_count_yc$Ni <- NA
myxo_count_yc$Zn <- NA

# Loop through each row of 'myxo_count_yc'
for (i in 1:nrow(myxo_count_yc)) {
  # Extract the collection date for the individual
  collection_date <- as.Date(myxo_count_yc$collection_date[i])
  
  # Define the start and end date for one year before the collection date
  start_date <- collection_date - 365
  end_date <- collection_date
  
  # Filter "elements_water_yearscontrolw" to get the data within that one-year range
  filtered_data <- elements_water_yearscontrolw[elements_water_yearscontrolw$measurement_date >= start_date & elements_water_yearscontrolw$measurement_date < end_date, ]
  
  # Calculate the mean for each variable and store it in 'myxo_count_yc'
  if (nrow(filtered_data) > 0) {
    myxo_count_yc$As[i] <- mean(filtered_data$Arsenic, na.rm = TRUE)
    myxo_count_yc$Ba[i] <- mean(filtered_data$Barium, na.rm = TRUE)
    myxo_count_yc$Cd[i] <- mean(filtered_data$Cadmium, na.rm = TRUE)
    myxo_count_yc$Cr[i] <- mean(filtered_data$Chromium, na.rm = TRUE)
    myxo_count_yc$Cu[i] <- mean(filtered_data$Copper, na.rm = TRUE)
    myxo_count_yc$Pb[i] <- mean(filtered_data$Lead, na.rm = TRUE)
    myxo_count_yc$Fe[i] <- mean(filtered_data$Iron, na.rm = TRUE)
    myxo_count_yc$Mg[i] <- mean(filtered_data$Magnesium, na.rm = TRUE)
    myxo_count_yc$Mn[i] <- mean(filtered_data$Manganese, na.rm = TRUE)
    myxo_count_yc$Hg[i] <- mean(filtered_data$Mercury, na.rm = TRUE)
    myxo_count_yc$Ni[i] <- mean(filtered_data$Nickel, na.rm = TRUE)
    myxo_count_yc$Zn[i] <- mean(filtered_data$Zinc, na.rm = TRUE)
    
  } else {
    myxo_count_yc$As[i] <- NA
    myxo_count_yc$Ba[i] <- NA
    myxo_count_yc$Cd[i] <- NA
    myxo_count_yc$Cr[i] <- NA
    myxo_count_yc$Cu[i] <- NA
    myxo_count_yc$Pb[i] <- NA
    myxo_count_yc$Fe[i] <- NA
    myxo_count_yc$Mg[i] <- NA
    myxo_count_yc$Mn[i] <- NA
    myxo_count_yc$Hg[i] <- NA
    myxo_count_yc$Ni[i] <- NA
    myxo_count_yc$Zn[i] <- NA
  }
}

## Remove non-numeric variables and elements which almost do not have data or have mostly zeroes (Beryllium, Chromium VI, Cobalt, Lithium, Molybdenum, Silver, Selenium )

elements_na_wide_num <- myxo_count_yc[-c(1:43)]
#elements_na_wide_num <- elements_na_wide[-c(1:4)] #new

elements_na_wide_scaled <- scale(elements_na_wide_num)

# Estimate the number of dimensions for the Principal Component Anal- ysis by cross-validation
# The number of components which leads to the smallest mean square error of prediction (MSEP) is retained. 
# For the Kfold cross-validation, pNA percentage of missing values is inserted and predicted with a PCA model using ncp.min to ncp.max dimensions.

library(missMDA)

nb <- estim_ncpPCA(elements_na_wide_scaled,method.cv = "Kfold", verbose = FALSE)
nb$ncp #5

plot(0:5, nb$criterion, xlab = "nb dim", ylab = "MSEP")

# Because we have a lot of NAs and PCAs can't handle that, we have to impute values based on our data
# The (regularized) iterative PCA algorithm (by default) first consists imputing missing values with initial values such as the mean of the variable.
# It is adviced to use the regularized version of the algorithm to avoid the overfitting problems which are very frequent when there are many missing values. 

library(FactoMineR)

res.comp <- imputePCA(elements_na_wide_scaled, ncp = nb$ncp) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set

imp <- cbind.data.frame(res.comp$completeObs)

res.pca <- PCA(imp, quanti.sup = 1, quali.sup = 12, ncp = 2, graph=FALSE) # generate PCA

plot(res.pca, choix="var") # visualize PCA

score <- as_tibble(factoextra::get_pca_ind(res.pca)$coord) #extract individual scores to be used in glm

myxo_count_ycm <- cbind(myxo_count_yc, score[1:2]) # merge scores with original data

## Create subset of myxo_count_ycm per fish species

pimvig_count_ycm <- subset(myxo_count_ycm, Fish_sp.x == "Pimephales vigilax")
gamaff_count_ycm <- subset(myxo_count_ycm, Fish_sp.x == "Gambusia affinis")
ictpun_count_ycm <- subset(myxo_count_ycm, Fish_sp.x == "Ictalurus punctatus")
notath_count_ycm <- subset(myxo_count_ycm, Fish_sp.x == "Notropis atherinoides")
carvel_count_ycm <- subset(myxo_count_ycm, Fish_sp.x == "Carpiodes velifer")
hybnuc_count_ycm <- subset(myxo_count_ycm, Fish_sp.x == "Hybognathus nuchalis")

# Set a theme for all your plots
apatheme= theme_bw(base_size = 14,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))


### Analysis of count data----

# Which random structure is best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID),
              data = carvel_count_ycm,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|CatalogNumber/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge


AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)


### Best link function

m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm,
              family = nbinom2(link="sqrt")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm,
              family = nbinom2(link="log")) # with log-link function, nbinom1 does not converge


AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m1)

# Flow with CI and psite_genus

mydf <- ggpredict(m1, terms= c("mean_nitrogen")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme

# Plot interactions

library(interactions)

interact_plot(m1, pred = mean_nitrogen, modx = mean_temperature,modx.values = c(18,20,22,24),
              allow.new.levels=TRUE,interval = TRUE,
              int.width = 0.95)+
  apatheme+  ggtitle("Full")


### Non-point source - Only for MYX.G----

carvel_count_ycm_myxg <- subset(carvel_count_ycm,psite_spp.x == "MYX.G")

# Which random structure is best? With lat_long instead of site
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|lat_long/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|lat_long/CatalogNumber/IndividualFishID),
              data = carvel_count_ycm_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|CatalogNumber/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|lat_long/CatalogNumber/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge


AIC(m1,m2,m3,m4)


# Which random structure is best?
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID),
              data = carvel_count_ycm_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|CatalogNumber/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge


AIC(m1,m2,m3,m4)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

### Which link function is best?

m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|lat_long/CatalogNumber/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|lat_long/CatalogNumber/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm_myxg,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|lat_long/CatalogNumber/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm_myxg,
              family = nbinom2(link="sqrt")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|lat_long/CatalogNumber/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm_myxg,
              family = nbinom2(link="log")) # with log-link function, nbinom1 does not converge

AIC(m1,m2,m3,m4)

# The best is 

m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)+
                offset(logTL_mm)+
                (1|lat_long/CatalogNumber/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)*Dim.1+
              scale(mean_streamflow)*scale(mean_nitrogen)*
                scale(mean_temperature)*Dim.2+
                offset(logTL_mm)+
                (1|lat_long/CatalogNumber/IndividualFishID)+
                (1|season),
              data = carvel_count_ycm_myxg,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m2)

# Flow with CI and psite_genus

mydf <- ggpredict(m1, terms= c("mean_nitrogen")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme

# Plot interactions

library(interactions)

interact_plot(m1, pred = mean_nitrogen, modx = mean_temperature,modx.values = c(18,20,22,24),
              allow.new.levels=TRUE,interval = TRUE,
              int.width = 0.95)+
  apatheme+  ggtitle("Full")


### Impact of point-source pollution ----



ggplot(carvel_count, aes(x= YearCollected,
                              y=psite_count,color=psite_spp.x))+
  geom_point()+apatheme+ylab("# myxozoan cysts/fish")+
  facet_wrap("Fish_sp.x")+
  ggtitle("Cyst number per fish")


m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|psite_spp.x)+
                (1|season),
              data = carvel_count,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|psite_spp.x)+
                (1|site/CatalogNumber/IndividualFishID),
              data = carvel_count,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID)+
                (1|psite_spp.x)+
                (1|season),
              data = carvel_count,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

AIC(m1,m2,m3) # best random structure according to AIC is m1: (1|site/IndividualFishID)+(1|season)

summary(m1)

#Evaluate residuals
s=simulateResiduals(fittedModel=m3,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

tab_model(m1)
plot_model(m3,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## Which is the best link function
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
          CI*scale(mean_temperature)*scale(mean_streamflow)+
          offset(logTL_mm)+
          (1|site/CatalogNumber/IndividualFishID)+
          (1|psite_spp.x)+
          (1|season),
        data = carvel_count,
        family = nbinom1(link="sqrt")) 

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID)+
                (1|psite_spp.x)+
                (1|season),
              data = carvel_count,
              family = nbinom1(link="log")) 

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID)+
                (1|psite_spp.x)+
                (1|season),
              data = carvel_count,
              family = nbinom2(link="log")) 

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID)+
                (1|psite_spp.x)+
                (1|season),
              data = carvel_count,
              family = nbinom2(link="sqrt")) 


AIC(m1,m2,m3,m4)

#best one
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID)+
                (1|psite_spp.x)+
                (1|season),
              data = carvel_count,
              family = nbinom2(link="log")) 

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m1, terms= c("YearCollected[n=100]","psite_spp.x")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+ylim(0,150)

mydf <- ggpredict(m1, terms= c("mean_streamflow[n=100]","CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme

mydf <- ggpredict(m1, terms= c("mean_temperature[n=100]","CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlim(18,24)+ylim(0,1200)

mydf <- ggpredict(m1, terms= c("mean_streamflow[n=100]","mean_temperature")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme

# Plot interactions

library(interactions)

interact_plot(m1, pred = mean_streamflow, modx = mean_temperature,
              allow.new.levels=TRUE,interval = FALSE,plot.points = FALSE,
              int.width = 0.95,partial.residuals = TRUE)+
  apatheme+  ggtitle("Full")

#modx.values = c(18,20,22,24),


### Impact of point-source pollution - MYX.G----

carvel_count_myxg <- subset(carvel_count,psite_spp.x == "MYX.G")

## Effect of time

m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = carvel_count_myxg,
              family = nbinom2(link="sqrt")) 

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

summary(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m1, terms= c("YearCollected[n=200]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlab("Year")+ylab("# of pseudocysts per fish")


### Testing the scaling of predictors -- another method
### This is another method that includes scaled the predictors beforehand. 
### However, it takes more coding effort and yields the same result as the method above.

carvel_count_myxg$scaledYear

## Effect of time

m1 <- glmmTMB(psite_count ~ scaledYear+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = carvel_count_myxg,
              family = nbinom2(link="sqrt")) 

summary(m1)

# Get the scaling parameters
mean_predictor <- attr(scale(carvel_count_myxg$YearCollected), "scaled:center")
sd_predictor <- attr(scale(carvel_count_myxg$YearCollected), "scaled:scale")

# Generate predictions with ggpredict
library(ggeffects)
pred <- ggpredict(m1, terms = "scaledYear")

# Back-transform to original scale
pred$YearCollected <- pred$x * sd_predictor + mean_predictor

# Plot the back-transformed predictions
library(ggplot2)

ggplot(pred, aes(x = YearCollected, y = predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(x = "Original Predictor", y = "Predicted Response")+ylim(0,1000)

###
m1 <- glmmTMB(psite_count ~ 
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected)+
                (1|season),
              data = carvel_count_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ 
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|YearCollected)+
              (1|season),
              data = carvel_count_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ 
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site)+
              (1|season),
              data = carvel_count_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ 
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected),
              data = carvel_count_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

AIC(m1,m2,m3,m4) # best random structure according to AIC is m1: (1|site/IndividualFishID)+(1|season)

summary(m1)

#Evaluate residuals
s=simulateResiduals(fittedModel=m4,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## Which is the best link function
m1 <- glmmTMB(psite_count ~ 
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected),
              data = carvel_count_myxg,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ 
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected),
              data = carvel_count_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ 
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected),
              data = carvel_count_myxg,
              family = nbinom2(link="log")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ 
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected),
              data = carvel_count_myxg,
              family = nbinom2(link="sqrt")) # with log-link function, nbinom1 does not converge


AIC(m1,m2,m3,m4)

#best one
m1 <- glmmTMB(psite_count ~ 
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site)+
                (1|YearCollected),
              data = carvel_count_myxg,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

### From Welicky et al. 2021 Front Ecol Evol

## Test for temporal autocorrelation
library(lmtest)
time <- unique(carvel_count_myxg$YearCollected)
dwtest(m1, order.by = time, alternative = "two.sided", exact = FALSE, tol = 1e-10)

## Test for spatial autocorrelation

#Jitter coordinates to remove repeats - required for autocorrelation analyses
carvel_count_myxg$latjitt<-jitter(carvel_count_myxg$Latitude, factor=0.1, amount=NULL)
carvel_count_myxg$longjitt<-jitter(carvel_count_myxg$Longitude, factor=0.1, amount=0)

simspatial.parasite<-simulateResiduals(fittedModel = m1)
testSpatialAutocorrelation(simulationOutput = simspatial.parasite,  x = carvel_count_myxg$longjitt, y = carvel_count_myxg$latjitt)


#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m2)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m1, terms= c("mean_streamflow[n=100]","CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlab("Streamflow (m3/sec)")+ylab("# of pseudocysts per fish")

mydf <- ggpredict(m1, terms= c("mean_temperature[n=100]","CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlab("Temperature (C)")+ylab("# of pseudocysts per fish")

mydf <- ggpredict(m1, terms= c("mean_streamflow[n=100]","mean_temperature"),type = 'fe') 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05)+
  apatheme


# Plot interactions

library(interactions)

interact_plot(m1, pred = mean_streamflow, modx = mean_temperature,
              allow.new.levels=TRUE,modx.values = c(18,20,22),interval = TRUE,plot.points = TRUE,
              int.width = 0.95,partial.residuals = FALSE)+
  apatheme


### Impact of point-source pollution - ICTPUNCT HENNEGUYA----

## Effect of time

m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site),
              data = ictpun_count,
              family = nbinom1(link="log")) 

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site)+
                (1|season),
              data = ictpun_count,
              family = nbinom1(link="log")) 

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|season),
              data = ictpun_count,
              family = nbinom1(link="log")) 

AIC(m1,m2,m3)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)


m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site),
              data = ictpun_count,
              family = nbinom2(link="sqrt")) 

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site),
              data = ictpun_count,
              family = nbinom2(link="log")) 

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site),
              data = ictpun_count,
              family = nbinom1(link="log")) 

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site),
              data = ictpun_count,
              family = nbinom1(link="sqrt")) 

AIC(m1,m2,m3,m4)

## 
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                offset(logTL_mm)+
                (1|site),
              data = ictpun_count,
              family = nbinom2(link="log")) 

performance::check_overdispersion(m1)

summary(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m1, terms= c("YearCollected[n=100]")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlab("Year")+ylab("# of pseudocysts per fish")



###
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = ictpun_count,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID),
              data = ictpun_count,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID)+
                (1|season),
              data = ictpun_count,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

AIC(m1,m2,m3) # best random structure according to AIC is m1: (1|site/IndividualFishID)+(1|season)

summary(m1)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## Which is the best link function
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = ictpun_count,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = ictpun_count,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = ictpun_count,
              family = nbinom2(link="sqrt")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = ictpun_count,
              family = nbinom2(link="log")) # with log-link function, nbinom1 does not converge


AIC(m1,m2,m3,m4)

#best one
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = ictpun_count,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge


### From Welicky et al. 2021 Front Ecol Evol

## Test for temporal autocorrelation
library(lmtest)
time <- unique(ictpun_count$YearCollected)
dwtest(m1, order.by = time, alternative = "two.sided", exact = FALSE, tol = 1e-10)

## Test for spatial autocorrelation

#Jitter coordinates to remove repeats - required for autocorrelation analyses
carvel_count_myxg$latjitt<-jitter(carvel_count_myxg$Latitude, factor=0.1, amount=NULL)
carvel_count_myxg$longjitt<-jitter(carvel_count_myxg$Longitude, factor=0.1, amount=0)

simspatial.parasite<-simulateResiduals(fittedModel = m1)
testSpatialAutocorrelation(simulationOutput = simspatial.parasite,  x = carvel_count_myxg$longjitt, y = carvel_count_myxg$latjitt)


#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m1)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m1, terms= c("YearCollected[all]","CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme

mydf <- ggpredict(m1, terms= c("mean_streamflow[n=100]","CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlim(100,450)

mydf <- ggpredict(m1, terms= c("mean_temperature[n=100]","CI")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme

mydf <- ggpredict(m1, terms= c("CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme

# Plot interactions

library(interactions)

interact_plot(m1, pred = mean_streamflow, modx = mean_temperature,
              allow.new.levels=TRUE,interval = FALSE,plot.points = FALSE,
              int.width = 0.95,partial.residuals = TRUE)+
  apatheme+  ggtitle("Full")

#modx.values = c(18,20,22,24),


### Impact of point-source pollution - NOTATH MYXOBOLUS----


###
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = notath_count,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID),
              data = notath_count,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/CatalogNumber/IndividualFishID)+
                (1|season),
              data = notath_count,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

AIC(m1,m2,m3) # best random structure according to AIC is m1: (1|site/IndividualFishID)+(1|season)

summary(m1)

#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

## Which is the best link function
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = notath_count,
              family = nbinom1(link="sqrt")) # with log-link function, nbinom1 does not converge

m2 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = notath_count,
              family = nbinom1(link="log")) # with log-link function, nbinom1 does not converge

m3 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = notath_count,
              family = nbinom2(link="sqrt")) # with log-link function, nbinom1 does not converge

m4 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = notath_count,
              family = nbinom2(link="log")) # with log-link function, nbinom1 does not converge


AIC(m1,m2,m3,m4)

#best one
m1 <- glmmTMB(psite_count ~ scale(YearCollected)+
                CI*scale(mean_temperature)*scale(mean_streamflow)+
                offset(logTL_mm)+
                (1|site/IndividualFishID)+
                (1|season),
              data = notath_count,
              family = nbinom2(link="log"))  # with log-link function, nbinom1 does not converge


### From Welicky et al. 2021 Front Ecol Evol

## Test for temporal autocorrelation
library(lmtest)
time <- unique(ictpun_count$YearCollected)
dwtest(m1, order.by = time, alternative = "two.sided", exact = FALSE, tol = 1e-10)

## Test for spatial autocorrelation

#Jitter coordinates to remove repeats - required for autocorrelation analyses
carvel_count_myxg$latjitt<-jitter(carvel_count_myxg$Latitude, factor=0.1, amount=NULL)
carvel_count_myxg$longjitt<-jitter(carvel_count_myxg$Longitude, factor=0.1, amount=0)

simspatial.parasite<-simulateResiduals(fittedModel = m1)
testSpatialAutocorrelation(simulationOutput = simspatial.parasite,  x = carvel_count_myxg$longjitt, y = carvel_count_myxg$latjitt)


#Evaluate residuals
s=simulateResiduals(fittedModel=m1,n=250)
s$scaledResiduals
plot(s)

performance::check_overdispersion(m1)

tab_model(m1)
plot_model(m1,type = "est")+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)
r.squaredGLMM(m1)

## visualize model

# Flow with CI and psite_genus

mydf <- ggpredict(m1, terms= c("YearCollected[all]")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme

mydf <- ggpredict(m1, terms= c("mean_streamflow[n=100]","CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme+xlim(100,450)

mydf <- ggpredict(m1, terms= c("mean_temperature[n=100]","CI")) 

plot(mydf,show_data=TRUE,show_residuals=FALSE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme

mydf <- ggpredict(m1, terms= c("CI")) 

plot(mydf,show_data=FALSE,show_residuals=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme

# Plot interactions

library(interactions)

interact_plot(m1, pred = mean_streamflow, modx = mean_temperature,
              allow.new.levels=TRUE,interval = FALSE,plot.points = FALSE,
              int.width = 0.95,partial.residuals = TRUE)+
  apatheme+  ggtitle("Full")

#modx.values = c(18,20,22,24),

