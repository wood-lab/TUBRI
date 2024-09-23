## This script curates chemical data from USGS including nutrients and element concentrations
# Daki Diaz-Morales
# dmdiaz@uw.edu; diazdakeishla@gmail.com


### Nutrients: nitrate and nitrite----
## Import data for nutrients

nutrients <- read_csv("data/Physicochemical/nutrients_resultphyschem.csv")

# Read USGS-gage site metadata

siteinfo <-  read_csv("data/Physicochemical/station_info.csv")

# Add info on sampling sites (latitude, longitude, CI, etc.)

nutrients_withmetadata <- merge(nutrients, siteinfo, by.x = "MonitoringLocationIdentifier", by.y = "MonitoringLocationIdentifier", all.x = TRUE)

# Summarize by taking the mean per year

nutrients_means <- nutrients_withmetadata %>% group_by(YearCollected,CI,Measure) %>% 
  summarise(Result=mean(Result,na.rm=TRUE),
            .groups = 'drop')

# Nitrogen species: nitrogen available to plants occurs in the form of 
# ammonia (NH3), ammonium (NH4+), nitrate (NO3−), nitrite (NO2-), and 
# organic nitrogen. However, the data that we have most complete is that of inorganic N
# (nitrate and nitrite). Nitrate and nitrite have data for most years except between 1977-1984. 
# However, we have data for that window of time as "inorganic nitrogen (nitrat and nitrite).
# Therefore, we use that to fill the gap. 


# Take only nitrate, nitrite, and their sum (inorganic nitrogen (nitrate and nitrite))

nitrate_nitrite <- subset(nutrients_means, Measure == "Nitrite" |
                            Measure == "Nitrate"|
                            Measure == "Inorganic nitrogen (nitrate and nitrite)")


# From the following plot you can see that nitrate and nitrite have the longest time span
# Inorganic nitrogen was measure in most recent years. So we will the gap in nitrate and nitrite with inorganic N.

# Evaluate dynamic of nutrients across years

ggplot(nutrients_means, aes(x= as.factor(YearCollected),
                            y=Result, color=CI))+
  geom_point()+apatheme+
  ggtitle("Nutrients per year")+
  facet_wrap("Measure")



# We remove data before 1977 and after 1984 to not be reduntant on nitrate and nitrite concentatrions individually measured for those years.
inorgN_a <- nitrate_nitrite %>% filter(!(Measure == "Inorganic nitrogen (nitrate and nitrite)" &
                                               YearCollected < 1977))

inorganicN <- inorgN_a %>% filter(!(Measure == "Inorganic nitrogen (nitrate and nitrite)" &
                                           YearCollected > 1984))

# OK, now we have to still sum individual measurements of nitrate and nitrite to have total inorganic (nitrate and nitrite)

inorganicN$Measure <- as.factor(inorganicN$Measure)

inorganicN <- inorganicN %>% 
  mutate(Measure_type = case_when(Measure == "Nitrate" ~ "nitrate_nitrite",
                                Measure == "Nitrite" ~ "nitrate_nitrite",
                              Measure == "Inorganic nitrogen (nitrate and nitrite)"~ 
                              "inorg_N"))


inorganicN_sum <- inorganicN %>% group_by(YearCollected,CI,Measure_type) %>% 
  summarise(Result=sum(Result,na.rm=TRUE),
            .groups = 'drop')


# Re-evaluate dynamic of nitrates and nitrites across years

ggplot(inorganicN_sum, aes(x= as.factor(YearCollected),
                           y=Result, color=CI))+
  geom_point()+apatheme+
  ggtitle("Total nitrate and nitrite concentration per year")

### Phosphates did not change much over time so we keep it out to reduce dimensionality
### Nutrients: mixed forms of N----
## Import data for nutrients

nutrients <- read_csv("data/Physicochemical/nutrients_resultphyschem.csv")

# Read USGS-gage site metadata

siteinfo <-  read_csv("data/Physicochemical/station_info.csv")

# Add info on sampling sites (latitude, longitude, CI, etc.)

nutrients_withmetadata <- merge(nutrients, siteinfo, by.x = "MonitoringLocationIdentifier", by.y = "MonitoringLocationIdentifier", all.x = TRUE)

# Summarize by taking the mean per year

nutrients_means <- nutrients_withmetadata %>% group_by(YearCollected,CI,Measure) %>% 
  summarise(meanN=mean(Result,na.rm=TRUE),
            .groups = 'drop')

# Nitrogen species: nitrogen available to plants occurs in the form of 
# ammonia (NH3), ammonium (NH4+), nitrate (NO3−), nitrite (NO2-), and 
# organic nitrogen. Therefore, we used the mixed nitrogen measurement for our analyses.
# However it only spans between 1973 and 1994. So we will perform the analyses 
# only for that time period and for the control because I cannot assume that the 
# downstram impact sites are going to be represented by this data.


# Take only nitrate, nitrite, and their sum (inorganic nitrogen (nitrate and nitrite))

mixN <- subset(nutrients_means, Measure == "Nitrogen, mixed forms (NH3), (NH4), organic, (NO2) and (NO3)" &
                                YearCollected > 1972 &
                                YearCollected < 1995 )

# Restrict only to control

mixN_control <- subset(mixN, CI == "control")

# Evaluate dynamic of nutrients across years

ggplot(nutrients_means, aes(x= as.factor(YearCollected),
                            y=Result, color=CI))+
  geom_point()+apatheme+
  ggtitle("Nutrients per year")+
  facet_wrap("Measure")

### Phosphates did not change much over time so we keep it out to reduce dimensionality
### Elements---- 

## Import data for metals

inorganics <- read_csv("data/Physicochemical/inorganics_resultphyschem.csv")

# Add info on sampling sites (latitude, longitude, CI, etc.)

inorganics_withmetadata <- merge(inorganics, siteinfo, by.x = "MonitoringLocationIdentifier", by.y = "MonitoringLocationIdentifier", all.x = TRUE)

# Summarize by taking the mean per year

inorganics_means <- inorganics_withmetadata %>% group_by(YearCollected,CI,Unit,Measure) %>% 
  summarise(Result=mean(Result),
            .groups = 'drop') # We do this prior to the PCA because we are interested to see what happens with the metals year by year


inorganics_means <- inorganics_withmetadata %>% group_by(YearCollected,MonthCollected,CI,Unit,Measure) %>% 
  summarise(Result=mean(Result),
            .groups = 'drop') # new, summarizing by month to have more data

# Make a subset for element concentration in water

elements_water <- subset(inorganics_means, Unit=="ug/l"|
                           Unit=="mg/l")

# Restrict to the years that we are going to work with
elements_water_years <- subset(elements_water, 
                                 YearCollected > 1972 &
                                 YearCollected < 1995 )

# Restrict to the years that we are going to work with
elements_water_yearscontrol <- subset(elements_water_years, 
                               CI == "control" )

elements_water_control <- subset(elements_water, 
                                      CI == "control" )#new

## Analyse with PCA

# Remove NAs
elements_na <- subset(elements_water_yearscontrol, !is.na(Result)) # remove NAs
elements_na <- subset(elements_water_control, !is.na(Result)) # new
elements_na <- subset(elements_water, !is.na(Result)) # new

# Keep those elements that we are interested on and that have enough data to be analyze

elements_sub <- subset(elements_na, Measure %in% c("Aluminum","Arsenic","Barium",
                                   "Cadmium","Chloride","Chromium",
                                   "Copper","Lead","Iron","Magnesium",
                                   "Manganese","Mercury",
                                   "Nickel","Potassium","Sodium","Silica",
                                   "Zinc"))

elements_sub <- subset(elements_na, Measure %in% c("Arsenic","Barium",
                                                   "Cadmium","Chromium",
                                                   "Copper","Lead","Iron","Magnesium",
                                                   "Manganese","Mercury",
                                                   "Nickel",
                                                   "Zinc"))


# Make data in wide format
elements_na_wide <- spread(elements_sub, Measure,Result)

## Remove non-numeric variables and elements which almost do not have data or have mostly zeroes (Beryllium, Chromium VI, Cobalt, Lithium, Molybdenum, Silver, Selenium )

elements_na_wide_num <- elements_na_wide[-c(1:3)]
elements_na_wide_num <- elements_na_wide[-c(1:4)] #new

elements_na_wide_scaled <- scale(elements_na_wide_num)


# Estimate the number of dimensions for the Principal Component Anal- ysis by cross-validation
# The number of components which leads to the smallest mean square error of prediction (MSEP) is retained. 
# For the Kfold cross-validation, pNA percentage of missing values is inserted and predicted with a PCA model using ncp.min to ncp.max dimensions.

library(missMDA)

nb <- estim_ncpPCA(elements_na_wide_scaled,method.cv = "Kfold", verbose = FALSE)
nb$ncp #1

plot(0:5, nb$criterion, xlab = "nb dim", ylab = "MSEP")

# Because we have a lot of NAs and PCAs can't handle that, we have to impute values based on our data
# The (regularized) iterative PCA algorithm (by default) first consists imputing missing values with initial values such as the mean of the variable.
# It is adviced to use the regularized version of the algorithm to avoid the overfitting problems which are very frequent when there are many missing values. 

library(FactoMineR)

res.comp <- imputePCA(elements_na_wide_scaled, ncp = nb$ncp) # iterativePCA algorithm
res.comp$completeObs[1:3,] # the imputed data set

mydata.cor = cor(res.comp$completeObs[1:3,])

library(corrplot)
corrplot(mydata.cor)


imp <- cbind.data.frame(res.comp$completeObs)

res.pca <- PCA(imp, quanti.sup = 1, quali.sup = 12, ncp = 2, graph=FALSE) # generate PCA

plot(res.pca, choix="var") # visualize PCA

score <- as_tibble(factoextra::get_pca_ind(res.pca)$coord) #extract individual scores to be used in glm

mod <- cbind(elements_na_wide[1], score[1:2]) # merge scores with original data

mod 

### Merge with your data----

## Import parasite dataset

full_dataset <- read_csv("data/processed/Full_dataset_physical_2024.09.16.csv")
full_dataset <- read_csv("data/processed/Full_dataset_physical_2024.09.23.csv")

## Restrict to year and control

# Restrict to the years that we are going to work with
full_dataset_years <- subset(full_dataset, 
                               YearCollected > 1972 &
                                 YearCollected < 1995 )

# Restrict to the years that we are going to work with
full_dataset_yearscontrol <- subset(full_dataset_years, 
                                      CI == "control" )

# Merge with parasite data

Full_dataset_elements <- merge(full_dataset, mod, by.x = "YearCollected", by.y = "YearCollected", all.x = TRUE)
Full_dataset_elements <- merge(full_dataset_yearscontrol, mod, by.x = "YearCollected", by.y = "YearCollected", all.x = TRUE)


# Merge with nitrogen data

Full_dataset_elemN <- merge(Full_dataset_elements, mixN_control, by.x = "YearCollected", by.y = "YearCollected", all.x = TRUE)

## Make a subset for myxozoans

full_dataset_myxo <- subset(Full_dataset_elemN, Parasite_taxonomic_group == "Myxozoa")
full_dataset_myxo <- subset(Full_dataset_elements, Parasite_taxonomic_group == "Myxozoa")

## Create column for presence and absence of myxozoans
full_dataset_myxo$psite_presence <- ifelse(full_dataset_myxo$psite_count > 0,                # condition
                                           1,    # what if condition is TRUE
                                           0)       # what if condition is FALSE

## Remove columns with myxozoans that aren't confirmed to be myxozoans and therefore are not identified
full_dataset_myxo <- subset(full_dataset_myxo, !is.na(Parasite_genus))
full_dataset_myxo <- subset(full_dataset_myxo, psite_spp.x!="MYX.GM") # According to Stephen, this is likely a contamination and also it is in the gall bladder so does not count for abundace

#Revise variable types
full_dataset_myxo$IndividualFishID <- as.factor(full_dataset_myxo$IndividualFishID)
full_dataset_myxo$CI <- as.factor(full_dataset_myxo$CI)
full_dataset_myxo$Fish_sp.x <- as.factor(full_dataset_myxo$Fish_sp.x)
full_dataset_myxo$Parasite_genus <- as.factor(full_dataset_myxo$Parasite_genus)
full_dataset_myxo$before_after <- as.factor(full_dataset_myxo$before_after)
full_dataset_myxo$MonthCollected <- as.factor(full_dataset_myxo$MonthCollected)
full_dataset_myxo$season <- as.factor(full_dataset_myxo$season)

## Revise reference levels

#full_dataset_myxo$before_after <- relevel(full_dataset_myxo$before_after, ref = "before")
#full_dataset_myxo$CI <- relevel(full_dataset_myxo$before_after, ref = "control")


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

# Scale variables

myxo_count$logTL_mm <- log(myxo_count$TotalLength_mm)

myxo_count$scaledmeanFlow <- scale(myxo_count$meanFlow)

myxo_count$scaledmeanN <- scale(myxo_count$meanN)

myxo_count$scaledmeanTemp <- scale(myxo_count$mean_temperature)

myxo_count$scaledYear <- scale(myxo_count$YearCollected)

# Scale variables

myxo_count$logTL_mm <- log(myxo_count$TotalLength_mm)

myxo_count$scaledmean_streamflow <- scale(myxo_count$mean_streamflow)

myxo_count$scaledmean_nitrogen <- scale(myxo_count$mean_nitrogen)

myxo_count$scaledmean_temperature <- scale(myxo_count$mean_temperature)

myxo_count$scaledYear <- scale(myxo_count$YearCollected)



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


### Impact of multiple stressors on the abundance of myxobolus----


# Evaluate dynamic of nutrients across years

ggplot(myxo_count, aes(x= mean_temperature,
                            y=psite_count))+
  geom_point()+apatheme+
  facet_wrap("Fish_sp.x")

#ggtitle("Nutrients per year")+

# Evaluate dynamic of nutrients across years

ggplot(myxo_count, aes(x= mean_nitrogen,
                       y=psite_count))+
  geom_point()+apatheme+
  facet_wrap("Fish_sp.x")

# Evaluate dynamic of nutrients across years

ggplot(myxo_count, aes(x= mean_streamflow,
                       y=psite_count))+
  geom_point()+apatheme+
  facet_wrap("Fish_sp.x")


## GLM for all fish

# Both CatalogNumber and season as as random effects
library(splines)

m1 <- glmmTMB(psite_count ~ scaledYear+
                scaledmean_streamflow*scaledmean_nitrogen+
                scaledmean_temperature*scaledmean_nitrogen+
                scaledmean_temperature*scaledmean_streamflow+
                Dim.1*scaledmean_streamflow+
                Dim.2*scaledmean_streamflow+
                Dim.1*scaledmean_temperature+
                Dim.2*scaledmean_temperature+
                Dim.1*scaledmean_nitrogen+
                Dim.2*scaledmean_nitrogen+
                (1|Fish_sp.x/IndividualFishID)+
                (1|season),
              data = myxo_count,
              family = nbinom1)

m2 <- glmmTMB(psite_count ~ scaledYear+
                scaledmean_streamflow*scaledmean_nitrogen+
                scaledmean_temperature*scaledmean_nitrogen+
                scaledmean_temperature*scaledmean_streamflow+
                Dim.1*scaledmean_streamflow+
                Dim.2*scaledmean_streamflow+
                Dim.1*scaledmean_temperature+
                Dim.2*scaledmean_temperature+
                Dim.1*scaledmean_nitrogen+
                Dim.2*scaledmean_nitrogen+
                offset(scaled_TL_mm)+
                (1|Fish_sp.x/IndividualFishID)+
                (1|season),
              data = myxo_count,
              family = nbinom2)


m2 <- glmmTMB(psite_count ~ scaledYear+
                scaledmean_streamflow*scaledmean_nitrogen+
                scaledmean_temperature*scaledmean_nitrogen+
                scaledmean_temperature*scaledmean_streamflow+
                offset(scaled_TL_mm)+
                (1|Fish_sp.x/IndividualFishID)+
                (1|season),
              data = myxo_count,
              family = nbinom2)



AIC(m1,m2)         
summary(m2)

#Evaluate residuals
s=simulateResiduals(fittedModel=m2,n=250)
s$scaledResiduals
plot(s)

tab_model(m2)
plot_model(m2)+apatheme+geom_hline(yintercept=1, linetype="dashed", color = "black", size=0.5)

# Flow with CI and psite_genus

mydf <- ggpredict(m2, terms= c("scaledmean_nitrogen[all]")) 

plot(mydf,rawdata=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme

# Flow with CI and psite_genus

mydf <- ggpredict(m2, c("scaledmean_temperature[all]","Fish_sp.x")) 

plot(mydf,rawdata=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  apatheme

### Impact of multiple stressors on the abundance of myxobolus----

## GLM for CARVEL

# Transform total length

carvel_count$logTL_mm <- log(carvel_count$TotalLength_mm)

carvel_count$scaledmeanFlow <- scale(carvel_count$meanFlow )

carvel_count$scaledmeanN <- scale(carvel_count$meanN)

carvel_count$scaledmeanTemp <- scale(carvel_count$meanTemp)

carvel_count$scaledYear <- scale(carvel_count$YearCollected)


# Both CatalogNumber and season as as random effects
glm_count2 <- glmmTMB(psite_count ~ 
                        scaledmeanTemp*scaledmeanFlow*scaledmeanN+
                        offset(logTL_mm)+
                        (1|IndividualFishID)+ 
                        (1|season),
                      data = carvel_count,family=nbinom1()) # This model is working good, good dianostic and fit.

glm_count2 <- glmmTMB(psite_count ~ scaledYear+
                        scaledmeanTemp*scaledmeanFlow*scaledmeanN+
                        offset(logTL_mm)+
                        (1|IndividualFishID)+ 
                        (1|season),
                      data = carvel_count,family=nbinom2()) # rank deficient

summary(glm_count2)

#Evaluate residuals
s=simulateResiduals(fittedModel=glm_count2,n=250)
s$scaledResiduals
plot(s)

tab_model(glm_count2)

# Flow with CI and psite_genus

mydf <- ggpredict(glm_count2, c("scaledmeanN[all]")) 

plot(mydf,rawdata=TRUE,jitter=0.05,color=c("#5aae61","#762a83"))+
  labs(x = 'Stream flow (m3/sec)', y = 'Abundance of myxozoans',title=NULL)+
  apatheme

