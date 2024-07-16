### Jolee Code REU Project ###
### created 7/11/24        ###
### by Jolee Thirtyacre    ###
### and Connor Whalen      ###

# install packages
install.packages("tidyverse")
install.packages("janitor")
install.packages("ggeffects")

# load packages
library(ggeffects)
library(tidyverse)
library(readr)
library(dbplyr)
library(ggarrange)
library(lubridate)
library(ggplot2)


# Connor's wd assignment and data assignments
P.vigilax.data <- read.csv("data/processed/Pimephales_vigilax_processed_machine_readable_UPDATED_2024.07.10.csv") %>% 
  janitor::clean_names()
I.punctatus.data <- read.csv("data/processed/Ictalurus_punctatus_processed_machine_readable_UPDATED_2024.07.10.csv")  %>% 
  janitor::clean_names()
N.atherinoides.data <- read.csv("data/processed/Notropis_atherinoides_processed_machine_readable.csv") %>% 
  janitor::clean_names()


# template for ggplot --> fill in the blank areas with data
plot_1 <- ggplot(df,
                 aes(x = , y = )) + geom_boxplot() + xlab("S")  + 
  ylab("")  + theme_bw()




# Jolee's code from another file
setwd("C:/Users/thirt/Documents/TUBRI_JT/data/processed")

### edited these to also use the tidyverse -- if you add the pipe followed by the clean names commmand, you can
### call objects directly without using df$column --> this isn't true fully for our datasets since they use the same 
### column names in each, but in areas where you are specifying a df (like in ggplot), using the tidyverse will 
### be more straightforwardm -- we can talk more about this on Friday! 


P.vigilax.data <- read.csv("Pimephales_vigilax_processed_machine_readable_UPDATED_2024.07.10.csv")
I.punctatus.data <- read.csv("Ictalurus_punctatus_processed_machine_readable_UPDATED_2024.07.10.csv")
N.atherinoides.data <- read.csv("Notropis_atherinoides_processed_machine_readable.csv")


# pull out weights of each species
weight.vig <- P.vigilax.data$Weight_mg
weight.pun <- I.punctatus.data$Weight_mg
weight.ath <- N.atherinoides.data$Weight_mg


# pull out total lengths of each species
tot.length.vig <- P.vigilax.data$TotalLength_mm
tot.length.pun <- I.punctatus.data$TotalLength_mm
tot.length.ath <- N.atherinoides.data$TotalLength_mm


# pull out standard lengths of each species
stan.length.vig <- P.vigilax.data$StandardLength_mm
stan.length.pun <- I.punctatus.data$StandardLength_mm
stan.length.ath <- N.atherinoides.data$StandardLength_mm

# calculating Folton's condition factor (Mozsar et al. 2014)
    #P. vigilax
body.index.vig <- (weight.vig / (stan.length.vig ^ 3)) * 100 

plot(body.index.vig)
boxplot(body.index.vig)



   #I. punctatus
body.index.pun <- (weight.pun / (stan.length.pun ^ 3)) * 100

plot(body.index.pun)
boxplot(body.index.pun)


    #N. atherinoides
body.index.ath <- (weight.ath / (stan.length.ath ^ 3)) * 100

plot(body.index.ath)
boxplot(body.index.ath)

### plots to look at associations to see if there's anything (preliminary) ###

# plot condition factor by time
plot(x = P.vigilax.data$YearCollected, 
     y = body.index.vig)
plot(x = I.punctatus.data$YearCollected, 
     y = body.index.pun)
plot(x = N.atherinoides.data$YearCollected, 
     y = body.index.ath)

# plot body size by time
    # standard length for body size
plot(x = P.vigilax.data$YearCollected, 
     y = P.vigilax.data$StandardLength_mm) #this is how we designed the study though - will show nothing

plot(x = I.punctatus.data$YearCollected, 
     y = I.punctatus.data$StandardLength_mm) #this is how we designed the study though - will show nothing

plot(x = N.atherinoides.data$YearCollected, 
     y = N.atherinoides.data$StandardLength_mm)

    # weight for body size
plot(x = P.vigilax.data$YearCollected, 
     y = P.vigilax.data$Weight_mg)

plot(x = I.punctatus.data$YearCollected, 
     y = I.punctatus.data$Weight_mg)

plot(x = N.atherinoides.data$YearCollected, 
     y = N.atherinoides.data$Weight_mg)

# plot total parasite burden by time
plot(x = P.vigilax.data$YearCollected, 
     y = P.vigilax.data$psite_count)

plot(x = I.punctatus.data$YearCollected, 
     y = I.punctatus.data$psite_count)

plot(x = N.atherinoides.data$YearCollected, 
     y = N.atherinoides.data$psite_count) # more variability than the others

# plot condition factor vs parasite burden in total - may be interesting
plot(x = body.index.vig,
     y = P.vigilax.data$psite_count)

plot(x = body.index.pun,
     y = I.punctatus.data$psite_count)

plot(x = body.index.ath,
     y = N.atherinoides.data$psite_count) # more variability than the others

# Giving adult and young binary IDs so we can look at life stage compared to parasite burden
library(dplyr)

## P.VIG - add a column to identify species as an adult or larval stage
Pvig.life.stage <- P.vigilax.data %>%
  mutate(LifeStage = ifelse(psite_spp %in% c("ACANTH.BIGB", 
                                               "ACANTH.BKR",
                                               "NEM.CONA",
                                               "NEM.DICH",
                                               "NEM.LARV",
                                               "NEM.TRBK",
                                               "NEM.UNK",
                                               "CEST.GODZ",
                                               "CEST.NTG",
                                               "CEST.L",
                                               "CEST.VIT",
                                               "CEST.UNK",
                                               "MONO.GYRO",
                                               "MONO.DACT",
                                               "MONO.LG",
                                               "TREM.RNDA",
                                               "TREM.UNK"),
                            "Adult", "Larval"))

    # removing the species I'm not using
Pvig.life.stage <- subset(Pvig.life.stage, !(psite_spp %in% c("COPE.CL",
                                                      "COPE.GRAB",
                                                      "CYST.EMPTY",
                                                      "CYST.NEH",
                                                      "CYST.NP",
                                                      "WORM.A")))


    # adding a binary code column that codes "Adult" as 1 and "Larval" as 0
Pvig.life.stage <- Pvig.life.stage %>%
  mutate(BinaryCode = ifelse(LifeStage == "Adult", 1, 0))


    # plot condition factor vs binary code for life stage (aka "stage-specific parasite stuff")
boxplot(Pvig.life.stage$BinaryCode, body.index.vig)


# QUESTION: Is there a way to control for the variability of sample sizes for these data (aka. there are 
# way more larval than adult worms). Is this important?

## I.PUN - adding columns for life stage and binary code

Ipun.life.stage <- I.punctatus.data %>%
  mutate(LifeStage = ifelse(psite_spp %in% c("NEM.CYST",
                                             "TREM.LG",
                                             "TREM.META", # is this the same as META.GG - in data as TREM.META.GG
                                             "META.GG", # is this the same as TREM.META - in data as TREM.META.GG
                                             "META.UNK", # is this TREM.META.UNK in the data?
                                             "MYX.TAIL",
                                             "MYX.GEOM",
                                             "CYST.MX",
                                             "CYST.BD",
                                             "CYST.BLACK",
                                             "TREM.META.UNK"),
                            "Larval", "Adult")) # switching because there are more adults than larvae - less code to write

    # removing the species I'm not using
Ipun.life.stage <- subset(Ipun.life.stage, !(psite_spp %in% c("NMORPH",
                                                              "TREM.LARV",
                                                              "TREM.MID",
                                                              "CYST.UNKN", # correct code? This is what is in the slides
                                                              "WORM.UNK.B",
                                                              "TREM.2", # exclude because it's a once-off?
                                                              "LEECH.KUT"
                                                              )))


    # adding a binary code column that codes "Adult" as 1 and "Larval" as 0
Ipun.life.stage <- Ipun.life.stage %>%
  mutate(BinaryCode = ifelse(LifeStage == "Adult", 1, 0))


    # plot condition factor vs binary code for life stage (aka "stag-specific parasite stuff")
boxplot(Ipun.life.stage$BinaryCode, body.index.pun)


## N.ATH - adding columns for life stage and binary code

Nath.life.stage <- N.atherinoides.data %>%
  mutate(LifeStage = ifelse(psite_spp %in% c("NEM.CYST", # go through and check to see if dataset has any other not in the slides
                                             "TREM.LARV",
                                             "TREM.META",
                                             "MYX.SP", # MYX.SP.GILL?
                                             "MYX.DIAT",
                                             "MYX.ROUND",
                                             "CYST.UNA", # include?
                                             "CYST.UNK", # include?
                                             "TREM.CYST"), 
                            
                            "Larval", "Adult")) # switching because there are more adults than larvae - less code to write

# removing the species I'm not using
Nath.life.stage <- subset(Nath.life.stage, !(psite_spp %in% c("CYST.NEH",
                                                              "COPE.A",
                                                              "CYST.SAC"
)))


# adding a binary code column that codes "Adult" as 1 and "Larval" as 0
Nath.life.stage <- Nath.life.stage %>%
  mutate(BinaryCode = ifelse(LifeStage == "Adult", 1, 0))


# plot condition factor vs binary code for life stage (aka "stag-specific parasite stuff")
boxplot(Nath.life.stage$BinaryCode, body.index.ath)