### Jolee Code REU Project ###
### created 7/11/24        ###
### by Jolee Thirtyacre    ###
### and Connor Whalen      ###

# load packages
install.packages("tidyverse")
install.packages("janitor")
install.packages("ggeffects")
library(ggeffects)
library(tidyverse)
library(readr)
library(dbplyr)
library(ggarrange)
library(lubridate)
library(ggplot2)


plot_1 <- ggplot(df,
                 aes(x = , y = )) + geom_boxplot() + xlab("S")  + 
  ylab("")  + theme_bw()



# Jolee's code from another file
setwd("C:/Users/thirt/Documents/TUBRI_JT/data/processed")

P.vigilax.data <- read.csv("Pimephales_vigilax_processed_machine_readable_UPDATED_2024.07.10.csv")
I.punctatus.data <- read.csv("Ictalurus_punctatus_processed_machine_readable_UPDATED_2024.07.10.csv")


# pull out weights of each species
weight.vig <- P.vigilax.data$Weight_mg
weight.pun <- I.punctatus.data$Weight_mg

# pull out total lengths of each species
tot.length.vig <- P.vigilax.data$TotalLength_mm
tot.length.pun <- I.punctatus.data$TotalLength_mm

# pull out standard lengths of each species
stan.length.vig <- P.vigilax.data$StandardLength_mm
stan.length.pun <- I.punctatus.data$StandardLength_mm

# calculating Folton's condition factor (Mozsar et al. 2014)
    #P. vigilax
body.index.vig <- (weight.vig / (stan.length.vig ^ 3)) * 100

plot(body.index.vig)
boxplot(body.index.vig)

   #I. punctatus
body.index.pun <- (weight.pun / (stan.length.pun ^ 3)) * 100

plot(body.index.pun)
boxplot(body.index.pun)


### plots to look at associations to see if there's anything (preliminary) ###

# plot condition factor by time
plot(x = P.vigilax.data$YearCollected, 
     y = body.index.vig)
plot(x = I.punctatus.data$YearCollected, 
     y = body.index.pun)

# plot body size by time
    # standard length for body size
plot(x = P.vigilax.data$YearCollected, 
     y = P.vigilax.data$StandardLength_mm) #this is how we designed the study though - will show nothing

plot(x = I.punctatus.data$YearCollected, 
     y = I.punctatus.data$StandardLength_mm) #this is how we designed the study though - will show nothing

    # weight for body size
plot(x = P.vigilax.data$YearCollected, 
     y = P.vigilax.data$Weight_mg)

plot(x = I.punctatus.data$YearCollected, 
     y = I.punctatus.data$Weight_mg)

# plot total parasite burden by time
plot(x = P.vigilax.data$YearCollected, 
     y = P.vigilax.data$psite_count)
plot(x = I.punctatus.data$YearCollected, 
     y = I.punctatus.data$psite_count)

# plot condition factor vs parasite burden in total - may be interesting
plot(x = body.index.vig,
     y = P.vigilax.data$psite_count)
plot(x = body.index.pun,
     y = I.punctatus.data$psite_count)

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


    # plot condition factor vs binary code for life stage (aka "stag-specific parasite stuff")
boxplot(Pvig.life.stage$BinaryCode, body.index.vig)


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
                                                              "CYST.UNKN", # correct code? What's in the slides
                                                              "WORM.UNK.B",
                                                              "TREM.2", # exclude because it's a once-off?
                                                              "LEECH.KUT"
                                                              )))


# adding a binary code column that codes "Adult" as 1 and "Larval" as 0
Ipun.life.stage <- Ipun.life.stage %>%
  mutate(BinaryCode = ifelse(LifeStage == "Adult", 1, 0))


# plot condition factor vs binary code for life stage (aka "stag-specific parasite stuff")
boxplot(Ipun.life.stage$BinaryCode, body.index.pun)
