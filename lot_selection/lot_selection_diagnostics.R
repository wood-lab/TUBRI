# Monitoring lot selection during dissections at TUBRI
# Chelsea Wood (chelwood@uw.edu)
# June-August 2024

# The purpose of this script is to allow dissectors to quickly assess levels of replication in each category
# (i.e., before-impact, after-impact, before-control, after-control) during dissections, so that we can 
# adjust lot selection on the fly if need be. We also want to make sure that fish body size range is constant
# across decades, that latitudes sampled before 1973 are the same latitudes sampled after, and that we are
# on track to hit our dissection goals.


# Load required packages
library(googledrive)
library(skimr)
library(lubridate)
library(tidyverse)


# First, bring in all of the meta-data you need to unite with your dissection data

meta_data<-read.csv("lot_selection/all_lots_okay_to_dissect.csv")
meta_data$combo<-paste(meta_data$CI,meta_data$decade,sep="_")



### PIMVIG

# Download latest file # GET THIS WORKING LATER
drive_download(as_id("16taxgLRu1-d_t8lT9mHVWOtxZ4MfQPvU7UBWHvLRfaI"), 
               path = "data/raw/pim_vig.csv", overwrite = TRUE)


# You can also do it the old-fashioned way

pim_vig_today<-read.csv("data/raw/Pimephales_vigilax_Datasheet_2024.06.25.csv")
length(pim_vig_today$CatalogNumber)


# Now merge dissection data with meta-data (all.x makes sure that you keep all individual fish from the same lot, even though
# they have the same catalog number).

pim_vig_with_metadata<-merge(pim_vig_today, meta_data, by.x = "CatalogNumber", by.y = "CatalogNumber", all.x = TRUE)
length(pim_vig_with_metadata$CatalogNumber)


# Plot to see where the dissected fish fall

plot(jitter(pim_vig_with_metadata$Latitude.y,30)~jitter(pim_vig_with_metadata$YearCollected.y,5))+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)


# Desired level of replication

pim_vig_desired<-read.csv(file="lot_selection/desired_replication/pim_vig_goal.csv")
pim_vig_desired<-pim_vig_desired[,-1]
colnames(pim_vig_desired)[2]<-"n_lots_desired"
pim_vig_desired$individual_fish_desired<-pim_vig_desired$n_lots_desired*4


# Actual level of replication

pim_vig_actual<-pim_vig_with_metadata %>%
  group_by(combo) %>%
  summarize(actual = n())


# Merge desired and actual

pim_vig_remaining<-merge(pim_vig_desired, pim_vig_actual, by.x = "combo", by.y = "combo", all.x = TRUE, all.y = TRUE)


# Calculate how many you have left

pim_vig_remaining$to_go<-(pim_vig_remaining$individual_fish_desired-ifelse(is.na(pim_vig_remaining$actual),0,pim_vig_remaining$actual))
pim_vig_remaining<-na.omit(pim_vig_remaining)
sum(pim_vig_remaining$actual,na.rm=T)
sum(pim_vig_remaining$to_go,na.rm=T)


# Check that body size is invariant through time

summary(glm(pim_vig_with_metadata$TotalLength_mm~pim_vig_with_metadata$YearCollected.x))
plot(pim_vig_with_metadata$TotalLength_mm~pim_vig_with_metadata$YearCollected.x)
summary(glm(pim_vig_with_metadata$TotalLength_mm~pim_vig_with_metadata$Latitude.y))
plot(pim_vig_with_metadata$TotalLength_mm~pim_vig_with_metadata$Latitude.y)


# Create a matrix that shows how many fish you are finishing per day

pim_vig_dailies<-pim_vig_with_metadata %>%
  group_by(DissectionDate) %>%
  summarize(actual = n())



### ICTPUN

# Download latest file # GET THIS WORKING LATER
drive_download(as_id("16taxgLRu1-d_t8lT9mHVWOtxZ4MfQPvU7UBWHvLRfaI"), 
               path = "data/raw/pim_vig.csv", overwrite = TRUE)


# You can also do it the old-fashioned way

ict_pun_today<-read.csv("data/raw/Ictalurus_punctatus_Datasheet_2024.06.28.csv")
length(ict_pun_today$CatalogNumber)


# Now merge dissection data with meta-data (all.x makes sure that you keep all individual fish from the same lot, even though
# they have the same catalog number).

ict_pun_with_metadata<-merge(ict_pun_today, meta_data, by.x = "CatalogNumber", by.y = "CatalogNumber", all.x = TRUE)
length(ict_pun_with_metadata$CatalogNumber)


# Plot to see where the dissected fish fall

plot(jitter(ict_pun_with_metadata$Latitude.y,30)~jitter(ict_pun_with_metadata$YearCollected.y,5))+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)


# Desired level of replication

ict_pun_desired<-read.csv(file="lot_selection/desired_replication/ict_pun_goal.csv")
ict_pun_desired<-ict_pun_desired[,-1]
colnames(ict_pun_desired)[2]<-"n_lots_desired"
ict_pun_desired$individual_fish_desired<-ict_pun_desired$n_lots_desired*2
sum(ict_pun_desired$individual_fish_desired)


# Actual level of replication

ict_pun_actual<-ict_pun_with_metadata %>%
  group_by(combo) %>%
  summarize(actual = n())


# Merge desired and actual

ict_pun_remaining<-merge(ict_pun_desired, ict_pun_actual, by.x = "combo", by.y = "combo", all.x = TRUE, all.y = TRUE)


# Calculate how many you have left

ict_pun_remaining$to_go<-(ict_pun_remaining$individual_fish_desired-ifelse(is.na(ict_pun_remaining$actual),0,ict_pun_remaining$actual))
ict_pun_remaining<-na.omit(ict_pun_remaining)
sum(ict_pun_remaining$actual,na.rm=T)
sum(ict_pun_remaining$to_go,na.rm=T)


# Check that body size is invariant through time

summary(glm(ict_pun_with_metadata$TotalLength_mm~ict_pun_with_metadata$YearCollected.x))
plot(ict_pun_with_metadata$TotalLength_mm~ict_pun_with_metadata$YearCollected.x)
summary(glm(ict_pun_with_metadata$TotalLength_mm~ict_pun_with_metadata$Latitude.y))
plot(ict_pun_with_metadata$TotalLength_mm~ict_pun_with_metadata$Latitude.y)


# Create a matrix that shows how many fish you are finishing per day

ict_pun_dailies<-ict_pun_with_metadata %>%
  group_by(DissectionDate) %>%
  summarize(actual = n())



### NOTATH

# You can also do it the old-fashioned way

not_ath_today<-read.csv("data/raw/Notropis_atherinoides_Datasheet_2024.01.04.csv")
length(not_ath_today$CatalogNumber)


# Now merge dissection data with meta-data (all.x makes sure that you keep all individual fish from the same lot, even though
# they have the same catalog number).

not_ath_with_metadata<-merge(not_ath_today, meta_data, by.x = "CatalogNumber", by.y = "CatalogNumber", all.x = TRUE)
length(not_ath_with_metadata$CatalogNumber)


# Plot to see where the dissected fish fall

plot(jitter(not_ath_with_metadata$Latitude.y,30)~jitter(not_ath_with_metadata$YearCollected.y,5))+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)


# Desired level of replication

not_ath_desired<-read.csv(file="lot_selection/desired_replication/not_ath_goal.csv")
not_ath_desired<-not_ath_desired[,-1]
colnames(not_ath_desired)[2]<-"n_lots_desired"
not_ath_desired$individual_fish_desired<-not_ath_desired$n_lots_desired*4
sum(not_ath_desired$individual_fish_desired)


# Actual level of replication

not_ath_actual<-not_ath_with_metadata %>%
  group_by(combo) %>%
  summarize(actual = n())


# Merge desired and actual

not_ath_remaining<-merge(not_ath_desired, not_ath_actual, by.x = "combo", by.y = "combo", all.x = TRUE, all.y = TRUE)


# Calculate how many you have left

not_ath_remaining$to_go<-(not_ath_remaining$individual_fish_desired-ifelse(is.na(not_ath_remaining$actual),0,not_ath_remaining$actual))
not_ath_remaining<-na.omit(not_ath_remaining)
sum(not_ath_remaining$actual,na.rm=T)
sum(not_ath_remaining$to_go,na.rm=T)


# Check that body size is invariant through time

summary(glm(not_ath_with_metadata$TotalLength_mm~not_ath_with_metadata$YearCollected.x))
plot(not_ath_with_metadata$TotalLength_mm~not_ath_with_metadata$YearCollected.x)
summary(glm(not_ath_with_metadata$TotalLength_mm~not_ath_with_metadata$Latitude.y))
plot(not_ath_with_metadata$TotalLength_mm~not_ath_with_metadata$Latitude.y)


# Create a matrix that shows how many fish you are finishing per day

not_ath_dailies<-not_ath_with_metadata %>%
  group_by(DissectionDate) %>%
  summarize(actual = n())
View(not_ath_dailies)


# Calculate how many fish you need to do per day to hit your target

working_days_remaining<-41
read.csv("lot_selection/desired_replication/all_spp_goal.csv")

200+144+228+100+224+136+65

1097/41
7*4
26/4


# Replacing lots that we can't find in the collection
# Look at the first three digits of the lots you are missing.
# Find lots with those digits (or close) on the "in bounds" shelves in the collection.
# Of all the close lots, check that they are in the category you're looking for and, if you have more than one,
# randomly select the ones you'll use.

