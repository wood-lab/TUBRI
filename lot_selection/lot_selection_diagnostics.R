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

pim_vig_today<-read.csv("data/Pimephales_vigilax_Datasheet_2024.06.13.csv")
length(pim_vig_today$CatalogNumber)


# Now merge dissection data with meta-data (all.x makes sure that you keep all individual fish from the same lot, even though
# they have the same catalog number).

pim_vig_with_metadata<-merge(pim_vig_today, meta_data, by.x = "CatalogNumber", by.y = "CatalogNumber", all.x = TRUE)
length(pim_vig_with_metadata$CatalogNumber)


# Plot to see where the dissected fish fall

plot(jitter(pim_vig_with_metadata$Latitude.y,1)~pim_vig_with_metadata$YearCollected.y)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)


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
sum(pim_vig_remaining$actual,na.rm=T)
sum(pim_vig_remaining$to_go,na.rm=T)


# Check that body size is invariant through time

summary(glm(pim_vig_with_metadata$TotalLength_mm~pim_vig_with_metadata$YearCollected.x))
plot(pim_vig_with_metadata$TotalLength_mm~pim_vig_with_metadata$YearCollected.x)
summary(glm(pim_vig_with_metadata$TotalLength_mm~pim_vig_with_metadata$Latitude.y))
plot(pim_vig_with_metadata$TotalLength_mm~pim_vig_with_metadata$Latitude.y)


# Calculate how many fish you need to do per day to hit your target

working_days_remaining<-45
read.csv("lot_selection/desired_replication/all_spp_goal.csv")

200+144+228+100+224+136+153

1185/45




