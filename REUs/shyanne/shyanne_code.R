# Shyanne's data exploration, plotting, and analysis
# by Shyanne Christner and Chelsea Wood
# christner.shyanne03@gmail.com and chelwood@uw.edu
# Started 4 July 2024

# R is really useful for doing quick tallies and plots.
# Let's do some together!

library(tidyverse)

ict_pun_data<-read.csv("data/processed/Ictalurus_punctatus_processed_human_readable.csv")


# Let's do some quick data tallies

a<-"shyanne"

tally_up <- ict_pun_data %>%
  group_by(combo) %>%
  summarize(mean.MONO.IP = mean(MONO.IP))

tally_up <- ict_pun_data %>%
  group_by(YearCollected) %>%
  summarize(mean.MONO.IP = mean(MONO.IP))
View(tally_up)


# Let's do some quick plots

plot(ict_pun_data$MONO.IP~ict_pun_data$Latitude)+abline(v=30.76)


# How about some quick stats?

summary(lm(pim_vig_data$MONO.IP~pim_vig_data$Latitude))


