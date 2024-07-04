# How to use R
# A workshop
# Written by Chelsea Wood (chelwood@uw.edu)

getwd()

setwd("/Users/chelsealwood/Dropbox/Vault/University of Washington/Projects/CAREER/Data/TUBRI")

pim_vig_data<-read.csv("data/processed/Ictalurus_punctatus_processed_human_readable.csv")

View(pim_vig_data)

pim_vig_data

pim_vig_data$YearCollected

pim_vig_data$MONO.IP

pim_vig_data[,]

sum(pim_vig_data$MONO.IP)

mean(pim_vig_data$MONO.IP)

hist(pim_vig_data$MONO.IP)

plot(pim_vig_data$MONO.IP~pim_vig_data$YearCollected)

summary(lm(pim_vig_data$MONO.IP~pim_vig_data$YearCollected))
