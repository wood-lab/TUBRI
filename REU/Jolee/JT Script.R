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