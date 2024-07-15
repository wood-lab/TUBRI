#Script for Imani to select random monogeneans from Notropis atherinoides
#Revised July 12, 2024
#Authors: Imani Jones, Dakeishla Diaz, Chelsea Wood


NOTATH <- read_csv("data/processed/Notropis_atherinoides_processed_machine_readable.csv")


#Make a subset of only monogeneans

NOTATH_mon <- subset(NOTATH, psite_spp=="MONO.ALL"|psite_spp=="MONO.UNK") 

#Set a theme for all your plots
apatheme= theme_bw(base_size = 11,base_family = "sans")+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())

#Plot

library(ggplot2)
NOTATH_plot <-ggplot(NOTATH_mon,aes(YearCollected,psite_count,color = CI))+
  geom_point(size=4)+apatheme

NOTATH_plot

#Make a selection of monogeneans with at least 5 worms per fish. Then we can mount 5 worms per fish to get an idea of how the diversity looks per fish. 
NOTATH_monb <- subset(NOTATH_mon, psite_count>5) 

#Count how many vouchers per combo
library(dplyr)
NOTATH_mon %>% count(combo)

#Look at the distribution of voucher
library(ggplot2)
NOTATH_plot <-ggplot(NOTATH_monb,aes(YearCollected,psite_count,color = CI))+
  geom_point(size=4)+apatheme

NOTATH_plot

#Create data frame with final vouchers 

NOTATH_finalmon <- cbind.data.frame(NOTATH_monb$IndividualFishID,NOTATH_monb$YearCollected,
                                    NOTATH_monb$combo,NOTATH_monb$psite_count)


write.csv(NOTATH_finalmon, 
          file="REU/Imani J/NOTATH_finalmon.csv")




