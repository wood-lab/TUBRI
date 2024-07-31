### Lot selection - samples to be sent to Germany
### Dakeishla M. Diaz Morales
### Updated July 30, 2024

### This script is meant to randomly select samples to be sent to Germany. Initially we agreed with Bernd and Milen to select 20 samples per muscle/liver of fish per control/impact category and before/after the clean water act. This will result in 560 samples of liver and 560 samples of muscle.



## Lot selection for PIMVIG----

#Upload dataset
library(readr)
library(dplyr)

PIMVIG <- read_csv("data/processed/Pimephales_vigilax_processed_machine_readable_UPDATED_2024.07.10.csv")
View(PIMVIG)
str(PIMVIG)

library(dplyr)

# Summarize data by jar/lot
PIMVIGb <- PIMVIG %>%
  group_by(CatalogNumber,combo)%>%
summarize(CI = any(CI),.groups = 'drop')%>% #This does not work well as of now but gets the job done: having a list of jars from which 4 fish were dissected
  as.data.frame()

#View(PIMVIGb)

# Count how many jars were selected per combo
PIMVIGb %>% count(combo)


#Select random lots per combo
LotPV_control_1954to1963<-PIMVIGb[ sample( PIMVIGb$CatalogNumber, which( PIMVIGb$combo == "control_1954-1963"), 1, replace = F), ]
LotPV_control_1964to1973<-PIMVIGb[ sample( which( PIMVIGb$combo == "control_1964-1973"), 10, replace = F), ]

LotPV_control_1974to1983<-PIMVIGb[ sample( which( PIMVIGb$combo == "control_1974-1983"), 3, replace = F), ]
LotPV_control_1984to1993<-PIMVIGb[ sample( which( PIMVIGb$combo == "control_1984-1993"), 1, replace = F), ]#not enough for 5
LotPV_control_1994to2003<-PIMVIGb[ sample( which( PIMVIGb$combo == "control_1994-2003"), 2, replace = F), ]#not enough for 5
LotPV_control_2004to2013<-PIMVIGb[ sample( which( PIMVIGb$combo == "control_2004-2013"), 2, replace = F), ]#not enough for 5

LotPV_impact_1954to1963<-PIMVIGb[ sample( which( PIMVIGb$combo == "impact_1954-1963"), 5, replace = F), ]
LotPV_impact_1964to1973<-PIMVIGb[ sample( which( PIMVIGb$combo == "impact_1964-1973"), 5, replace = F), ]
LotPV_impact_1974to1983<-PIMVIGb[ sample( which( PIMVIGb$combo == "impact_1974-1983"), 5, replace = F), ]
LotPV_impact_1984to1993<-PIMVIGb[ sample( which( PIMVIGb$combo == "impact_1984-1993"), 5, replace = F), ]
LotPV_impact_1994to2003<-PIMVIGb[ sample( which( PIMVIGb$combo == "impact_1994-2003"), 5, replace = F), ]
LotPV_impact_2004to2013<-PIMVIGb[ sample( which( PIMVIGb$combo == "impact_2004-2013"), 4, replace = F), ]#don't bother because we don't have this decade for control


# Put data together
PIMVIG_postDE <-rbind.data.frame(LotPV_control_1954to1963,
                                 LotPV_control_1964to1973,
                                 LotPV_control_1974to1983,
                                 LotPV_control_1984to1993,
                                 LotPV_control_1994to2003,
                                 LotPV_control_2004to2013,
                                 LotPV_impact_1954to1963,
                                 LotPV_impact_1964to1973,
                                 LotPV_impact_1974to1983,
                                 LotPV_impact_1984to1993,
                                 LotPV_impact_1994to2003,
                                 LotPV_impact_2004to2013)

# Export table 
write.csv(per_vig_matrix,file="lot_selection/PIMVIG_postDE.csv")



