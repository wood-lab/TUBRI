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
PIMVIG %>% count(combo)

PIMVIGb <- PIMVIG %>%
  group_by(CatalogNumber,combo)%>%
summarize(CI = any(CI),.groups = 'drop')%>%
  as.data.frame()

View(PIMVIGb)

LotPV_control_1954to1963<-PIMVIGb[ sample( PIMVIGb$CatalogNumber, which( PIMVIGb$combo == "control_1954-1963"), 1, replace = F), ]
LotPV_control_1964to1973<-PIMVIG[ sample( which( PIMVIG$combo == "control_1964-1973"), 10, replace = F), ]

control_1974to1983<-PIMVIG[ sample( which( PIMVIG$combo == "control_1974-1983"), 3, replace = F), ]
control_1984to1993<-PIMVIG[ sample( which( PIMVIG$combo == "control_1984-1993"), 1, replace = F), ]#not enough for 5
control_1994to2003<-PIMVIG[ sample( which( PIMVIG$combo == "control_1994-2003"), 2, replace = F), ]#not enough for 5
control_2004to2013<-PIMVIG[ sample( which( PIMVIG$combo == "control_2004-2013"), 2, replace = F), ]#not enough for 5


impact_1954to1963<-PIMVIG[ sample( which( PIMVIG$combo == "impact_1954-1963"), 5, replace = F), ]
impact_1964to1973<-PIMVIG[ sample( which( PIMVIG$combo == "impact_1964-1973"), 5, replace = F), ]
impact_1974to1983<-PIMVIG[ sample( which( PIMVIG$combo == "impact_1974-1983"), 5, replace = F), ]
impact_1984to1993<-PIMVIG[ sample( which( PIMVIG$combo == "impact_1984-1993"), 5, replace = F), ]
impact_1994to2003<-PIMVIG[ sample( which( PIMVIG$combo == "impact_1994-2003"), 5, replace = F), ]
impact_2004to2013<-PIMVIG[ sample( which( PIMVIG$combo == "impact_2004-2013"), 4, replace = F), ]#don't bother because we don't have this decade for control




