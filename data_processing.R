# Processing data from raw to analyzable form
# Chelsea Wood (chelwood@uw.edu)
# June-August 2024

# The purpose of this script is to clean the raw dataset, strip out unneeded fields, sum parasite counts across 
# organs, and generally get the thing into analyzable form. If you ever need things from the raw dataset, they 
# will still be there for you to access!  But the processed dataset will have what is needed for most analyses.


# Load required packages
library(googledrive)
library(skimr)
library(lubridate)
library(tidyverse)


# Here's the workflow that you need to follow for each fish species:
# 1. List all parasite species detected and figure out which are valid (i.e., lumping some, 
#    discarding those that aren't parasites).
# 2. Sum across organs within valid parasite species and double counts for those organs that need to be doubled.
# 3. Create a dataset containing only the data you want (get rid of organ and voucher headers).
# 4. Melt the dataset so that parasite counts are in one column and parasite names are in another column.

# Additional processing:
# 5. Trim out any parasite species that occur at <5% prevalence WITHIN an archipelago.
# 6. Add fish-level meta-data (from Stuart/BZ/Bev), including date, time, coordinates, depth, habitat type, current, and conditions.
# 7. Add island-level meta-data (from various sources), including chl-a, island size.
# 8. DON'T FORGET that there are Palmyra and Kiritimati data too (from Dani and Stuart). Include them?


##### Check with KATIE - Which organs need to be doubled?


# First, bring in all of the meta-data you need to unite with your dissection data

meta_data<-read.csv("lot_selection/all_lots_okay_to_dissect.csv")
meta_data$combo<-paste(meta_data$CI,meta_data$decade,sep="_")



### PIMVIG

# Download latest file # GET THIS WORKING LATER
drive_download(as_id("16taxgLRu1-d_t8lT9mHVWOtxZ4MfQPvU7UBWHvLRfaI"), 
               path = "data/raw/pim_vig.csv", overwrite = TRUE)


# You can also do it the old-fashioned way

pim_vig_today<-read.csv("data/Pimephales_vigilax_Datasheet_2024.06.18.csv")
length(pim_vig_today$CatalogNumber)


# Now merge dissection data with meta-data (all.x makes sure that you keep all individual fish from the same lot, even though
# they have the same catalog number).

pim_vig_with_metadata<-merge(pim_vig_today, meta_data, by.x = "CatalogNumber", by.y = "CatalogNumber", all.x = TRUE)
length(pim_vig_with_metadata$CatalogNumber)


# Plot to see where the dissected fish fall

plot(jitter(pim_vig_with_metadata$Latitude.y,1)~pim_vig_with_metadata$YearCollected.y)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

str(pim_vig_with_metadata)

pim_vig_headers<-ls(pim_vig_with_metadata)
pim_vig_headers<-as.factor(pim_vig_headers)

# Add parasites across organs

CEST.GODZ<-cest.godz.intestine
CEST.L<-CEST.L.FLUSH
CEST.PEA<-CEST.PEA.Intestine 
CEST.UNK
CEST.UNK.Intestine              
CEST.VIT.INTESTINE
COP.GRAB.GILL
CYST.EMPTY.ConnectiveTissue     
CYST.EMPTY.Stomach               
CYST.LB.EYE 
CYST.LB.GILL
CYST.UNK.AnalFin                
CYST.UNK.CaudalFin
CYST.UNK.CONNECTIVETISSUES 
CYST.UNK.Eye 
CYST.UNK.PectoralFin             
cyst.unkn.gill
META.HS.CaudalFin                
META.UNK.AnalFin                 
META.UNK.CAUDALFIN              
META.UNK.DorsalFin               
META.UNK.EYE                     
META.UNK.GILL                   
META.UNK.LIVER                   
meta.unkn.connectivetissue      
meta.unkn.intestine              
MONO.DACT.GILL                   
MONO.GYRO.Flush                  
MONO.GYRO.Gill
MYX.GO.CONNECTIVETISSUE         
MYX.GO.GONAD                     
MYX.GO.VOUCH                     
MYXO.BC.GILL                    
MYXO.BC.VOUCH                    
MYXO.GO.Eye                      
MYXO.GO.Kidney                  
MYXO.SBAD.GILL                   
MYXO.THEL.CaudalFin             
MYXO.THEL.ConnectiveTissue       
MYXO.THEL.dorsalFIN              
MYXO.THEL.GILL                  
myxo.thel.Skin                   
MYXO.UNK.Connective.Tissue      
MYXO.UNK.GILL                    
NEM.CONA.CONNECTIVETISSUES       
NEM.CONA.Intestine              
NEM.CYST.CONNECTIVETISSUE        
NEM.CYST.VOUCH                  
NEM.LARV.Flush                   
NEM.LARV.VOUCH                   
NEM.TRBK.intestine              
trem.buc.AnalFin                
trem.buc.CaudalFin
TREM.BUC.PectoralFin             
TREM.MET.AF.PectoralFin          
TREM.MET.BodyCavity              
TREM.MET.CaudalFin              
TREM.MET.ConnectiveTissue        
TREM.MET.FLUSH                   
TREM.MET.GILL                   
TREM.MET.GS.Intestine            
TREM.MET.Intestine               
TREM.MET.Kidney                 
TREM.MET.LIVER                   
TREM.META.UCC.Flush             
TREM.METAF.AnalFin               
trem.metaf.caudalfin            
TREM.METAF.GILL                  
TREM.METAF.PelvicFin             
TREM.METAGS.CaudalFin            
TREM.METAGS.CONNECTIVETISSUE     
TREM.METAGS.FLUSH               
TREM.METAGS.GILL                 
TREM.METAGS.Kidney               
TREM.METAGS.Liver               
TREM.METAGS.SKIN                 
TREM.METAGS.Stomach              
TREM.MS.FLUSH                    
WORM.A.Flush


per_vig_processed_data<-rbind(pim_vig_with_metadata$CatalogNumber,pim_vig_with_metadata$YearCollected.x,
                              pim_vig_with_metadata$MonthCollected.x,pim_vig_with_metadata$DayCollected.x,
                              pim_vig_with_metadata$IndividualFishID,pim_vig_with_metadata$Dissector_and_Examiner,
                              pim_vig_with_metadata$DissectionDate,pim_vig_with_metadata$Sex,
                              pim_vig_with_metadata$TotalLength_mm,pim_vig_with_metadata$StandardLength_mm,
                              pim_vig_with_metadata$Weight_mg,pim_vig_with_metadata$combo,
                              pim_vig_with_metadata$decade.x,Latitude.y,Longitude.y,)
#name columns

