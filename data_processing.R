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
library(reshape2)


# Here's the workflow that you need to follow for each fish species:
# 1. List all parasite species detected and figure out which are valid (i.e., lumping some, 
#    discarding those that aren't parasites).
# 2. Sum across organs within valid parasite species and double counts for those organs that need to be doubled.
# 3. Create a dataset containing only the data you want (get rid of organ and voucher headers).
# 4. Melt the dataset so that parasite counts are in one column and parasite names are in another column.


# First, bring in all of the meta-data you need to unite with your dissection data

meta_data<-read.csv("lot_selection/all_lots_okay_to_dissect.csv")
meta_data$combo<-paste(meta_data$CI,meta_data$decade,sep="_")



### PIMVIG

# Download latest file # GET THIS WORKING LATER
drive_download(as_id("16taxgLRu1-d_t8lT9mHVWOtxZ4MfQPvU7UBWHvLRfaI"), 
               path = "data/raw/pim_vig.csv", overwrite = TRUE)


# You can also do it the old-fashioned way

pim_vig_today<-read.csv("data/Pimephales_vigilax_Datasheet_2024.06.25.csv")
length(pim_vig_today$CatalogNumber)


# Now merge dissection data with meta-data (all.x makes sure that you keep all individual fish from the same lot, even though
# they have the same catalog number).

pim_vig_with_metadata<-merge(pim_vig_today, meta_data, by.x = "CatalogNumber", by.y = "CatalogNumber", all.x = TRUE)
length(pim_vig_with_metadata$CatalogNumber)


# Plot to see where the dissected fish fall

plot(jitter(pim_vig_with_metadata$Latitude.y,30)~jitter(pim_vig_with_metadata$YearCollected.y,5))+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

str(pim_vig_with_metadata)

pim_vig_headers<-ls(pim_vig_with_metadata)
pim_vig_headers<-as.factor(pim_vig_headers)

# Add parasites across organs and double what needs to be doubled.

ACANTH.BIGB<-pim_vig_with_metadata$ACANTH.BIGB.Intestine            
ACANTH.BKR<-pim_vig_with_metadata$ACANTH.BKR.Intestine            
CEST.COMP<-(as.numeric(pim_vig_with_metadata$CEST.COMP.FLUSH)+pim_vig_with_metadata$CEST.COMP.Intestine)              
CEST.GODZ<-pim_vig_with_metadata$cest.godz.intestine
CEST.L<-pim_vig_with_metadata$CEST.L.FLUSH                   
CEST.NTG<-pim_vig_with_metadata$CEST.NTG.INTESTINE         
CEST.PEA<-pim_vig_with_metadata$CEST.PEA.Intestine               
CEST.UNK<-(pim_vig_with_metadata$CEST.UNK.Intestine)           
CEST.VIT<-(pim_vig_with_metadata$CEST.VIT.INTESTINE)
COPE.CL<-2*(pim_vig_with_metadata$COPE.CL.GILL)               
COPE.GRAB<-2*(pim_vig_with_metadata$COPE.GRAB.GILL)                   
# not a parasite: CYST.EMPTY<-(pim_vig_with_metadata$CYST.EMPTY.ConnectiveTissue+as.numeric(pim_vig_with_metadata$CYST.EMPTY.Stomach))              
# not a parasite: CYST.LB<-2*(pim_vig_with_metadata$CYST.LB.EYE+pim_vig_with_metadata$CYST.LB.GILL)                  
# not a parasite: CYST.UNK<-(pim_vig_with_metadata$CYST.UNK.AnalFin+pim_vig_with_metadata$CYST.UNK.CaudalFin+
             pim_vig_with_metadata$CYST.UNK.CONNECTIVETISSUES+pim_vig_with_metadata$CYST.UNK.DorsalFin+
             (2*pim_vig_with_metadata$CYST.UNK.Eye)+pim_vig_with_metadata$CYST.UNK.Flush+
             pim_vig_with_metadata$CYST.UNK.intestines+as.numeric(pim_vig_with_metadata$CYST.UNK.Kidney)+
             pim_vig_with_metadata$CYST.UNK.LIVER+(2*pim_vig_with_metadata$CYST.UNK.PectoralFin)+
             as.numeric(pim_vig_with_metadata$CYST.UNK.PelvicFin)+(2*pim_vig_with_metadata$cyst.unkn.gill))                   
META.HS<-pim_vig_with_metadata$META.HS.CaudalFin                
META.UNK<-(pim_vig_with_metadata$META.UNK.AnalFin+pim_vig_with_metadata$META.UNK.CAUDALFIN+
             as.numeric(pim_vig_with_metadata$META.UNK.DorsalFin)+(2*pim_vig_with_metadata$META.UNK.EYE)+
             pim_vig_with_metadata$META.UNK.FLUSH+(2*pim_vig_with_metadata$META.UNK.GILL)+
             as.numeric(pim_vig_with_metadata$META.UNK.LIVER)+(2*pim_vig_with_metadata$META.UNK.PectoralFin)+
             (2*as.numeric(pim_vig_with_metadata$META.UNK.PelvicFin))+pim_vig_with_metadata$meta.unkn.connectivetissue+
             pim_vig_with_metadata$meta.unkn.intestine)
MONO.DACT<-(2*pim_vig_with_metadata$MONO.DACT.GILL)              
MONO.GYRO<-pim_vig_with_metadata$MONO.GYRO.Flush+(2*pim_vig_with_metadata$MONO.GYRO.Gill)              
MONO.LG<-2*(pim_vig_with_metadata$MONO.LG.Gill+pim_vig_with_metadata$MONO.LG.PectoralFin)           
MYX.GEOM<-pim_vig_with_metadata$MYX.GEOM.Flush+pim_vig_with_metadata$MYX.GEOM.Intestines  
MYX.GO<-(pim_vig_with_metadata$MYX.GO.CONNECTIVETISSUE+(2*pim_vig_with_metadata$MYX.GO.Eye)+
  as.numeric(pim_vig_with_metadata$MYX.GO.GONAD)+as.numeric(pim_vig_with_metadata$MYX.GO.Kidney))             
MYX.THEL<-(2*pim_vig_with_metadata$MYX.THEL.PectoralFin)+pim_vig_with_metadata$MYXO.THEL.CaudalFin+
  pim_vig_with_metadata$MYXO.THEL.ConnectiveTissue+as.numeric(pim_vig_with_metadata$MYXO.THEL.dorsalFIN)+
  (2*pim_vig_with_metadata$MYXO.THEL.GILL)+pim_vig_with_metadata$myxo.thel.Skin   
MYXO.BC<-(2*pim_vig_with_metadata$MYXO.BC.GILL)                  
MYXO.SBAD<-(pim_vig_with_metadata$MYXO.SBAD.ConnectiveTissue+(2*pim_vig_with_metadata$MYXO.SBAD.GILL))                   
MYXO.UNK<-pim_vig_with_metadata$MYXO.UNK.Connective.Tissue+(2*pim_vig_with_metadata$MYXO.UNK.GILL)                  
MYXO.UP<-(2*pim_vig_with_metadata$MYXO.UP.Eye)               
NEM.CAP<-pim_vig_with_metadata$NEM.CAP.Intestine                
NEM.CONA<-pim_vig_with_metadata$NEM.CONA.CONNECTIVETISSUES+pim_vig_with_metadata$NEM.CONA.Intestine               
NEM.DICH<-pim_vig_with_metadata$NEM.DICH.Intestine    
NEM.LARV<-pim_vig_with_metadata$NEM.LARV.Flush+pim_vig_with_metadata$NEM.LARV.Intestine+
  pim_vig_with_metadata$NEM.CYST.CONNECTIVETISSUE
NEM.TRBK<-pim_vig_with_metadata$NEM.TRBK.intestine               
NEM.UNK<-pim_vig_with_metadata$NEM.UNK.Intestine                
TREM.BUC<-pim_vig_with_metadata$trem.buc.AnalFin+pim_vig_with_metadata$trem.buc.CaudalFin+
  pim_vig_with_metadata$TREM.BUC.DorsalFin+(2*pim_vig_with_metadata$TREM.BUC.PectoralFin)            
TREM.MET<-pim_vig_with_metadata$TREM.MET.BodyCavity+pim_vig_with_metadata$TREM.MET.CaudalFin+
  pim_vig_with_metadata$TREM.MET.ConnectiveTissue+(2*pim_vig_with_metadata$TREM.MET.EYE)+
  pim_vig_with_metadata$TREM.MET.FLUSH+(2*pim_vig_with_metadata$TREM.MET.GILL)+
  pim_vig_with_metadata$TREM.MET.Intestine+as.numeric(pim_vig_with_metadata$TREM.MET.Kidney)+
  as.numeric(pim_vig_with_metadata$TREM.MET.LIVER)+as.numeric(pim_vig_with_metadata$TREM.MET.Spleen)+
  pim_vig_with_metadata$TREM.MET.Stomach   
TREM.METAF<-(2*pim_vig_with_metadata$TREM.MET.AF.PectoralFin)+
  pim_vig_with_metadata$TREM.METAF.AnalFin+pim_vig_with_metadata$trem.metaf.caudalfin+
  (2*pim_vig_with_metadata$TREM.METAF.EYE)+pim_vig_with_metadata$TREM.METAF.FLUSH+
  (2*pim_vig_with_metadata$TREM.METAF.GILL)+(2*as.numeric(pim_vig_with_metadata$TREM.METAF.PelvicFin))
TREM.METAGS<-pim_vig_with_metadata$TREM.METAGS.AnalFin+pim_vig_with_metadata$TREM.METAGS.CaudalFin+
  pim_vig_with_metadata$TREM.METAGS.CONNECTIVETISSUE+pim_vig_with_metadata$TREM.METAGS.FLUSH+
  (2*pim_vig_with_metadata$TREM.METAGS.GILL)+pim_vig_with_metadata$TREM.METAGS.Intestine+
  as.numeric(pim_vig_with_metadata$TREM.METAGS.Kidney)+as.numeric(pim_vig_with_metadata$TREM.METAGS.Liver)+
  pim_vig_with_metadata$TREM.METAGS.SKIN+pim_vig_with_metadata$TREM.METAGS.Stomach 
TREM.MS<-pim_vig_with_metadata$TREM.MS.FLUSH                    
TREM.NEA<-pim_vig_with_metadata$TREM.NEA.CONNECTIVETISSUE+pim_vig_with_metadata$TREM.NEA.FLUSH+
  pim_vig_with_metadata$TREM.NEA.Intestine+pim_vig_with_metadata$TREM.NEA.LIVER                   
TREM.RNDA<-pim_vig_with_metadata$TREM.RNDA.FLUSH+pim_vig_with_metadata$TREM.RNDA.Intestine              
TREM.SSS<-(2*pim_vig_with_metadata$TREM.SSS.GILL)
TREM.UNK<-pim_vig_with_metadata$TREM.UNK.CONNECTIVETISSUE+pim_vig_with_metadata$TREM.UNK.FLUSH
# not a parasite: WORM.A<-pim_vig_with_metadata$WORM.A.Flush                     



per_vig_processed_data<-cbind.data.frame(pim_vig_with_metadata$CatalogNumber,pim_vig_with_metadata$YearCollected.x,
                              pim_vig_with_metadata$MonthCollected.x,pim_vig_with_metadata$DayCollected.x,
                              pim_vig_with_metadata$IndividualFishID,pim_vig_with_metadata$Dissector_and_Examiner,
                              pim_vig_with_metadata$DissectionDate,pim_vig_with_metadata$Sex,
                              pim_vig_with_metadata$TotalLength_mm,pim_vig_with_metadata$StandardLength_mm,
                              pim_vig_with_metadata$Weight_mg,pim_vig_with_metadata$combo,
                              pim_vig_with_metadata$Latitude.y,
                              pim_vig_with_metadata$Longitude.y,ACANTH.BIGB,ACANTH.BKR,            
                              CEST.COMP,CEST.GODZ,CEST.L,CEST.NTG,CEST.PEA,CEST.UNK,CEST.VIT,COPE.CL,COPE.GRAB,
                              META.HS,META.UNK,MONO.DACT,MONO.GYRO,MONO.LG,MYX.GEOM,MYX.GO,MYX.THEL,
                              MYXO.BC,MYXO.SBAD,MYXO.UNK,MYXO.UP,NEM.CAP,NEM.CONA,NEM.DICH,NEM.CYST,NEM.LARV,NEM.TRBK,
                              NEM.UNK,TREM.BUC,TREM.MET,TREM.METAF,TREM.METAGS,TREM.MS,TREM.NEA,TREM.RNDA,
                              TREM.SSS,TREM.UNK)
# Name the columns

colnames(per_vig_processed_data)[1]<-"CatalogNumber"
colnames(per_vig_processed_data)[2]<-"YearCollected"
colnames(per_vig_processed_data)[3]<-"MonthCollected"
colnames(per_vig_processed_data)[4]<-"DayCollected"
colnames(per_vig_processed_data)[5]<-"IndividualFishID"
colnames(per_vig_processed_data)[6]<-"Dissector_and_Examiner"
colnames(per_vig_processed_data)[7]<-"DissectionDate"
colnames(per_vig_processed_data)[8]<-"Sex"
colnames(per_vig_processed_data)[9]<-"TotalLength_mm"
colnames(per_vig_processed_data)[10]<-"StandardLength_mm"
colnames(per_vig_processed_data)[11]<-"Weight_mg"
colnames(per_vig_processed_data)[12]<-"combo"
colnames(per_vig_processed_data)[13]<-"Latitude"
colnames(per_vig_processed_data)[14]<-"Longitude"


# Make the dataset analyzable

per_vig_processed_data_longer<-melt(per_vig_processed_data,id=c("CatalogNumber", "YearCollected", 
                                       "MonthCollected", "DayCollected", 
                                       "IndividualFishID", "Dissector_and_Examiner",
                                       "DissectionDate", "Sex", 
                                       "TotalLength_mm", "StandardLength_mm",
                                       "Weight_mg","combo","Latitude","Longitude"))

colnames(per_vig_processed_data_longer)[15]<-"psite_spp"
colnames(per_vig_processed_data_longer)[16]<-"psite_count"


# Export both sheets

write.csv(per_vig_processed_data_longer, file="data/processed/Pimephales_vigilax_processed_machine_readable.csv")
write.csv(per_vig_processed_data, file="data/processed/Pimephales_vigilax_processed_human_readable.csv")
