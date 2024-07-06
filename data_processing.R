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

pim_vig_today<-read.csv("data/raw/Pimephales_vigilax_Datasheet_2024.06.25.csv")
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
             #pim_vig_with_metadata$CYST.UNK.CONNECTIVETISSUES+pim_vig_with_metadata$CYST.UNK.DorsalFin+
             #(2*pim_vig_with_metadata$CYST.UNK.Eye)+pim_vig_with_metadata$CYST.UNK.Flush+
             #pim_vig_with_metadata$CYST.UNK.intestines+as.numeric(pim_vig_with_metadata$CYST.UNK.Kidney)+
             #pim_vig_with_metadata$CYST.UNK.LIVER+(2*pim_vig_with_metadata$CYST.UNK.PectoralFin)+
             #as.numeric(pim_vig_with_metadata$CYST.UNK.PelvicFin)+(2*pim_vig_with_metadata$cyst.unkn.gill))                   
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


# Need to fix the masses of fish, because they are one decimal point off.

pim_vig_with_metadata$TotalLength_mm<-(10*pim_vig_with_metadata$TotalLength_mm)

per_vig_processed_data<-cbind.data.frame(pim_vig_with_metadata$CatalogNumber,pim_vig_with_metadata$YearCollected.x,
                              pim_vig_with_metadata$MonthCollected.x,pim_vig_with_metadata$DayCollected.x,
                              pim_vig_with_metadata$IndividualFishID,pim_vig_with_metadata$Dissector_and_Examiner,
                              pim_vig_with_metadata$DissectionDate,pim_vig_with_metadata$Sex,
                              pim_vig_with_metadata$TotalLength_mm,pim_vig_with_metadata$StandardLength_mm,
                              pim_vig_with_metadata$Weight_mg, pim_vig_with_metadata$CI.y,
                              pim_vig_with_metadata$combo,
                              pim_vig_with_metadata$Latitude.y,
                              pim_vig_with_metadata$Longitude.y,ACANTH.BIGB,ACANTH.BKR,            
                              CEST.COMP,CEST.GODZ,CEST.L,CEST.NTG,CEST.PEA,CEST.UNK,CEST.VIT,COPE.CL,COPE.GRAB,
                              META.HS,META.UNK,MONO.DACT,MONO.GYRO,MONO.LG,MYX.GEOM,MYX.GO,MYX.THEL,
                              MYXO.BC,MYXO.SBAD,MYXO.UNK,MYXO.UP,NEM.CAP,NEM.CONA,NEM.DICH,NEM.LARV,NEM.TRBK,
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
colnames(per_vig_processed_data)[12]<-"CI"
colnames(per_vig_processed_data)[13]<-"combo"
colnames(per_vig_processed_data)[14]<-"Latitude"
colnames(per_vig_processed_data)[15]<-"Longitude"


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

write.csv(per_vig_processed_data_longer, file="data/processed/Pimephales_vigilax_processed_machine_readable_UPDATED_2024.07.06.csv")
write.csv(per_vig_processed_data, file="data/processed/Pimephales_vigilax_processed_human_readable_UPDATED_2024.07.06.csv")



### ICTPUN


# You can also do it the old-fashioned way

ict_pun_today<-read.csv("data/raw/Ictalurus_punctatus_Datasheet_2024.06.28.csv")
length(ict_pun_today$CatalogNumber)


# Now merge dissection data with meta-data (all.x makes sure that you keep all individual fish from the same lot, even though
# they have the same catalog number).

ict_pun_with_metadata<-merge(ict_pun_today, meta_data, by.x = "CatalogNumber", by.y = "CatalogNumber", all.x = TRUE)
length(ict_pun_with_metadata$CatalogNumber)


# Plot to see where the dissected fish fall

plot(jitter(ict_pun_with_metadata$Latitude.y,30)~jitter(ict_pun_with_metadata$YearCollected.y,5))+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)
str(ict_pun_with_metadata)
ict_pun_headers<-ls(ict_pun_with_metadata)
ict_pun_headers<-as.factor(ict_pun_headers)


# Add parasites across organs and double what needs to be doubled.

ACAN.OSTR<-ict_pun_with_metadata$ACAN.OSTR..Stomach+ict_pun_with_metadata$ACANTH.OSTR.INTESTINE            
CEST.COMP<-ict_pun_with_metadata$CEST.COMP.Intestine              
CEST.MEG<-ict_pun_with_metadata$CEST.MEG.Intestine
# not a parasite, ignore: CYST.BD<-ict_pun_with_metadata$CYST.BD.GILL                    
# not a parasite, ignore: CYST.MX<-ict_pun_with_metadata$CYST.MX.ConnectiveTissue+ict_pun_with_metadata$CYST.MX.GILL                     

# There are many, many cysts in ICTPUN, some of which contain no discernable structure within - no folded up worm,
# no myxozoans. We originally classified them as larval trematodes, but upon further inspection could not confirm
# that they were trematodes or, indeed, anything. So we created a new classification: CYST.NEH.  Anything that was
# originally classified as TREM.LARV.FIN or TREM.MID.FIN should be categorized as positive for CYST.NEH.
# Since we cannot confirm that these are actually parasites, I'm not tallying them here, but someone else could
# if they were interested in this category.
# CYST.NEH.ANAL.FIN               
# TREM.LARV.AnalFin                
# TREM.LARV.CaudalFin             
# TREM.LARV.DorsalFin       
# TREM.LARV.PectoralFin            
# TREM.LARV.PelvicFin   
# TREM.MID.ANALFIN                 
# TREM.MID.CAUDALFIN              
# TREM.MID.DORSALFIN               
# TREM.MID.PECTORALFIN             
# TREM.MID.PELVICFIN      

# not a parasite, ignore: CYST.UNK<-ict_pun_with_metadata$CYST.UNK..URINARY.BLADDER+
# ict_pun_with_metadata$CYST.UNK.Caudal+ict_pun_with_metadata$CYST.UNK.ConnectiveTissue+
# (2*ict_pun_with_metadata$CYST.UNK.Eye)+ict_pun_with_metadata$CYST.UNK.Flush+
# ict_pun_with_metadata$CYST.UNK.Intestine+ict_pun_with_metadata$CYST.UNK.Kidney+
# ict_pun_with_metadata$CYST.UNK.Liver+(2*ict_pun_with_metadata$CYST.UNK.Pectoral)+
# ict_pun_with_metadata$CYST.UNKN.GILL  

# Just to make sure that we're accounting properly: In the beginning, we found cysts in the kidney that we thought 
# were parasites, but now we know they are host tissue. We were also having trouble disambiguating kidney from
# liver tissue, so some of these might have been recorded as CYST.UNK.Liver instead of CYST.UNK.Kidney.
# We figured it out by the second day of dissection at TUBRI. So, make sure to put in NA for all CYST.UNK.Liver 
# and CYST.UNK.Kidney.

ict_pun_with_metadata$CYST.UNK.Kidney[ict_pun_with_metadata$DissectionDate<"6/26/2024"]<-NA
ict_pun_with_metadata$CYST.UNK.Liver[ict_pun_with_metadata$DissectionDate<"6/26/2024"]<-NA

LEECH.KUT<-ict_pun_with_metadata$LEECH.KUT.FLUSH                  
MONO.IP<-ict_pun_with_metadata$MONO.IP.FLUSH+(2*ict_pun_with_metadata$MONO.IP.GILL)                    
MONO.UNK<-(2*ict_pun_with_metadata$MONO.UNKN.GILL)           
MYX.GEOM<-ict_pun_with_metadata$MYX.GEOM.Intestine               
NEM.CYST<-ict_pun_with_metadata$CYST.NEM.ConnectiveTissue+ict_pun_with_metadata$NEM.CYST.Liver
NEM.NMTS<-ict_pun_with_metadata$NEM.NMTS.FLUSH+ict_pun_with_metadata$NEM.NMTS.Intestine+
  ict_pun_with_metadata$NEM.NMTS.Stomach                 
NEM.PHAR<-ict_pun_with_metadata$NEM.PHAR.Flush+(2*ict_pun_with_metadata$NEM.PHAR.GILL)+
  ict_pun_with_metadata$NEM.PHAR.INTESTINE+ict_pun_with_metadata$NEM.PHAR.Stomach                
NEM.PTY<-ict_pun_with_metadata$NEM.PTY.INTESTINE                
NEM.SP<-ict_pun_with_metadata$NEM.SP.INTESTINE                 
NEM.TAF<-ict_pun_with_metadata$NEM.TAF.Flush                   
NEM.UNK<-(2*ict_pun_with_metadata$NEM.UNK.gill)+ict_pun_with_metadata$NEM.UNK.Intestine                

# The first nematomorph we saw was probably a nematomorph. Then we saw another one that appeared to be much
# smaller and insdie of an insect. We marked it as a nematomorph. But later, we saw a similar insect, and its legs
# looked just like the "nematomorph" we thought we saw. The first, large nematomorph is therefore valid, but the 
# second, smaller one isn't. I've got to get rid of that one.

NMORPH<-ict_pun_with_metadata$NMORPH.INT.Intestine    

stuff<-ict_pun_with_metadata %>%
  filter(NMORPH.INT.Intestine>0)

# You want to replace IndividualFishID 157316_01

ict_pun_with_metadata$NMORPH.INT.Intestine[ict_pun_with_metadata$IndividualFishID=="157316_01"]<-"0"

stuff<-ict_pun_with_metadata %>%
  filter(NMORPH.INT.Intestine>0)

TREM.2<-ict_pun_with_metadata$TREM.2.FLUSH                    
TREM.ALLO<-ict_pun_with_metadata$TREM.ALLO.Flush+ict_pun_with_metadata$TREM.ALLO.Intestine+
  ict_pun_with_metadata$TREM.BLE.ConnectiveTissue       
TREM.BLE<-ict_pun_with_metadata$TREM.BLE.Intestine               
TREM.CREP<-ict_pun_with_metadata$TREM.CREP.INTESTINE              
TREM.F<-ict_pun_with_metadata$TREM.F.ConnectiveTIssue+ict_pun_with_metadata$TREM.F.STOMACH                   
TREM.META.GG<-(2*ict_pun_with_metadata$META.GG.gill)                

# There are some weird values in TREM.META.UNK - way too high for this parasite. I looked back at the datasheet,
# and these values should have been recorded for MYXO.TAIL.GILL. Fix it.

ict_pun_with_metadata$MYX.TAIL.GILL<-0

stuff<-ict_pun_with_metadata %>%
  filter(META.UNK.GILL>160)

ict_pun_with_metadata$META.UNK.GILL[ict_pun_with_metadata$IndividualFishID=="31542_01"]<-"0"
ict_pun_with_metadata$MYX.TAIL.GILL[ict_pun_with_metadata$IndividualFishID=="31542_01"]<-"161"

stuff<-ict_pun_with_metadata %>%
  filter(MYX.TAIL.GILL>160)

MYX.TAIL<-ict_pun_with_metadata$MYX.TAIL.Caudal+ict_pun_with_metadata$MYX.TAIL.Dorsal+
  (2*as.numeric(ict_pun_with_metadata$MYX.TAIL.GILL))

TREM.META.UNK<-(2*as.numeric(ict_pun_with_metadata$META.UNK.GILL))+ict_pun_with_metadata$TREM.META.FLUSH+
  ict_pun_with_metadata$TREM.UNK.ConnectiveTissue

stuff<-ict_pun_with_metadata %>%
  filter(META.UNK.GILL>60)

ict_pun_with_metadata$META.UNK.GILL[ict_pun_with_metadata$IndividualFishID=="182143_02"]<-"0"
ict_pun_with_metadata$MYX.TAIL.GILL[ict_pun_with_metadata$IndividualFishID=="182143_02"]<-"62"

stuff<-ict_pun_with_metadata %>%
  filter(META.UNK.GILL>0)

ict_pun_with_metadata$META.UNK.GILL[ict_pun_with_metadata$IndividualFishID=="157299_02"]<-"0"
ict_pun_with_metadata$MYX.TAIL.GILL[ict_pun_with_metadata$IndividualFishID=="157299_02"]<-"14"

stuff<-ict_pun_with_metadata %>%
  filter(META.UNK.GILL>0)

ict_pun_with_metadata$META.UNK.GILL[ict_pun_with_metadata$IndividualFishID=="198297_01"]<-"0"
ict_pun_with_metadata$MYX.TAIL.GILL[ict_pun_with_metadata$IndividualFishID=="198297_01"]<-"9"

stuff<-ict_pun_with_metadata %>%
  filter(META.UNK.GILL>0)

ict_pun_with_metadata$META.UNK.GILL[ict_pun_with_metadata$IndividualFishID=="40616_01"]<-"0"
ict_pun_with_metadata$MYX.TAIL.GILL[ict_pun_with_metadata$IndividualFishID=="40616_01"]<-"4"

stuff<-ict_pun_with_metadata %>%
  filter(META.UNK.GILL>0)

ict_pun_with_metadata$META.UNK.GILL[ict_pun_with_metadata$IndividualFishID=="153841_01"]<-"0"
ict_pun_with_metadata$MYX.TAIL.GILL[ict_pun_with_metadata$IndividualFishID=="153841_01"]<-"4"

stuff<-ict_pun_with_metadata %>%
  filter(META.UNK.GILL>0)

ict_pun_with_metadata$META.UNK.GILL[ict_pun_with_metadata$IndividualFishID=="99517_01"]<-"0"
ict_pun_with_metadata$MYX.TAIL.GILL[ict_pun_with_metadata$IndividualFishID=="99517_01"]<-"3"

stuff<-ict_pun_with_metadata %>%
  filter(META.UNK.GILL>0)

ict_pun_with_metadata$META.UNK.GILL[ict_pun_with_metadata$IndividualFishID=="160020_02"]<-"0"
ict_pun_with_metadata$MYX.TAIL.GILL[ict_pun_with_metadata$IndividualFishID=="160020_02"]<-"3"

stuff<-ict_pun_with_metadata %>%
  filter(META.UNK.GILL>0)

# The last remaining one is real!

TREM.META.UNK<-(2*as.numeric(ict_pun_with_metadata$META.UNK.GILL))+ict_pun_with_metadata$TREM.META.FLUSH+
  ict_pun_with_metadata$TREM.UNK.ConnectiveTissue
TREM.LG<-ict_pun_with_metadata$TREM.LG.DORSALFIN+ict_pun_with_metadata$TREM.LG.FIN.Anal.Fin+
  ict_pun_with_metadata$TREM.LG.FIN.Caudal+(2*ict_pun_with_metadata$TREM.LG.FIN.Pelvic.Fin)+
  (2*ict_pun_with_metadata$TREM.LG.PECTORALFIN)            
TREM.LIP<-ict_pun_with_metadata$TREM.LIP.Intestine              
TREM.MEGICT<-ict_pun_with_metadata$TREM.MEGICT.FLUSH+ict_pun_with_metadata$TREM.MEGICT.Intestine+
  ict_pun_with_metadata$TREM.MEGICT.Stomach             
TREM.SMILE<-ict_pun_with_metadata$TREM.SMILE.Intestine             
TREM.UNK<-ict_pun_with_metadata$TREM.UNKN.Intestine                    
#not a parasite, ignore: WORM.UNK<-ict_pun_with_metadata$WORM.UNK.B.intestine+ict_pun_with_metadata$WORM.UNK.Stomach


ict_pun_processed_data<-cbind.data.frame(ict_pun_with_metadata$CatalogNumber,ict_pun_with_metadata$YearCollected.x,
                                         ict_pun_with_metadata$MonthCollected.x,ict_pun_with_metadata$DayCollected.x,
                                         ict_pun_with_metadata$IndividualFishID,ict_pun_with_metadata$Dissector,
                                         ict_pun_with_metadata$DissectionDate,ict_pun_with_metadata$Sex,
                                         ict_pun_with_metadata$TotalLength_mm,ict_pun_with_metadata$StandardLength_mm,
                                         ict_pun_with_metadata$weight_mg,ict_pun_with_metadata$combo,
                                         ict_pun_with_metadata$Latitude.y,
                                         ict_pun_with_metadata$Longitude.y,ACAN.OSTR, CEST.COMP, CEST.MEG, 
                                         LEECH.KUT, MONO.IP, MONO.UNK, MYX.GEOM, NEM.CYST, NEM.NMTS, NEM.PHAR, 
                                         NEM.PTY, NEM.SP, NEM.TAF, NEM.UNK, NMORPH, TREM.2, TREM.ALLO, TREM.BLE, 
                                         TREM.CREP, TREM.F, TREM.META.GG, MYX.TAIL, TREM.META.UNK, 
                                         TREM.LG, TREM.LIP, TREM.MEGICT, TREM.SMILE, TREM.UNK)

# Name the columns

colnames(ict_pun_processed_data)[1]<-"CatalogNumber"
colnames(ict_pun_processed_data)[2]<-"YearCollected"
colnames(ict_pun_processed_data)[3]<-"MonthCollected"
colnames(ict_pun_processed_data)[4]<-"DayCollected"
colnames(ict_pun_processed_data)[5]<-"IndividualFishID"
colnames(ict_pun_processed_data)[6]<-"Dissector_and_Examiner"
colnames(ict_pun_processed_data)[7]<-"DissectionDate"
colnames(ict_pun_processed_data)[8]<-"Sex"
colnames(ict_pun_processed_data)[9]<-"TotalLength_mm"
colnames(ict_pun_processed_data)[10]<-"StandardLength_mm"
colnames(ict_pun_processed_data)[11]<-"Weight_mg"
colnames(ict_pun_processed_data)[12]<-"combo"
colnames(ict_pun_processed_data)[13]<-"Latitude"
colnames(ict_pun_processed_data)[14]<-"Longitude"

View(ict_pun_processed_data)


# It looks like meta-data are missing for some of the fish.

stuff<-ict_pun_processed_data %>%
  filter(IndividualFishID=="0")

stuff<-ict_pun_processed_data %>%
  filter(CatalogNumber=="159967")

ict_pun_processed_data$IndividualFishID[ict_pun_processed_data$CatalogNumber=="159967" & 
                                          ict_pun_processed_data$TotalLength_mm=="64"]<-"159967_01"

stuff<-ict_pun_processed_data %>%
  filter(is.na(YearCollected))

ict_pun_processed_data$YearCollected[ict_pun_processed_data$IndividualFishID=="63780_01"]<-"1970"

stuff<-ict_pun_processed_data %>%
  filter(is.na(YearCollected))

stuff<-ict_pun_processed_data %>%
  filter(is.na(combo))

stuff<-meta_data %>%
  filter(CatalogNumber=="62099")

# Someone recorded the Catalog Number incorrectly

ict_pun_processed_data$CatalogNumber[ict_pun_processed_data$CatalogNumber=="62099"]<-"62094"

stuff<-meta_data %>%
  filter(CatalogNumber=="62099")

ict_pun_processed_data$Latitude[ict_pun_processed_data$CatalogNumber=="62094"]<-"30.77861"
ict_pun_processed_data$Longitude[ict_pun_processed_data$CatalogNumber=="62094"]<-"-89.82972"
ict_pun_processed_data$combo[ict_pun_processed_data$CatalogNumber=="62094"]<-"control_1964_1973"
ict_pun_processed_data$IndividualFishID[ict_pun_processed_data$CatalogNumber=="62094" & 
                                          ict_pun_processed_data$TotalLength_mm=="77"]<-"62094_02"

stuff<-ict_pun_processed_data %>%
  filter(CatalogNumber=="62094")

stuff<-ict_pun_processed_data %>%
  filter(is.na(combo))

stuff<-meta_data %>%
  filter(CatalogNumber=="157229")

# Someone recorded the Catalog Number incorrectly

ict_pun_processed_data$CatalogNumber[ict_pun_processed_data$CatalogNumber=="157229"]<-"157299"

stuff<-meta_data %>%
  filter(CatalogNumber=="157299")

ict_pun_processed_data$Latitude[ict_pun_processed_data$CatalogNumber=="157299"]<-"30.77583"
ict_pun_processed_data$Longitude[ict_pun_processed_data$CatalogNumber=="157299"]<-"-89.82722"
ict_pun_processed_data$combo[ict_pun_processed_data$CatalogNumber=="157299"]<-"control_1984-1993"

stuff<-ict_pun_processed_data %>%
  filter(CatalogNumber=="157299")

stuff<-ict_pun_processed_data %>%
  filter(is.na(combo))





# Make the dataset analyzable

ict_pun_processed_data_longer<-melt(ict_pun_processed_data,id=c("CatalogNumber", "YearCollected", 
                                                                "MonthCollected", "DayCollected", 
                                                                "IndividualFishID", "Dissector_and_Examiner",
                                                                "DissectionDate", "Sex", 
                                                                "TotalLength_mm", "StandardLength_mm",
                                                                "Weight_mg","combo","Latitude","Longitude"))

colnames(ict_pun_processed_data_longer)[15]<-"psite_spp"
colnames(ict_pun_processed_data_longer)[16]<-"psite_count"


# Export both sheets

write.csv(ict_pun_processed_data_longer, file="data/processed/Ictalurus_punctatus_processed_machine_readable.csv")
write.csv(ict_pun_processed_data, file="data/processed/Ictalurus_punctatus_processed_human_readable.csv")
