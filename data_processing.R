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



### PIMVIG: Pimephales vigilax----

# Download latest file # GET THIS WORKING LATER
# drive_download(as_id("16taxgLRu1-d_t8lT9mHVWOtxZ4MfQPvU7UBWHvLRfaI"), 
#               path = "data/raw/pim_vig.csv", overwrite = TRUE)


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

### MYX.GO were found in many organs. Sometimes, those were organs that were not found in all fish
### (e.g., gonads). However, even if the gonads (or spleen, or heart) were too small to ID, we feel confident that
### we would have seen myxozoan cysts in those organs, even if the organs themselves were just perceived as part 
### of the visceral mush inside of fish. Therefore, it wouldn't be appropriate to propagate NAs across the entire
### MYX.GO count within a fish (across organs). The code below allows the sum of MYX.GO within an individual
### fish ignore the NAs, and therefore provide a total number of MYX.GO cysts regardless of whether one of the organs
### was "missing".

pim_vig_with_metadata$MYX.GO.Eye.DUB<-2*pim_vig_with_metadata$MYX.GO.Eye
pim_vig_with_metadata$MYX.GO.GONAD<-as.numeric(pim_vig_with_metadata$MYX.GO.GONAD)
pim_vig_with_metadata$MYX.GO.Kidney<-as.numeric(pim_vig_with_metadata$MYX.GO.Kidney)

MYX.GO<-rowSums(pim_vig_with_metadata[,c("MYX.GO.CONNECTIVETISSUE","MYX.GO.Eye.DUB",
                                            "MYX.GO.GONAD","MYX.GO.Kidney")], na.rm=TRUE)
            
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


### TREM.MET were found in many organs. Sometimes, those were organs that were not found in all fish
### (e.g., gonads). However, even if the gonads (or spleen, or heart) were too small to ID, we feel confident that
### we would have seen metacercariae in those organs, even if the organs themselves were just perceived as part 
### of the visceral mush inside of fish. Therefore, it wouldn't be appropriate to propagate NAs across the entire
### TREM.MET count within a fish (across organs). The code below allows the sum of TREM.MET within an individual
### fish ignore the NAs, and therefore provide a total number of metacercariae regardless of whether one of the organs
### was "missing".

pim_vig_with_metadata$TREM.MET.EYE.DUB<-2*pim_vig_with_metadata$TREM.MET.EYE
pim_vig_with_metadata$TREM.MET.GILL.DUB<-2*pim_vig_with_metadata$TREM.MET.GILL
pim_vig_with_metadata$TREM.MET.Kidney<-as.numeric(pim_vig_with_metadata$TREM.MET.Kidney)
pim_vig_with_metadata$TREM.MET.LIVER<-as.numeric(pim_vig_with_metadata$TREM.MET.LIVER)
pim_vig_with_metadata$TREM.MET.Spleen<-as.numeric(pim_vig_with_metadata$TREM.MET.Spleen)

TREM.MET<-rowSums(pim_vig_with_metadata[,c("TREM.MET.BodyCavity","TREM.MET.CaudalFin",
                                         "TREM.MET.ConnectiveTissue","TREM.MET.EYE.DUB",
                                         "TREM.MET.FLUSH","TREM.MET.GILL","TREM.MET.Intestine",
                                         "TREM.MET.Kidney","TREM.MET.LIVER","TREM.MET.Spleen",
                                         "TREM.MET.Stomach")], na.rm=TRUE)

TREM.METAF<-(2*pim_vig_with_metadata$TREM.MET.AF.PectoralFin)+
  pim_vig_with_metadata$TREM.METAF.AnalFin+pim_vig_with_metadata$trem.metaf.caudalfin+
  (2*pim_vig_with_metadata$TREM.METAF.EYE)+pim_vig_with_metadata$TREM.METAF.FLUSH+
  (2*pim_vig_with_metadata$TREM.METAF.GILL)+(2*as.numeric(pim_vig_with_metadata$TREM.METAF.PelvicFin))

### TREM.METAGS were found in many organs. Sometimes, those were organs that were not found in all fish
### (e.g., gonads). However, even if the gonads (or spleen, or heart) were too small to ID, we feel confident that
### we would have seen metacercariae in those organs, even if the organs themselves were just perceived as part 
### of the visceral mush inside of fish. Therefore, it wouldn't be appropriate to propagate NAs across the entire
### TREM.METAGS count within a fish (across organs). The code below allows the sum of TREM.METAGS within an individual
### fish to ignore the NAs, and therefore provide a total number of metacercariae regardless of whether one of the organs
### was "missing".

pim_vig_with_metadata$TREM.METAGS.GILL.DUB<-2*pim_vig_with_metadata$TREM.METAGS.GILL
pim_vig_with_metadata$TREM.METAGS.Kidney<-as.numeric(pim_vig_with_metadata$TREM.METAGS.Kidney)
pim_vig_with_metadata$TREM.METAGS.Liver<-as.numeric(pim_vig_with_metadata$TREM.METAGS.Liver)

TREM.METAGS<-rowSums(pim_vig_with_metadata[,c("TREM.METAGS.AnalFin","TREM.METAGS.CaudalFin",
                                           "TREM.METAGS.CONNECTIVETISSUE","TREM.METAGS.FLUSH",
                                           "TREM.METAGS.GILL.DUB","TREM.METAGS.Intestine","TREM.METAGS.Kidney",
                                           "TREM.METAGS.Liver","TREM.METAGS.SKIN","TREM.METAGS.Stomach")], na.rm=TRUE)

TREM.MS<-pim_vig_with_metadata$TREM.MS.FLUSH                    
TREM.NEA<-pim_vig_with_metadata$TREM.NEA.CONNECTIVETISSUE+pim_vig_with_metadata$TREM.NEA.FLUSH+
  pim_vig_with_metadata$TREM.NEA.Intestine+pim_vig_with_metadata$TREM.NEA.LIVER                   
TREM.RNDA<-pim_vig_with_metadata$TREM.RNDA.FLUSH+pim_vig_with_metadata$TREM.RNDA.Intestine              
TREM.SSS<-(2*pim_vig_with_metadata$TREM.SSS.GILL)
TREM.UNK<-pim_vig_with_metadata$TREM.UNK.CONNECTIVETISSUE+pim_vig_with_metadata$TREM.UNK.FLUSH
# not a parasite: WORM.A<-pim_vig_with_metadata$WORM.A.Flush                     


# Need to fix the masses of fish, because they are one decimal point off.

pim_vig_with_metadata$Weight_mg<-(10*as.numeric(pim_vig_with_metadata$Weight_mg))

# Evaluate relationship between total length and weight.
ggplot(pim_vig_with_metadata,aes(TotalLength_mm,Weight_mg))+
  geom_point(size=4)

# Evaluate ratio between total length and weight. Divide length by weight
pim_vig_with_metadata$LenWeig <- pim_vig_with_metadata$TotalLength_mm/pim_vig_with_metadata$Weight_mg

# According to the following histogram the ratio between length and weight is somewhere between 0-0.30
hist(pim_vig_with_metadata$LenWeig)

# Now put it all together

pim_vig_processed_data<-cbind.data.frame(pim_vig_with_metadata$CatalogNumber,pim_vig_with_metadata$YearCollected.x,
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

colnames(pim_vig_processed_data)[1]<-"CatalogNumber"
colnames(pim_vig_processed_data)[2]<-"YearCollected"
colnames(pim_vig_processed_data)[3]<-"MonthCollected"
colnames(pim_vig_processed_data)[4]<-"DayCollected"
colnames(pim_vig_processed_data)[5]<-"IndividualFishID"
colnames(pim_vig_processed_data)[6]<-"Dissector_and_Examiner"
colnames(pim_vig_processed_data)[7]<-"DissectionDate"
colnames(pim_vig_processed_data)[8]<-"Sex"
colnames(pim_vig_processed_data)[9]<-"TotalLength_mm"
colnames(pim_vig_processed_data)[10]<-"StandardLength_mm"
colnames(pim_vig_processed_data)[11]<-"Weight_mg"
colnames(pim_vig_processed_data)[12]<-"CI"
colnames(pim_vig_processed_data)[13]<-"combo"
colnames(pim_vig_processed_data)[14]<-"Latitude"
colnames(pim_vig_processed_data)[15]<-"Longitude"

View(pim_vig_processed_data)


# Looks like meta-data (CI and combo) are missing for four of the fish

# It looks like meta-data are missing for some of the fish.

stuff<-pim_vig_processed_data %>%
  filter(is.na(combo))

pim_vig_processed_data$CI[pim_vig_processed_data$CatalogNumber=="157327"]<-"impact"
pim_vig_processed_data$combo[pim_vig_processed_data$CatalogNumber=="157327"]<-"impact_1984-1993"

stuff<-pim_vig_processed_data %>%
  filter(is.na(Latitude))

# Found the lat/long in iDigBio

pim_vig_processed_data$Latitude[pim_vig_processed_data$CatalogNumber=="157327"]<-30.19778
pim_vig_processed_data$Longitude[pim_vig_processed_data$CatalogNumber=="157327"]<--89.65611



# Make the dataset analyzable

pim_vig_processed_data_longer<-melt(pim_vig_processed_data,id=c("CatalogNumber", "YearCollected", 
                                       "MonthCollected", "DayCollected", 
                                       "IndividualFishID", "Dissector_and_Examiner",
                                       "DissectionDate", "Sex", 
                                       "TotalLength_mm", "StandardLength_mm",
                                       "Weight_mg","CI","combo","Latitude","Longitude"))

colnames(pim_vig_processed_data_longer)[16]<-"psite_spp"
colnames(pim_vig_processed_data_longer)[17]<-"psite_count"


# Export both sheets

write.csv(pim_vig_processed_data_longer, 
          file="data/processed/Pimephales_vigilax_processed_machine_readable_UPDATED_2024.07.10.csv")
write.csv(pim_vig_processed_data, 
          file="data/processed/Pimephales_vigilax_processed_human_readable_UPDATED_2024.07.10.csv")




### ICTPUN: Ictalurus punctatus----


# You can also do it the old-fashioned way

ict_pun_today<-read.csv("data/raw/Ictalurus_punctatus_Datasheet_2024.07.19.csv")
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
TREM.PHYLS<-ict_pun_with_metadata$TREM.PHYLS.URINARYBLADDER
TREM.SMILE<-ict_pun_with_metadata$TREM.SMILE.Intestine             
TREM.UNK<-ict_pun_with_metadata$TREM.UNKN.Intestine                    
#not a parasite, ignore: WORM.UNK<-ict_pun_with_metadata$WORM.UNK.B.intestine+ict_pun_with_metadata$WORM.UNK.Stomach


# Need to fix the masses of fish, because they are one decimal point off.

ict_pun_with_metadata$weight_mg<-(10*as.numeric(ict_pun_with_metadata$weight_mg))

#Revise relationship between fish size and weight
ggplot(ict_pun_with_metadata,aes(TotalLength_mm,weight_mg))+geom_point(size=4)

# Evaluate ratio between total length and weight. Divide length by weight
ict_pun_with_metadata$LenWeig <- ict_pun_with_metadata$TotalLength_mm/ict_pun_with_metadata$weight_mg

# According to the following histogram the ratio between length and weight is somewhere between 0-0.015
hist(ict_pun_with_metadata$LenWeig)

# Now put it all together

ict_pun_processed_data<-cbind.data.frame(ict_pun_with_metadata$CatalogNumber,ict_pun_with_metadata$YearCollected.x,
                                         ict_pun_with_metadata$MonthCollected.x,ict_pun_with_metadata$DayCollected.x,
                                         ict_pun_with_metadata$IndividualFishID,ict_pun_with_metadata$Dissector,
                                         ict_pun_with_metadata$DissectionDate,ict_pun_with_metadata$Sex,
                                         ict_pun_with_metadata$TotalLength_mm,ict_pun_with_metadata$StandardLength_mm,
                                         ict_pun_with_metadata$weight_mg,ict_pun_with_metadata$CI.y,
                                         ict_pun_with_metadata$combo,
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
colnames(ict_pun_processed_data)[12]<-"CI"
colnames(ict_pun_processed_data)[13]<-"combo"
colnames(ict_pun_processed_data)[14]<-"Latitude"
colnames(ict_pun_processed_data)[15]<-"Longitude"

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


stuff<-ict_pun_processed_data %>%
  filter(is.na(CI))

ict_pun_processed_data$CI[ict_pun_processed_data$CatalogNumber=="62094"]<-"control"
ict_pun_processed_data$CI[ict_pun_processed_data$CatalogNumber=="157299"]<-"control"

stuff<-ict_pun_processed_data %>%
  filter(is.na(CI))


# Make the dataset analyzable

ict_pun_processed_data_longer<-melt(ict_pun_processed_data,id=c("CatalogNumber", "YearCollected", 
                                                                "MonthCollected", "DayCollected", 
                                                                "IndividualFishID", "Dissector_and_Examiner",
                                                                "DissectionDate", "Sex", 
                                                                "TotalLength_mm", "StandardLength_mm",
                                                                "Weight_mg","CI","combo","Latitude","Longitude"))

colnames(ict_pun_processed_data_longer)[16]<-"psite_spp"
colnames(ict_pun_processed_data_longer)[17]<-"psite_count"


# Export both sheets

write.csv(ict_pun_processed_data_longer, 
          file="data/processed/Ictalurus_punctatus_processed_machine_readable_UPDATED_2024.07.19.csv")
write.csv(ict_pun_processed_data, 
          file="data/processed/Ictalurus_punctatus_processed_human_readable.csv_UPDATED_2024.07.19.csv")




### NOTATH: Notropis atherinoides----


# You can also do it the old-fashioned way

not_ath_today<-read.csv("data/raw/Notropis_atherinoides_Datasheet_2024.01.10.csv")
length(not_ath_today$CatalogNumber)


# Now merge dissection data with meta-data (all.x makes sure that you keep all individual fish from the same lot, even though
# they have the same catalog number).

not_ath_with_metadata<-merge(not_ath_today, meta_data, by.x = "CatalogNumber", by.y = "CatalogNumber", all.x = TRUE)
length(not_ath_with_metadata$CatalogNumber)


# Plot to see where the dissected fish fall

plot(jitter(not_ath_with_metadata$Latitude.y,30)~jitter(not_ath_with_metadata$YearCollected.y,5))+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)
str(not_ath_with_metadata)
not_ath_headers<-ls(not_ath_with_metadata)
not_ath_headers<-as.factor(not_ath_headers)


# Add parasites across organs and double what needs to be doubled.

CEST.COMP<-not_ath_with_metadata$CEST.COMP.INTESTINE              
# not a parasite: CEST.CYRC.INTESTINE              
CEST.HYD<-not_ath_with_metadata$CEST.HYD.INTESTINE               
CEST.LVS<-not_ath_with_metadata$CEST.LVS.INTESTINE              
CEST.TRI<-not_ath_with_metadata$cest.tri.intestine               
COPE.A<-(2*not_ath_with_metadata$COPE.A.GILL)                  
# not a parasite: CYST.SAC.Skin                    
# not a parasite: CYST.UNK.ANALFIN                 
# not a parasite: CYST.UNK.CAUDALFIN              
# not a parasite: CYST.UNK.GILL                    
# not a parasite: CYST.UNK.Kidney                  
# not a parasite: CYST.UNKN.ConnectiveTissue      
# not a parasite: CYST.UNKN.DORSALFIN              
# not a parasite: CYST.UNKN.LIVER                  
# not a parasite: CYST.UNKN.MOUNT                 
# not a parasite: CYST.UNKN.PECTORALFIN            
# not a parasite: CYST.UNKN.PelvicFin              
# not a parasite: cyst.unkn.skin                  
# not a parasite: CYST.UNKN.VOUCH                                 

### This monogene group is provisional. DACT and NA are definitely the same thing. BULL is 
### also probably NA, but we saw it at a weird angle? Mounting will distinguish among them. 
MONO.ALL<-(not_ath_with_metadata$MONO.BULL.FLUSH+(2*not_ath_with_metadata$MONO.BULL.GILL)+
             (2*not_ath_with_metadata$MONO.DACT.GILL)+(2*not_ath_with_metadata$MONO.NA.Gill))

MONO.UNK<-(2*not_ath_with_metadata$MONO.UNKN.GILL)
MYX.ROUND<-(2*not_ath_with_metadata$MYX.ROUND.skin)             
MYX.SP<-(2*not_ath_with_metadata$MYX.SP.GILL)                
NEM.BRAN<-(not_ath_with_metadata$NEM.BRAN.FLUSH+not_ath_with_metadata$NEM.BRAN.INTESTINE)              
NEM.CAP<-not_ath_with_metadata$NEM.CAP.INTESTINE                
NEM.CYST<-(not_ath_with_metadata$NEM.CYST.CONNECTIVETISSUE+not_ath_with_metadata$NEM.CYST.Intestine+
             not_ath_with_metadata$nem.cystL.Liver)
NEM.SPKY<-not_ath_with_metadata$NEM.SPKY.INTESTINE              
NEM.TFS<-not_ath_with_metadata$NEM.TFS.INTESTINE                
NEM.TWAS<-not_ath_with_metadata$NEM.TWAS.FLUSH                   
NEM.UNKN<-not_ath_with_metadata$NEM.UNKN.INTESTINE               
TREM.CYST<-(not_ath_with_metadata$TREM.CYST.CONNECTIVETISSUE+(2*not_ath_with_metadata$TREM.CYST.EYE)+
              not_ath_with_metadata$TREM.CYST.FLUSH+(2*not_ath_with_metadata$TREM.CYST.GILL))               
TREM.GB<-not_ath_with_metadata$TREM.GB.INTESTINE               

### Katie confirms that TREM.L was only seen once, but is definitely it's own thing.
TREM.L<-not_ath_with_metadata$TREM.L.INTESTINE                

TREM.LARV<-(not_ath_with_metadata$TREM.LARV.AnalFin+not_ath_with_metadata$trem.larv.CaudalFin+
              not_ath_with_metadata$TREM.LARV.DorsalFin+(2*not_ath_with_metadata$TREM.LARV.PectoralFin))           

### Larval trematodes were found in many organs. Sometimes, those were organs that were not found in all fish
### (e.g., gonads). However, even if the gonads (or spleen, or heart) were too small to ID, we feel confident that
### we would have seen metacercariae in those organs, even if the organs themselves were just perceived as part 
### of the visceral mush inside of fish. Therefore, it wouldn't be appropriate to propagate NAs across the entire
### TREM.META count within a fish (across organs). The code below allows the sum of TREM.META within an individual
### fish ignore the NAs, and therefore provide a total number of TREM.META regardless of whether one of the organs
### was "missing".

not_ath_with_metadata$TREM.META_4.EYE.DUB<-2*not_ath_with_metadata$TREM.META_4.EYE
not_ath_with_metadata$TREM.META_4.GILL.DUB<-2*not_ath_with_metadata$TREM.META_4.GILL
not_ath_with_metadata$TREM.META_4.SKIN.DUB<-2*not_ath_with_metadata$TREM.META_4.SKIN
not_ath_with_metadata$TREM.META_5.EYE.DUB<-2*not_ath_with_metadata$TREM.META_5.EYE
not_ath_with_metadata$TREM.META_6.Eye.DUB<-2*not_ath_with_metadata$TREM.META_6.Eye
not_ath_with_metadata$TREM.META_7.EYE.DUB<-2*not_ath_with_metadata$TREM.META_7.EYE
not_ath_with_metadata$TREM.META.8.EYE.DUB<-2*not_ath_with_metadata$TREM.META.8.EYE
not_ath_with_metadata$TREM.META.EYE.DUB<-2*not_ath_with_metadata$TREM.META.EYE
not_ath_with_metadata$TREM.META.GILL.DUB<-2*not_ath_with_metadata$TREM.META.GILL
not_ath_with_metadata$TREM.META.PELVICFIN.DUB<-2*not_ath_with_metadata$TREM.META.PELVICFIN
not_ath_with_metadata$TREM.META.GONAD<-as.numeric(not_ath_with_metadata$TREM.META.GONAD)
not_ath_with_metadata$TREM.META.Kidney<-as.numeric(not_ath_with_metadata$TREM.META.Kidney)

TREM.META<-rowSums(not_ath_with_metadata[,c("TREM.META_1.FLUSH","TREM.META_2.CONNECTIVETISSUES",
                                            "TREM.META_3.CONNECTIVETISSUES",
                                            "TREM.META_4.EYE.DUB","TREM.META_4.GILL.DUB","TREM.META_4.SKIN.DUB",
                                            "TREM.META_5.EYE.DUB","TREM.META_6.Eye.DUB","TREM.META_7.EYE.DUB",
                                            "TREM.META_9.CONNECTIVETISSUE","TREM.META.8.EYE.DUB",
                                            "TREM.META.CAUDALFIN",
                                            "TREM.META.CONNECTIVETISSUE","TREM.META.EYE.DUB","TREM.META.FLUSH",
                                            "TREM.META.GILL.DUB","TREM.META.GONAD","TREM.META.Kidney",
                                            "TREM.META.LIVER",
                                            "TREM.META.PELVICFIN.DUB")], na.rm=TRUE)
            
# not a parasite: CYSTNEH<-(2*not_ath_with_metadata$TREM.SMLARV.PECTORALFIN)       
# not a parasite: UNK.GC.INTESTINE                 
# not a parasite: UNK.SPK.INTESTINE


# Now put it all together

not_ath_processed_data<-cbind.data.frame(not_ath_with_metadata$CatalogNumber,not_ath_with_metadata$YearCollected.x,
                                         not_ath_with_metadata$MonthCollected.x,not_ath_with_metadata$DayCollected.x,
                                         not_ath_with_metadata$IndividualFishID,not_ath_with_metadata$Dissector,
                                         not_ath_with_metadata$DissectionDate,not_ath_with_metadata$Sex,
                                         not_ath_with_metadata$TotalLength_mm,not_ath_with_metadata$StandardLength_mm,
                                         not_ath_with_metadata$weight_mg,not_ath_with_metadata$CI.y,
                                         not_ath_with_metadata$combo,
                                         not_ath_with_metadata$Latitude.y,
                                         not_ath_with_metadata$Longitude.y,CEST.COMP,CEST.HYD,CEST.LVS,CEST.TRI,
                                         COPE.A,MONO.ALL,MONO.UNK,MYX.ROUND,MYX.SP,NEM.BRAN,NEM.CAP,NEM.CYST,
                                         NEM.SPKY,NEM.TFS,
                                         NEM.TWAS,NEM.UNKN,TREM.CYST,TREM.GB,TREM.L,TREM.LARV,TREM.META)

# Name the columns

colnames(not_ath_processed_data)[1]<-"CatalogNumber"
colnames(not_ath_processed_data)[2]<-"YearCollected"
colnames(not_ath_processed_data)[3]<-"MonthCollected"
colnames(not_ath_processed_data)[4]<-"DayCollected"
colnames(not_ath_processed_data)[5]<-"IndividualFishID"
colnames(not_ath_processed_data)[6]<-"Dissector_and_Examiner"
colnames(not_ath_processed_data)[7]<-"DissectionDate"
colnames(not_ath_processed_data)[8]<-"Sex"
colnames(not_ath_processed_data)[9]<-"TotalLength_mm"
colnames(not_ath_processed_data)[10]<-"StandardLength_mm"
colnames(not_ath_processed_data)[11]<-"Weight_mg"
colnames(not_ath_processed_data)[12]<-"CI"
colnames(not_ath_processed_data)[13]<-"combo"
colnames(not_ath_processed_data)[14]<-"Latitude"
colnames(not_ath_processed_data)[15]<-"Longitude"

View(not_ath_processed_data)

#Remove #VALUE! entries from Weight_mg with "NA"
not_ath_processed_data <- not_ath_processed_data %>%
  mutate(Weight_mg = recode(Weight_mg,
                           "#VALUE!" = "NA"))

#Evaluate relationship between TotalLength_mm and Weight_mg. There is evidently a mistake. Some values must be multiplied by 10.
ggplot(not_ath_processed_data,aes(TotalLength_mm,Weight_mg))+geom_point(size=4)

# Evaluate ratio between total length and weight to inform decision about which numbers must be multiplied by ten. For this, divide length by weight
not_ath_processed_data$LenWeig <- as.numeric(not_ath_processed_data$Weight_mg)/not_ath_processed_data$TotalLength_mm


# All incorrect weights that must be multiplied by 10 have a ratio below 3 which does not make sense. See the data below to confirm. The weights have a decimal point which should not be the case.
lenweighb <- subset(not_ath_processed_data, LenWeig<3) 
View(lenweighb)

# Therefore, multiply all weights with LenWeig higher than 3 by 10

not_ath_processed_data$Weight_mg <- ifelse(not_ath_processed_data$LenWeig < 3,                # condition
                                                       as.numeric(not_ath_processed_data$Weight_mg)*10,    # what if condition is TRUE
                                                       as.numeric(not_ath_processed_data$Weight_mg)       # what if condition is FALSE
)


#Re-evaluate relationship between TotalLength_mm and Weight_mg. Now it is fixed.
ggplot(not_ath_processed_data,aes(TotalLength_mm,Weight_mg))+geom_point(size=4)


# Make the dataset analyzable

not_ath_processed_data_longer<-melt(not_ath_processed_data,id=c("CatalogNumber", "YearCollected", 
                                                                "MonthCollected", "DayCollected", 
                                                                "IndividualFishID", "Dissector_and_Examiner",
                                                                "DissectionDate", "Sex", 
                                                                "TotalLength_mm", "StandardLength_mm",
                                                                "Weight_mg","CI","combo","Latitude","Longitude"))

colnames(not_ath_processed_data_longer)[16]<-"psite_spp"
colnames(not_ath_processed_data_longer)[17]<-"psite_count"

# Export both sheets

write.csv(not_ath_processed_data_longer, file="data/processed/Notropis_atherinoides_processed_machine_readable.csv")
write.csv(not_ath_processed_data, file="data/processed/Notropis_atherinoides_processed_human_readable.csv")




### HYBNUC: Hybognathus nuchalis----


# You can also do it the old-fashioned way

hyb_nuc_today<-read.csv("data/raw/Hybognathus_nuchalis_Datasheet_2024.07.25.csv")
length(hyb_nuc_today$CatalogNumber)


# Now merge dissection data with meta-data (all.x makes sure that you keep all individual fish from the same lot, even though
# they have the same catalog number).

hyb_nuc_today_with_metadata<-merge(hyb_nuc_today, meta_data, by.x = "CatalogNumber", by.y = "CatalogNumber", all.x = TRUE)
length(hyb_nuc_today$CatalogNumber)


# Plot to see where the dissected fish fall

plot(jitter(hyb_nuc_today_with_metadata$Latitude,30)~jitter(hyb_nuc_today_with_metadata$YearCollected.y,5))+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)
str(hyb_nuc_today_with_metadata)
hyb_nuc_headers<-ls(hyb_nuc_today_with_metadata)
hyb_nuc_headers<-as.factor(hyb_nuc_headers)


# Add parasites across organs and double what needs to be doubled.

ACAN.AD<-hyb_nuc_today_with_metadata$ACAN.AD.INTestine                
CEST.BB<-hyb_nuc_today_with_metadata$CEST.BB.INTESTINE               
CEST.UNK<-hyb_nuc_today_with_metadata$CEST.UNK.INTESTINE
COPE.POOD<-(2*hyb_nuc_today_with_metadata$COPE.POOD.Gill)                   
# not a parasite: CYST.UNK.AnalFin                 
# not a parasite: CYST.UNK.CaudalFin              
# not a parasite: CYST.UNK.CONNECTIVETISSUE        
# not a parasite: CYST.UNK.DorsalFin               
# not a parasite: CYST.UNK.EYE                    
# not a parasite: CYST.UNK.GILL                    
# not a parasite: CYST.UNK.Intestine               
# not a parasite: cyst.unk.kidney                 
# not a parasite: cyst.unk.liver                   
# not a parasite: CYST.UNK.PelvicFin               
# not a parasite: CYST.UNK.SKIN                   
# not a parasite: CYST.UNKF.VOUCH                 
# not a parasite: FIBER.UNK.Intestine             
# not a parasite: FIBER.UNK.STOMACH                
# not a parasite: FIBER.UNK.VOUCH                  
META.UNK<-(2*hyb_nuc_today_with_metadata$META.UNK.EYE)+hyb_nuc_today_with_metadata$META.UNK.Flush+
  hyb_nuc_today_with_metadata$meta.unk.INTESTINE               
MONO.DACT<-(2*hyb_nuc_today_with_metadata$MONO.DACT.Gill)                   
MONO.GDAC<-hyb_nuc_today_with_metadata$MONO.GDAC.PelvicFin+(2*hyb_nuc_today_with_metadata$MONO.UNK.GILL)                    
MYX.AK<-as.numeric(hyb_nuc_today_with_metadata$MYX.AK.GallBladder)              
MYX.EYE<-(2*hyb_nuc_today_with_metadata$MYX.EYE.EYE)+hyb_nuc_today_with_metadata$MYX.EYE.intestine                
MYX.FI<-hyb_nuc_today_with_metadata$MYX.FI.CaudalFin                 
MYX.GC<-(2*hyb_nuc_today_with_metadata$myx.gc.gill)                     
MYX.GL<-(2*hyb_nuc_today_with_metadata$MYX.GL.gill)                     
MYX.OT<-(2*hyb_nuc_today_with_metadata$myx.ot.gill)+(2*hyb_nuc_today_with_metadata$MYX.OT.Skin)                      
MYX.TIN<-(2*hyb_nuc_today_with_metadata$MYX.TIN.GILL)                     

# free-living?
NEM.BUD<-(2*hyb_nuc_today_with_metadata$NEM.BUD.GILL)+(hyb_nuc_today_with_metadata$NEM.BUD.INTESTINE)                
NEM.DAFT<-hyb_nuc_today_with_metadata$NEM.DAFT.CONNECTIVETISSUE       
NEM.UNK<-hyb_nuc_today_with_metadata$NEM.UNK.INTESTINE                
TREM.ASS<-hyb_nuc_today_with_metadata$TREM.ASS.ConnectiveTissue+hyb_nuc_today_with_metadata$TREM.ASS.Flush+
  hyb_nuc_today_with_metadata$TREM.ASS.INTESTINE+as.numeric(hyb_nuc_today_with_metadata$TREM.ASS.LIVER)                 
TREM.BC<-hyb_nuc_today_with_metadata$TREM.BC.CAUDALFIN+hyb_nuc_today_with_metadata$trem.bc.dorsalfin+
  (2*hyb_nuc_today_with_metadata$trem.bc.skin)                     


### Larval trematodes were found in many organs. Sometimes, those were organs that were not found in all fish
### (e.g., gonads). However, even if the gonads (or spleen, or heart) were too small to ID, we feel confident that
### we would have seen metacercariae in those organs, even if the organs themselves were just perceived as part 
### of the visceral mush inside of fish. Therefore, it wouldn't be appropriate to propagate NAs across the entire
### TREM.CM count within a fish (across organs). The code below allows the sum of TREM.CM within an individual
### fish ignore the NAs, and therefore provide a total number of TREM.META regardless of whether one of the organs
### was "missing".

hyb_nuc_today_with_metadata$TREM.CM.EYE.DUB<-2*hyb_nuc_today_with_metadata$TREM.CM.EYE
hyb_nuc_today_with_metadata$TREM.CM.GILL.DUB<-2*hyb_nuc_today_with_metadata$TREM.CM.GILL
hyb_nuc_today_with_metadata$TREM.CM.AnalFin<-as.numeric(hyb_nuc_today_with_metadata$TREM.CM.AnalFin)
hyb_nuc_today_with_metadata$TREM.CM.BodyCavity<-as.numeric(hyb_nuc_today_with_metadata$TREM.CM.BodyCavity)
hyb_nuc_today_with_metadata$TREM.CM.ConnectiveTissue<-as.numeric(hyb_nuc_today_with_metadata$TREM.CM.ConnectiveTissue)
hyb_nuc_today_with_metadata$TREM.CM.Flush<-as.numeric(hyb_nuc_today_with_metadata$TREM.CM.Flush)
hyb_nuc_today_with_metadata$TREM.CM.Gonad<-as.numeric(hyb_nuc_today_with_metadata$TREM.CM.Gonad)
hyb_nuc_today_with_metadata$TREM.CM.INTESTINE<-as.numeric(hyb_nuc_today_with_metadata$TREM.CM.INTESTINE)
hyb_nuc_today_with_metadata$TREM.CM.Kidney<-as.numeric(hyb_nuc_today_with_metadata$TREM.CM.Kidney)
hyb_nuc_today_with_metadata$TREM.CM.Liver<-as.numeric(hyb_nuc_today_with_metadata$TREM.CM.Liver)
hyb_nuc_today_with_metadata$TREM.CM.Spleen<-as.numeric(hyb_nuc_today_with_metadata$TREM.CM.Spleen)
hyb_nuc_today_with_metadata$TREM.CM.Stomach<-as.numeric(hyb_nuc_today_with_metadata$TREM.CM.Stomach)


TREM.CM<-rowSums(hyb_nuc_today_with_metadata[,c("TREM.CM.AnalFin","TREM.CM.BodyCavity",
                                            "TREM.CM.ConnectiveTissue",
                                            "TREM.CM.EYE.DUB","TREM.CM.Flush","TREM.CM.GILL.DUB",
                                            "TREM.CM.Gonad","TREM.CM.INTESTINE","TREM.CM.Kidney",
                                            "TREM.CM.Liver","TREM.CM.Spleen",
                                            "TREM.CM.Stomach")], na.rm=TRUE)

TREM.DIPLO<-(2*hyb_nuc_today_with_metadata$TREM.DIPLO.EYE)                   

# is this different from all the other TREM.METAs?
TREM.META<-hyb_nuc_today_with_metadata$TREM.META.ES.CAUDALFIN          

TREM.META.ES<-(2*hyb_nuc_today_with_metadata$TREM.META.ES.eye)                 
TREM.META.GO<-(2*hyb_nuc_today_with_metadata$TREM.META.GO.Eye)+(2*hyb_nuc_today_with_metadata$TREM.META.GO.GILL)+
  hyb_nuc_today_with_metadata$TREM.META.GO.INTESTINE    


### Larval trematodes were found in many organs. Sometimes, those were organs that were not found in all fish
### (e.g., gonads). However, even if the gonads (or spleen, or heart) were too small to ID, we feel confident that
### we would have seen metacercariae in those organs, even if the organs themselves were just perceived as part 
### of the visceral mush inside of fish. Therefore, it wouldn't be appropriate to propagate NAs across the entire
### TREM.META.HET count within a fish (across organs). The code below allows the sum of TREM.META.HET within an individual
### fish ignore the NAs, and therefore provide a total number of TREM.META regardless of whether one of the organs
### was "missing".

hyb_nuc_today_with_metadata$TREM.META.HET.pectoralfin.DUB<-2*hyb_nuc_today_with_metadata$TREM.META.HET.pectoralfin
hyb_nuc_today_with_metadata$TREM.META.HET.PELVICFIN.DUB<-2*hyb_nuc_today_with_metadata$TREM.META.HET.PELVICFIN
hyb_nuc_today_with_metadata$TREM.META.HET.SKIN.DUB<-2*hyb_nuc_today_with_metadata$TREM.META.HET.SKIN

TREM.META.HET<-rowSums(hyb_nuc_today_with_metadata[,c("TREM.META.HET.pectoralfin.DUB","TREM.META.HET.PELVICFIN.DUB",
                                                "TREM.META.HET.SKIN.DUB",
                                                "TREM.META.HET.ANALFIN","TREM.META.HET.CAUDALFIN")], na.rm=TRUE)

TREM.META.SP<-hyb_nuc_today_with_metadata$TREM.META.SP.AnalFin+hyb_nuc_today_with_metadata$TREM.META.SP.CaudalFin+
  hyb_nuc_today_with_metadata$TREM.META.SP.DorsalFin+(2*hyb_nuc_today_with_metadata$TREM.META.SP.pectoral.fin)+
  (2*hyb_nuc_today_with_metadata$TREM.META.SP.PelvicFin)           
TREM.META.UNK<-(2*hyb_nuc_today_with_metadata$TREM.META.UNK.EYE)                
TREM.NS<-hyb_nuc_today_with_metadata$TREM.NS.INTESTINE               
TREM.UNK<-(2*hyb_nuc_today_with_metadata$TREM.UNK.EYE)+hyb_nuc_today_with_metadata$TREM.UNK.Flush+
  (2*hyb_nuc_today_with_metadata$TREM.UNK.GILL)+as.numeric(hyb_nuc_today_with_metadata$TREM.UNK.Liver)


# Now put it all together

hyb_nuc_processed_data<-cbind.data.frame(hyb_nuc_today_with_metadata$CatalogNumber,hyb_nuc_today_with_metadata$YearCollected.x,
                                         hyb_nuc_today_with_metadata$MonthCollected.x,hyb_nuc_today_with_metadata$DayCollected.x,
                                         hyb_nuc_today_with_metadata$IndividualFishID,hyb_nuc_today_with_metadata$Dissector,
                                         hyb_nuc_today_with_metadata$DissectionDate,hyb_nuc_today_with_metadata$Sex,
                                         hyb_nuc_today_with_metadata$TotalLength_mm,hyb_nuc_today_with_metadata$StandardLength_mm,
                                         hyb_nuc_today_with_metadata$weight_mg,hyb_nuc_today_with_metadata$CI.y,
                                         hyb_nuc_today_with_metadata$combo,
                                         hyb_nuc_today_with_metadata$Latitude,
                                         hyb_nuc_today_with_metadata$Longitude,
                                         ACAN.AD,CEST.BB,CEST.UNK,COPE.POOD,META.UNK,MONO.DACT,MONO.GDAC,
                                         MYX.AK,MYX.EYE,MYX.FI,MYX.GC,MYX.GL,MYX.OT,MYX.TIN,NEM.BUD,NEM.DAFT,
                                         NEM.UNK,TREM.ASS,TREM.BC,TREM.CM,TREM.DIPLO,TREM.META,TREM.META.ES,
                                         TREM.META.GO,TREM.META.HET,TREM.META.SP,TREM.META.UNK,TREM.NS,TREM.UNK)

# Name the columns

colnames(hyb_nuc_processed_data)[1]<-"CatalogNumber"
colnames(hyb_nuc_processed_data)[2]<-"YearCollected"
colnames(hyb_nuc_processed_data)[3]<-"MonthCollected"
colnames(hyb_nuc_processed_data)[4]<-"DayCollected"
colnames(hyb_nuc_processed_data)[5]<-"IndividualFishID"
colnames(hyb_nuc_processed_data)[6]<-"Dissector_and_Examiner"
colnames(hyb_nuc_processed_data)[7]<-"DissectionDate"
colnames(hyb_nuc_processed_data)[8]<-"Sex"
colnames(hyb_nuc_processed_data)[9]<-"TotalLength_mm"
colnames(hyb_nuc_processed_data)[10]<-"StandardLength_mm"
colnames(hyb_nuc_processed_data)[11]<-"Weight_mg"
colnames(hyb_nuc_processed_data)[12]<-"CI"
colnames(hyb_nuc_processed_data)[13]<-"combo"
colnames(hyb_nuc_processed_data)[14]<-"Latitude"
colnames(hyb_nuc_processed_data)[15]<-"Longitude"

View(hyb_nuc_processed_data)

#Remove #VALUE! entries from Weight_mg with "NA"
#not_ath_processed_data <- not_ath_processed_data %>%
#  mutate(Weight_mg = recode(Weight_mg,
#                            "#VALUE!" = "NA"))

#Evaluate relationship between TotalLength_mm and Weight_mg. There is evidently a mistake. Some values must be multiplied by 10.
# ggplot(hyb_nuc_processed_data,aes(TotalLength_mm,Weight_mg))+geom_point(size=4)

# Evaluate ratio between total length and weight to inform decision about which numbers must be multiplied by ten. For this, divide length by weight
# not_ath_processed_data$LenWeig <- as.numeric(not_ath_processed_data$Weight_mg)/not_ath_processed_data$TotalLength_mm


# All incorrect weights that must be multiplied by 10 have a ratio below 3 which does not make sense. See the data below to confirm. The weights have a decimal point which should not be the case.
# lenweighb <- subset(not_ath_processed_data, LenWeig<3) 
# View(lenweighb)

# Therefore, multiply all weights with LenWeig higher than 3 by 10

#not_ath_processed_data$Weight_mg <- ifelse(not_ath_processed_data$LenWeig < 3,                # condition
#                                           as.numeric(not_ath_processed_data$Weight_mg)*10,    # what if condition is TRUE
#                                           as.numeric(not_ath_processed_data$Weight_mg)       # what if condition is FALSE
#)


#Re-evaluate relationship between TotalLength_mm and Weight_mg. Now it is fixed.
#ggplot(not_ath_processed_data,aes(TotalLength_mm,Weight_mg))+geom_point(size=4)


# It looks like meta-data are missing for some of the fish.

stuff<-hyb_nuc_processed_data %>%
  filter(is.na(CI))

stuff<-hyb_nuc_processed_data %>%
  filter(CatalogNumber=="17623")

hyb_nuc_processed_data$CatalogNumber[hyb_nuc_processed_data$CatalogNumber=="17623"]<-"175623"

stuff<-hyb_nuc_processed_data %>%
  filter(CatalogNumber=="175623")

hyb_nuc_processed_data$CI[hyb_nuc_processed_data$IndividualFishID=="175623_04"]<-"impact"
hyb_nuc_processed_data$combo[hyb_nuc_processed_data$IndividualFishID=="175623_04"]<-"impact_1994-2003"
hyb_nuc_processed_data$Latitude[hyb_nuc_processed_data$IndividualFishID=="175623_04"]<-"30.75417"
hyb_nuc_processed_data$Longitude[hyb_nuc_processed_data$IndividualFishID=="175623_04"]<-"-89.82694"

stuff<-hyb_nuc_processed_data %>%
  filter(CatalogNumber=="175623")

stuff<-hyb_nuc_processed_data %>%
  filter(is.na(CI))

hyb_nuc_processed_data$CatalogNumber[hyb_nuc_processed_data$IndividualFishID=="156660_01"]<-"156660"

stuff<-hyb_nuc_processed_data %>%
  filter(CatalogNumber=="156660")

hyb_nuc_processed_data$CI[hyb_nuc_processed_data$IndividualFishID=="156660_01"]<-"control"
hyb_nuc_processed_data$combo[hyb_nuc_processed_data$IndividualFishID=="156660_01"]<-"control_1984-1993"
hyb_nuc_processed_data$Latitude[hyb_nuc_processed_data$IndividualFishID=="156660_01"]<-"30.76528"
hyb_nuc_processed_data$Longitude[hyb_nuc_processed_data$IndividualFishID=="156660_01"]<-"-89.83222"

stuff<-hyb_nuc_processed_data %>%
  filter(CatalogNumber=="156660")

stuff<-hyb_nuc_processed_data %>%
  filter(is.na(CI))

hyb_nuc_processed_data$CatalogNumber[hyb_nuc_processed_data$IndividualFishID=="106828_03"]<-"106828"

stuff<-hyb_nuc_processed_data %>%
  filter(CatalogNumber=="106828")

hyb_nuc_processed_data$CI[hyb_nuc_processed_data$IndividualFishID=="106828_03"]<-"control"
hyb_nuc_processed_data$combo[hyb_nuc_processed_data$IndividualFishID=="106828_03"]<-"control_1974-1983"
hyb_nuc_processed_data$Latitude[hyb_nuc_processed_data$IndividualFishID=="106828_03"]<-"30.76805"
hyb_nuc_processed_data$Longitude[hyb_nuc_processed_data$IndividualFishID=="106828_03"]<-"-89.83083"

stuff<-hyb_nuc_processed_data %>%
  filter(CatalogNumber=="106828")

stuff<-hyb_nuc_processed_data %>%
  filter(is.na(CI))

hyb_nuc_processed_data$CatalogNumber[hyb_nuc_processed_data$CatalogNumber=="184748"]<-"184784"

stuff<-hyb_nuc_processed_data %>%
  filter(CatalogNumber=="184784")

hyb_nuc_processed_data$CI[hyb_nuc_processed_data$CatalogNumber=="184784"]<-"impact"
hyb_nuc_processed_data$combo[hyb_nuc_processed_data$CatalogNumber=="184784"]<-"impact_1994-2003"
hyb_nuc_processed_data$Latitude[hyb_nuc_processed_data$CatalogNumber=="184784"]<-"30.74195"
hyb_nuc_processed_data$Longitude[hyb_nuc_processed_data$CatalogNumber=="184784"]<-"-89.82528"

stuff<-hyb_nuc_processed_data %>%
  filter(CatalogNumber=="184784")

stuff<-hyb_nuc_processed_data %>%
  filter(is.na(CI))



# Make the dataset analyzable

hyb_nuc_processed_data_longer<-melt(hyb_nuc_processed_data,id=c("CatalogNumber", "YearCollected", 
                                                                "MonthCollected", "DayCollected", 
                                                                "IndividualFishID", "Dissector_and_Examiner",
                                                                "DissectionDate", "Sex", 
                                                                "TotalLength_mm", "StandardLength_mm",
                                                                "Weight_mg","CI","combo","Latitude","Longitude"))

colnames(hyb_nuc_processed_data_longer)[16]<-"psite_spp"
colnames(hyb_nuc_processed_data_longer)[17]<-"psite_count"


# Export both sheets

write.csv(hyb_nuc_processed_data_longer, file="data/processed/Hybognathus_nuchalis_processed_machine_readable.csv")
write.csv(hyb_nuc_processed_data, file="data/processed/Hybognathus_nuchalis_processed_human_readable.csv")




### PERVIG: Percina vigil----


# You can also do it the old-fashioned way

per_vig_today<-read.csv("data/raw/Percina_vigil_Datasheet - Sheet1_2024.07.29.csv")
length(per_vig_today$CatalogNumber)


# Now merge dissection data with meta-data (all.x makes sure that you keep all individual fish from the same lot, even though
# they have the same catalog number).

per_vig_today_with_metadata<-merge(per_vig_today, meta_data, by.x = "CatalogNumber", by.y = "CatalogNumber", all.x = TRUE)
length(per_vig_today_with_metadata$CatalogNumber)


# Plot to see where the dissected fish fall

plot(jitter(per_vig_today_with_metadata$Latitude.y,30)~jitter(per_vig_today_with_metadata$YearCollected.y,5))+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)
str(per_vig_today_with_metadata)
per_vig_headers<-ls(per_vig_today_with_metadata)
per_vig_headers<-as.factor(per_vig_headers)


# Add parasites across organs and double what needs to be doubled.

ACANTH.BB<-per_vig_today_with_metadata$ACANTH.BB.FLUSH+per_vig_today_with_metadata$ACANTH.BB.INTESTINE    
ACANTH.CYA<-per_vig_today_with_metadata$ACANTH.CYA.CONNECTIVETISSUE+
  per_vig_today_with_metadata$ACANTH.CYA.FLUSH+per_vig_today_with_metadata$acanth.cya.intestine            
ACANTH.THIN<-per_vig_today_with_metadata$ACANTH.THIN.CONNECTIVE           
CEST.BOT<-per_vig_today_with_metadata$CEST.BOT.INTESTINE+per_vig_today_with_metadata$CEST.BOT.Stomach                
# not known whether this is a parasite: CYST.LOBE.ConnectiveTissue       
# not known whether this is a parasite: CYST.LOBE.FLUSH                 
# not known whether this is a parasite: CYST.LOBE.INTESTINE              
# not known whether this is a parasite: CYST.LOBE.Liver                  
# this is cancer, not a parasite: CYST.PRO.AnalFin                 
# this is cancer, not a parasite: CYST.PRO.CaudalFin               
# this is cancer, not a parasite: CYST.PRO.DorsalFin              
# this is cancer, not a parasite: CYST.PRO.PectoralFin             
# this is cancer, not a parasite: CYST.PRO.PelvicFin               
# not known whether this is a parasite: CYST.UNK.CONNECTIVETISSUES       
# not known whether this is a parasite: CYST.UNK.EYE                     
# not known whether this is a parasite: CYST.UNK.GILL                   
# not known whether this is a parasite: CYST.UNKN.INTESTINE              
MONO.GYRO<-(2*per_vig_today_with_metadata$MONO.GYRO.Gill)                  
NEM.HEATH<-per_vig_today_with_metadata$NEM.HEATH.INTESTINE              
NEM.LARV<-per_vig_today_with_metadata$NEM.LARV.CONNECTIVETISSUE+(2*per_vig_today_with_metadata$nem.larv.eye)+
  per_vig_today_with_metadata$NEM.LARV.Flush+per_vig_today_with_metadata$NEM.LARV.INTESTINE+
  per_vig_today_with_metadata$NEM.LARV.Liver                   
NEM.SQUIG<-per_vig_today_with_metadata$NEM.SQUIG.INTESTINE              
NEM.UNK<-per_vig_today_with_metadata$nem.unk.connectivetissue+per_vig_today_with_metadata$NEM.UNK.INTESTINE                
TREM.CHONK<-per_vig_today_with_metadata$TREM.CHONK.BODYCAV              
TREM.CL<-per_vig_today_with_metadata$TREM.CL.CONNECTIVETISSUES        
TREM.ENDOUS<-per_vig_today_with_metadata$TREM.ENDOUS.CONNECTIVE+per_vig_today_with_metadata$TREM.ENDOUS.LIVER                
TREM.ISPOT<-(2*per_vig_today_with_metadata$TREM.ISPOT.GILL)                  

per_vig_today_with_metadata$TREM.PHLb.UrinaryBladder 

# We agreed that TREM.PHLB would be 0 even where urinary bladder wasn't specifically checked, since these trems
# are big and difficult to miss. Write a loop to fix.

per_vig_today_with_metadata$TREM.PHLB<-vector("numeric",length(per_vig_today_with_metadata$CatalogNumber))

for(i in 1:length(per_vig_today_with_metadata$CatalogNumber)) {
  
  if(is.na(per_vig_today_with_metadata$TREM.PHLb.UrinaryBladder[i])){
    per_vig_today_with_metadata$TREM.PHLB[i] <- 0
    
  } else {
    
    if(per_vig_today_with_metadata$TREM.PHLb.UrinaryBladder[i] > 0){
      per_vig_today_with_metadata$TREM.PHLB[i] <- per_vig_today_with_metadata$TREM.PHLb.UrinaryBladder[i]
      
      } else {
          
        per_vig_today_with_metadata$TREM.PHLB[i] <- 0
      }
  }
}

per_vig_today_with_metadata$TREM.PHLB<-as.numeric(per_vig_today_with_metadata$TREM.PHLB) 

per_vig_today_with_metadata$TREM.PHLB[is.na(per_vig_today_with_metadata$TREM.PHLB)]<-0
TREM.PHLB<-per_vig_today_with_metadata$TREM.PHLB

TREM.SK<-per_vig_today_with_metadata$TREM.SK.BodyCavity+(2*per_vig_today_with_metadata$TREM.SK.Eye)+
  (2*per_vig_today_with_metadata$TREM.Sk.SKIN)                    



# Now put it all together

per_vig_processed_data<-cbind.data.frame(per_vig_today_with_metadata$CatalogNumber,per_vig_today_with_metadata$YearCollected.x,
                                         per_vig_today_with_metadata$MonthCollected.x,per_vig_today_with_metadata$DayCollected.x,
                                         per_vig_today_with_metadata$IndividualFishID,per_vig_today_with_metadata$Dissector,
                                         per_vig_today_with_metadata$DissectionDate,per_vig_today_with_metadata$Sex,
                                         per_vig_today_with_metadata$TotalLength_mm,per_vig_today_with_metadata$StandardLength_mm,
                                         per_vig_today_with_metadata$weight_mg,per_vig_today_with_metadata$CI,
                                         per_vig_today_with_metadata$combo,
                                         per_vig_today_with_metadata$Latitude.y,
                                         per_vig_today_with_metadata$Longitude.y,
                                         ACANTH.BB,ACANTH.CYA,ACANTH.THIN,CEST.BOT,MONO.GYRO,NEM.HEATH,
                                         NEM.LARV,NEM.SQUIG,NEM.UNK,TREM.CHONK,TREM.CL,TREM.ENDOUS,
                                         TREM.ISPOT,TREM.PHLB,TREM.SK)

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

View(per_vig_processed_data)


# Make the dataset analyzable
library(reshape2)
per_vig_processed_data_longer<-melt(per_vig_processed_data,id=c("CatalogNumber", "YearCollected", 
                                                                "MonthCollected", "DayCollected", 
                                                                "IndividualFishID", "Dissector_and_Examiner",
                                                                "DissectionDate", "Sex", 
                                                                "TotalLength_mm", "StandardLength_mm",
                                                                "Weight_mg","CI","combo","Latitude","Longitude"))

colnames(per_vig_processed_data_longer)[16]<-"psite_spp"
colnames(per_vig_processed_data_longer)[17]<-"psite_count"


# Export both sheets

write.csv(per_vig_processed_data_longer, file="data/processed/Percina_vigil_processed_machine_readable.csv")
write.csv(per_vig_processed_data, file="data/processed/Percina_vigil_processed_human_readable.csv")




### Put all the sheets together----

# Start by scaling fish body size within fish species

pim_vig_processed_data_longer$scaled_TL<-scale(pim_vig_processed_data_longer$TotalLength_mm)
ict_pun_processed_data_longer$scaled_TL<-scale(ict_pun_processed_data_longer$TotalLength_mm)
not_ath_processed_data_longer$scaled_TL<-scale(not_ath_processed_data_longer$TotalLength_mm)
hyb_nuc_processed_data_longer$scaled_TL<-scale(hyb_nuc_processed_data_longer$TotalLength_mm)
per_vig_processed_data_longer$scaled_TL<-scale(per_vig_processed_data_longer$TotalLength_mm)


# Then make sure that there is a column for the fish species

pim_vig_processed_data_longer$Fish_sp<-c(rep("Pimephales vigilax",length(pim_vig_processed_data_longer$CatalogNumber)))
ict_pun_processed_data_longer$Fish_sp<-c(rep("Ictalurus punctatus",length(ict_pun_processed_data_longer$CatalogNumber)))
not_ath_processed_data_longer$Fish_sp<-c(rep("Notropis atherinoides",length(not_ath_processed_data_longer$CatalogNumber)))
hyb_nuc_processed_data_longer$Fish_sp<-c(rep("Hybognathus nuchalis",length(hyb_nuc_processed_data_longer$CatalogNumber)))
per_vig_processed_data_longer$Fish_sp<-c(rep("Percina vigil",length(per_vig_processed_data_longer$CatalogNumber)))

full_dataset<-rbind.data.frame(pim_vig_processed_data_longer,ict_pun_processed_data_longer,
                                          not_ath_processed_data_longer,hyb_nuc_processed_data_longer,
                               per_vig_processed_data_longer)

full_dataset$fish_psite_combo<-paste(full_dataset$Fish_sp,full_dataset$psite_spp,sep="_")

life_histories<-read.csv("data/raw/Parasite_Life_History_Strategies_2024.07.29.csv",header=T,sep=",")
life_histories$fish_psite_combo<-paste(life_histories$Fish_sp,life_histories$psite_spp,sep="_")


full_dataset_with_LH<-merge(full_dataset, life_histories, by.x = "fish_psite_combo", by.y = "fish_psite_combo", all.x = TRUE)

View(full_dataset_with_LH)

colnames(full_dataset_with_LH)[19]<-"scaled_TL_mm"
full_dataset_with_LH$YearCollected<-as.numeric(full_dataset_with_LH$YearCollected)


# Write a loop to create a before/after column

full_dataset_with_LH$before_after<-vector("character",length(full_dataset_with_LH$fish_psite_combo))

for(i in 1:length(full_dataset_with_LH$fish_psite_combo)) {
  
    if(full_dataset_with_LH$YearCollected[i] > 1972){
      full_dataset_with_LH$before_after[i] <- "after"
      
    } else {
      
      if(full_dataset_with_LH$YearCollected[i] < 1973){
        full_dataset_with_LH$before_after[i] <- "before"
        
      } else {
                    
        full_dataset_with_LH$before_after[i] <- "NA"
      }
    }
}
                 
full_dataset_with_LH$before_after     
              


# Export the sheet

write.csv(full_dataset_with_LH, file="data/processed/Full_dataset_with_psite_life_history_info_2024.07.29.csv")


# Playing around with the big analysis

model_draft<-glmer.nb(psite_count~CI*before_after+(1|psite_spp.x),data=full_dataset_with_LH,family="nbinom")
summary(model_draft)

offset(log(scaled_TL_mm))

