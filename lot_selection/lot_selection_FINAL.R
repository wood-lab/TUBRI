# Selecting specimens for TUBRI project
# Chelsea Wood
# Updated May 2024
#we luv worms

library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales) #for transparency
library(tidyverse)
library(ggmap)


# Our goal is to obtain a sample of fish with replication above and below pulp mills before and after
# 1973, when the Clean Water Act was implemented.

# "before mitigation" = 1963-1972
# "after mitigation" = 1973-2005

# fish host species = aiming for 10. Selecting from 16 species identified as among the most 
# common in the TUBRI collection. In order of declining trophic level, they are:
# Channel catfish Ictalurus punctatus 
# Spotted gar Lepisosteus oculatus
# Largemouth black bass Micropterus salmoides
# Longear sunfish Lepomis megalotis
# Crystal darter Crystallaria asprella
# Pearl darter Percina aurora
# Frecklebelly madtom Noturus munitus
# Saddleback darter Percina vigil
# Blacktail shiner Cyprinella venusta
# Gulf chub Macrorhybopsis tomellerii
# Bluegill Lepomis macrochirus
# Mosquitofish Gambusia affinis
# Bullhead minnow Pimephales vigilax
# Emerald shiner Notropis atherinoides
# Alabama shad Alosa alabamae
# Highfin carpsucker Carpiodes velifer
# Mississippi silvery minnow Hybognathus nuchalis

# total replication = aim to dissect 4,000 fish in total: 20 individual fish of each species 
# (n = 10 species) from each impact category (i.e., control versus polluted) at each mill 
# (n = 2 mills in total) in each decade (n = 5 decades)

# The two sites are MS_lower (Pearl River) and AL (Alabama River)

# Preliminary investigation revealed extremely low sampling above the mill on the Alabama River.
# We must proceed with just one mill - MS_lower.


# First step is to load the full dataset. This was derived by searching for everything within a geographic
# polygon specific to the MS_lower site via FishNet2.

lots_MS_lower<-read.csv("lot_selection/MS_lower_search_results.csv")


#Filter out just the ten species you want to focus on

lots_MS_lower <- lots_MS_lower %>%
  filter(ScientificName == "Pimephales vigilax" | 
           ScientificName == "Notropis atherinoides" | 
           ScientificName == "Gambusia affinis" |
           ScientificName == "Ictalurus punctatus" | 
           ScientificName == "Lepomis macrochirus" | 
           ScientificName == "Carpiodes velifer" |
           ScientificName == "Hybognathus nuchalis" | 
           ScientificName == "Micropterus salmoides" |
           ScientificName == "Lepomis megalotis" | 
           ScientificName == "Percina vigil" | 
           ScientificName == "Lepisosteus oculatus" | 
           ScientificName == "Crystallaria asprella" | 
           ScientificName == "Cyprinella venusta" | 
           ScientificName == "Macrhybopsis tomellerii" | 
           ScientificName == "Percina aurora"| 
           ScientificName == "Noturus munitus"| 
           ScientificName == "Alosa alabamae")

# Create a new column to indicate whether you're dealing with control or impact

lots_MS_lower$CI<-vector("character",length(lots_MS_lower$InstitutionCode))

for(i in 1:length(lots_MS_lower$InstitutionCode)) {
  
  if(is.na(lots_MS_lower$Latitude[i])){
    lots_MS_lower$CI[i] <- NA
  } else {
    
    if(lots_MS_lower$Latitude[i] > 30.76){
      lots_MS_lower$CI[i] <- "control"
      
    } else {
      
      lots_MS_lower$CI[i] <- "impact"
      
    }
  }
}
lots_MS_lower$CI


# Just make sure everything looks okay and is making sense

install.packages("devtools")
devtools::install_github("stadiamaps/ggmap") 
library("ggmap")
register_stadiamaps("45588e55-a703-4f6e-9ef0-ee6ed672b73a", write = TRUE)

bounds<-c(left=-89.87, bottom=30.67, right=-89.79, top=30.81)
map_check<-get_stadiamap(bounds, zoom=11, maptype = "stamen_terrain_background") %>% ggmap()+
  geom_point(data = lots_MS_lower, aes(x=Longitude,y=Latitude,fill=CI),shape=21)+
  scale_fill_manual(values=c("white","#bdbdbd"))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))
map_check


# Some weird points not located on the river ~30.72 latitude. Get rid of them.

lots_MS_lower<-lots_MS_lower %>%
  filter(Longitude > -89.85)

map_check<-get_stadiamap(bounds, zoom=11, maptype = "stamen_terrain_background") %>% ggmap()+
  geom_point(data = lots_MS_lower, aes(x=Longitude,y=Latitude,fill=CI),shape=21)+
  scale_fill_manual(values=c("white","#bdbdbd"))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=8),axis.text.y=element_text(size=8))
map_check


# How many of each species is there in this stretch of the Pearl River?

total_counts<-lots_MS_lower %>%
  group_by(ScientificName) %>%
  summarize(total_available = sum(IndividualCount,na.rm = T))
total_counts


# Some of these clearly are not going to work out, including 
# A. alabamae, C. asprella, L. oculatus, P. aurora, and M. tomellerii (which doesn't even appear in this
# stretch of the Pearl River).


# Create a separate dataset for each of the remaining fish species.

car_vel<-lots_MS_lower %>%
  filter(ScientificName =="Carpiodes velifer")

cyp_ven<-lots_MS_lower %>%
  filter(ScientificName =="Cyprinella venusta")

gam_aff<-lots_MS_lower %>%
  filter(ScientificName =="Gambusia affinis")

hyb_nuc<-lots_MS_lower %>%
  filter(ScientificName =="Hybognathus nuchalis")

ict_pun<-lots_MS_lower %>%
  filter(ScientificName =="Ictalurus punctatus")

lep_mac<-lots_MS_lower %>%
  filter(ScientificName =="Lepomis macrochirus")

lep_meg<-lots_MS_lower %>%
  filter(ScientificName =="Lepomis megalotis")

mic_sal<-lots_MS_lower %>%
  filter(ScientificName =="Micropterus salmoides")

not_ath<-lots_MS_lower %>%
  filter(ScientificName =="Notropis atherinoides")

not_mun<-lots_MS_lower %>%
  filter(ScientificName =="Noturus munitus")

per_vig<-lots_MS_lower %>%
  filter(ScientificName =="Percina vigil")

pim_vig<-lots_MS_lower %>%
  filter(ScientificName =="Pimephales vigilax")


# Make sure everything looks good and exclude species / latitudes where you don't have before and after data for
# every impact category or site.

# CAR_VEL

plot(jitter(car_vel$Latitude,10)~car_vel$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

#Need to cut off high and low latitudes, which don't extend across the entire time series

car_vel<-car_vel %>%
  filter(Latitude < 30.78 & Latitude > 30.72)

plot(jitter(car_vel$Latitude,10)~car_vel$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# CYP_VEN - EXCLUDE - Only collected very recently

plot(jitter(cyp_ven$Latitude,10)~cyp_ven$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# GAM_AFF - This species doesn't have enough pre-1973 data from above the mill. However, it IS an invasive
# species, so we should keep it for comparisons of native/non-native species

plot(jitter(gam_aff$Latitude,10)~gam_aff$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

#Need to cut off high latitudes, which don't extend across the entire time series

gam_aff<-gam_aff %>%
  filter(Latitude < 30.79)

plot(jitter(gam_aff$Latitude,10)~gam_aff$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# HYB_NUC

plot(jitter(hyb_nuc$Latitude,10)~hyb_nuc$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

#This looks great, just need to trim that one stray site in the south and the most northerly sites

hyb_nuc<-hyb_nuc %>%
  filter(Latitude > 30.7 & Latitude < 30.78)

plot(jitter(hyb_nuc$Latitude,10)~hyb_nuc$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# ICT_PUN

plot(jitter(ict_pun$Latitude,10)~ict_pun$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

#Need to trim the most northerly and southerly sites for this species

ict_pun<-ict_pun %>%
  filter(Latitude < 30.78 & Latitude > 30.7)

plot(jitter(ict_pun$Latitude,10)~ict_pun$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# LEP_MAC - EXCLUDE THIS SPECIES - Very little pre-1973 data from above the mill

plot(jitter(lep_mac$Latitude,10)~lep_mac$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# LEP_MEG - EXCLUDE THIS SPECIES - Very little pre-1973 data from above the mill

plot(jitter(lep_meg$Latitude,10)~lep_meg$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# MIC_SAL -EXCLUDE THIS SPECIES - Very little pre-1973 data from above the mill

plot(jitter(mic_sal$Latitude,10)~mic_sal$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# NOT_ATH

plot(jitter(not_ath$Latitude,10)~not_ath$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

#Just trim the northernmost site, which has very little pre-1973 replication

not_ath<-not_ath %>%
  filter(Latitude < 30.79)

plot(jitter(not_ath$Latitude,10)~not_ath$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# NOT_MUN - EXCLUDE - Limited to just one site downstream of the mill

plot(jitter(not_mun$Latitude,10)~not_mun$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# PER_VIG - BORDERLINE - Very little data from above the mill

plot(jitter(per_vig$Latitude,10)~per_vig$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# Trim northerly site

per_vig<-per_vig %>%
  filter(Latitude < 30.79)

plot(jitter(per_vig$Latitude,10)~per_vig$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# PIM_VIG

plot(jitter(pim_vig$Latitude,10)~pim_vig$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

# Trim northerly site

pim_vig<-pim_vig %>%
  filter(Latitude < 30.79)

plot(jitter(pim_vig$Latitude,10)~pim_vig$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)


#Okay, now we've trimmed each species to the appropriate latitudes and cut out some species. Re-create the
#dataset.

lots_MS_lower<-rbind.data.frame(car_vel,gam_aff,hyb_nuc,ict_pun,not_ath,per_vig,pim_vig)



# Create a new column to indicate what decade you're in

lots_MS_lower$decade<-vector("character",length(lots_MS_lower$InstitutionCode))


# Now loop to populate the decade column for the first mill

for(i in 1:length(lots_MS_lower$InstitutionCode)) {
  
  if(is.na(lots_MS_lower$YearCollected[i])){
    lots_MS_lower$decade[i] <- NA
    
  } else {
    
    if(lots_MS_lower$YearCollected[i]>2013){
      lots_MS_lower$decade[i] <- NA
      
    } else {
      
      if(lots_MS_lower$YearCollected[i] <1964){
        lots_MS_lower$decade[i] <- "1954-1963"
        
      } else {
        
        if(lots_MS_lower$YearCollected[i] <1974 & lots_MS_lower$YearCollected[i] > 1963){
          lots_MS_lower$decade[i] <- "1964-1973"
          
        } else {
          
          if(lots_MS_lower$YearCollected[i] <1984 & lots_MS_lower$YearCollected[i] > 1973){
            lots_MS_lower$decade[i] <- "1974-1983"
            
          } else {
            
            if(lots_MS_lower$YearCollected[i] <1994 & lots_MS_lower$YearCollected[i] > 1983){
              lots_MS_lower$decade[i] <- "1984-1993"
              
            } else {
              
              if(lots_MS_lower$YearCollected[i] <2004 & lots_MS_lower$YearCollected[i] > 1993){
                lots_MS_lower$decade[i] <- "1994-2003"
                
              } else {
                
                if(lots_MS_lower$YearCollected[i] <2014 & lots_MS_lower$YearCollected[i] > 2003){
                  lots_MS_lower$decade[i] <- "2004-2013"
                  
                }}}}}}}}}

lots_MS_lower$decade


# GAM_AFF
gam_aff<-lots_MS_lower %>%
  filter(ScientificName =="Gambusia affinis")

## there are so few of these in control-before, drop down the threshold for picking a jar from 8 to 4

gam_aff<-gam_aff[which(gam_aff$IndividualCount >= 4),]

gam_aff_matrix<-gam_aff %>%
  group_by(CI, decade) %>%
  summarize(total_available = n())

gam_aff$combo<-paste(gam_aff$CI,gam_aff$decade,sep="_")

#control_1954to1963<-gam_aff[ sample( which( gam_aff$combo == "control_1954-1963"), 3, replace = F), ]#not enough for 5
control_1964to1973<-gam_aff[ sample( which( gam_aff$combo == "control_1964-1973"), 2, replace = F), ]#not enough for 5
control_1974to1983<-gam_aff[ sample( which( gam_aff$combo == "control_1974-1983"), 5, replace = F), ]
control_1984to1993<-gam_aff[ sample( which( gam_aff$combo == "control_1984-1993"), 5, replace = F), ]
control_1994to2003<-gam_aff[ sample( which( gam_aff$combo == "control_1994-2003"), 5, replace = F), ]#not enough for 5
control_2004to2013<-gam_aff[ sample( which( gam_aff$combo == "control_2004-2013"), 3, replace = F), ]#not enough for 5
impact_1954to1963<-gam_aff[ which( gam_aff$combo == "impact_1954-1963"), ]#not enough for 5
impact_1964to1973<-gam_aff[ sample( which( gam_aff$combo == "impact_1964-1973"), 5, replace = F), ]
impact_1974to1983<-gam_aff[ sample( which( gam_aff$combo == "impact_1974-1983"), 5, replace = F), ]
impact_1984to1993<-gam_aff[ sample( which( gam_aff$combo == "impact_1984-1993"), 5, replace = F), ]
impact_1994to2003<-gam_aff[ sample( which( gam_aff$combo == "impact_1994-2003"), 5, replace = F), ]
impact_2004to2013<-gam_aff[ sample( which( gam_aff$combo == "impact_2004-2013"), 5, replace = F), ]#don't bother because we don't have this decade for control


gam_aff_selected<-rbind.data.frame(control_1964to1973,control_1974to1983,control_1984to1993,
                                   control_1994to2003,control_2004to2013,
                                   impact_1954to1963,impact_1964to1973,impact_1974to1983,impact_1984to1993,
                                   impact_1994to2003,impact_2004to2013)

gam_aff_matrix<-gam_aff_selected %>%
  group_by(combo) %>%
  summarize(total_request = n())

plot(jitter(gam_aff_selected$Latitude,5)~gam_aff_selected$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)


#HYB_NUC
hyb_nuc<-lots_MS_lower %>%
  filter(ScientificName =="Hybognathus nuchalis")

hyb_nuc$combo<-paste(hyb_nuc$CI,hyb_nuc$decade,sep="_")
hyb_nuc_next<-hyb_nuc[which(hyb_nuc$IndividualCount >= 8),]

hyb_nuc_matrix<-hyb_nuc_next %>%
  group_by(CI, decade) %>%
  summarize(total_request = n())

control_1954to1963<-hyb_nuc_next[ sample( which( hyb_nuc_next$combo == "control_1954-1963"), 4, replace = F), ]#not enough for 5
control_1964to1973<-hyb_nuc_next[ sample( which( hyb_nuc_next$combo == "control_1964-1973"), 5, replace = F), ]
control_1974to1983<-hyb_nuc_next[ sample( which( hyb_nuc_next$combo == "control_1974-1983"), 5, replace = F), ]
control_1984to1993<-hyb_nuc_next[ sample( which( hyb_nuc_next$combo == "control_1984-1993"), 5, replace = F), ]
control_1994to2003<-hyb_nuc_next[ sample( which( hyb_nuc_next$combo == "control_1994-2003"), 5, replace = F), ]
control_2004to2013<-hyb_nuc_next[ sample( which( hyb_nuc_next$combo == "control_2004-2013"), 5, replace = F), ]
impact_1954to1963<-hyb_nuc_next[ sample( which( hyb_nuc_next$combo == "impact_1954-1963"), 3, replace = F), ]#not enough for 5
impact_1964to1973<-hyb_nuc_next[ sample( which( hyb_nuc_next$combo == "impact_1964-1973"), 5, replace = F), ]
impact_1974to1983<-hyb_nuc_next[ sample( which( hyb_nuc_next$combo == "impact_1974-1983"), 5, replace = F), ]
impact_1984to1993<-hyb_nuc_next[ sample( which( hyb_nuc_next$combo == "impact_1984-1993"), 5, replace = F), ]
impact_1994to2003<-hyb_nuc_next[ sample( which( hyb_nuc_next$combo == "impact_1994-2003"), 5, replace = F), ]
impact_2004to2013<-hyb_nuc_next[ sample( which( hyb_nuc_next$combo == "impact_2004-2013"), 5, replace = F), ]

hyb_nuc_selected<-rbind.data.frame(control_1954to1963,control_1964to1973,control_1974to1983,control_1984to1993,
                                   control_1994to2003,control_2004to2013,
                                   impact_1954to1963,impact_1964to1973,impact_1974to1983,impact_1984to1993,
                                   impact_1994to2003,impact_2004to2013)

hyb_nuc_matrix<-hyb_nuc_selected %>%
  group_by(combo) %>%
  summarize(total_request = n())

plot(jitter(hyb_nuc_selected$Latitude,5)~hyb_nuc_selected$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)


#ICT_PUN
ict_pun<-lots_MS_lower %>%
  filter(ScientificName =="Ictalurus punctatus")
ict_pun$combo<-paste(ict_pun$CI,ict_pun$decade,sep="_")

# This is a big fish and there are fewer specimens per lot than others, so lower the threshold
ict_pun_next<-ict_pun[which(ict_pun$IndividualCount >= 4),]

ict_pun_matrix<-ict_pun_next %>%
  group_by(CI, decade) %>%
  summarize(total_request = n())

control_1954to1963<-ict_pun_next[ sample( which( ict_pun_next$combo == "control_1954-1963"), 2, replace = F), ]#not enough for 5
control_1964to1973<-ict_pun_next[ sample( which( ict_pun_next$combo == "control_1964-1973"), 5, replace = F), ]
control_1974to1983<-ict_pun_next[ sample( which( ict_pun_next$combo == "control_1974-1983"), 5, replace = F), ]
control_1984to1993<-ict_pun_next[ sample( which( ict_pun_next$combo == "control_1984-1993"), 3, replace = F), ]#not enough for 5
control_1994to2003<-ict_pun_next[ which( ict_pun_next$combo == "control_1994-2003"), ]#not enough for 5
impact_1954to1963<-ict_pun_next[ sample( which( ict_pun_next$combo == "impact_1954-1963"), 5, replace = F), ]
impact_1964to1973<-ict_pun_next[ sample( which( ict_pun_next$combo == "impact_1964-1973"), 5, replace = F), ]
impact_1974to1983<-ict_pun_next[ sample( which( ict_pun_next$combo == "impact_1974-1983"), 5, replace = F), ]
impact_1984to1993<-ict_pun_next[ sample( which( ict_pun_next$combo == "impact_1984-1993"), 5, replace = F), ]
impact_1994to2003<-ict_pun_next[ sample( which( ict_pun_next$combo == "impact_1994-2003"), 4, replace = F), ]#not enough for 5
impact_2004to2013<-ict_pun_next[ which( ict_pun_next$combo == "impact_2004-2013"), ]

ict_pun_selected<-rbind.data.frame(control_1954to1963,control_1964to1973,control_1974to1983,control_1984to1993,
                                   control_1994to2003,
                                   impact_1954to1963,impact_1964to1973,impact_1974to1983,impact_1984to1993,
                                   impact_1994to2003,impact_2004to2013)

ict_pun_matrix<-ict_pun_selected %>%
  group_by(combo) %>%
  summarize(total_request = n())

plot(jitter(ict_pun_selected$Latitude,5)~ict_pun_selected$YearCollected)+
  abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)


# For ICTPUN specifically: most of these fish are fine. For a few, Katie was able to go into the collection and figure out 
# whether the lots were present and contained fish of an appropriate size. For those, we needed to find replacements 
# in a few instances.

# Katie needs replacements for some ICTPUN specimens: 34760, 37883, 115142, 147855, 147873

# 34760 not included

ict_pun_34760 <- ict_pun_selected %>%
  filter(ScientificName == "Ictalurus punctatus" & CatalogNumber == "34760")

# 37883 not included

ict_pun_37883 <- ict_pun_selected %>%
  filter(ScientificName == "Ictalurus punctatus" & CatalogNumber == "37883")

# 115142 not included

ict_pun_115142 <- ict_pun_selected %>%
  filter(ScientificName == "Ictalurus punctatus" & CatalogNumber == "115142")

# 147855

ict_pun_147855 <- ict_pun_selected %>%
  filter(ScientificName == "Ictalurus punctatus" & CatalogNumber == "147855")

replace_147855 <- ict_pun_next %>%
  filter(ScientificName == "Ictalurus punctatus" & decade == "1984-1993" & CI == "control")

# Remove the bad lot

ict_pun_selected <- ict_pun_selected %>%
  filter(!(ScientificName == "Ictalurus punctatus" & CatalogNumber == "147855"))

# All other control_1984-1993 are used, so the bad lot cannot be replaced

check <- ict_pun_selected1%>%
  filter(ScientificName == "Ictalurus punctatus" & CatalogNumber == "159967")


# 147873 not included

ict_pun_147873 <- ict_pun_selected %>%
  filter(ScientificName == "Ictalurus punctatus" & CatalogNumber == "147873")



#NOT_ATH
not_ath<-lots_MS_lower %>%
  filter(ScientificName =="Notropis atherinoides")

not_ath$combo<-paste(not_ath$CI,not_ath$decade,sep="_")
not_ath_next<-not_ath[which(not_ath$IndividualCount >= 8),]

not_ath_matrix<-not_ath_next %>%
  group_by(CI, decade) %>%
  summarize(total_request = n())

control_1954to1963<-not_ath_next[ sample( which( not_ath_next$combo == "control_1954-1963"), 2, replace = F), ]#not enough for 5
control_1964to1973<-not_ath_next[ sample( which( not_ath_next$combo == "control_1964-1973"), 5, replace = F), ]
control_1974to1983<-not_ath_next[ sample( which( not_ath_next$combo == "control_1974-1983"), 5, replace = F), ]
control_1984to1993<-not_ath_next[ sample( which( not_ath_next$combo == "control_1984-1993"), 5, replace = F), ]
control_1994to2003<-not_ath_next[ sample( which( not_ath_next$combo == "control_1994-2003"), 5, replace = F), ]
control_2004to2013<-not_ath_next[ sample( which( not_ath_next$combo == "control_2004-2013"), 5, replace = F), ]
impact_1954to1963<-not_ath_next[ sample( which( not_ath_next$combo == "impact_1954-1963"), 2, replace = F), ]#not enough for 5
impact_1964to1973<-not_ath_next[ sample( which( not_ath_next$combo == "impact_1964-1973"), 5, replace = F), ]
impact_1974to1983<-not_ath_next[ sample( which( not_ath_next$combo == "impact_1974-1983"), 5, replace = F), ]
impact_1984to1993<-not_ath_next[ sample( which( not_ath_next$combo == "impact_1984-1993"), 5, replace = F), ]
impact_1994to2003<-not_ath_next[ sample( which( not_ath_next$combo == "impact_1994-2003"), 5, replace = F), ]
impact_2004to2013<-not_ath_next[ sample( which( not_ath_next$combo == "impact_2004-2013"), 5, replace = F), ]

not_ath_selected<-rbind.data.frame(control_1954to1963,control_1964to1973,control_1974to1983,control_1984to1993,
                                   control_1994to2003,control_2004to2013,
                                   impact_1954to1963,impact_1964to1973,impact_1974to1983,impact_1984to1993,
                                   impact_1994to2003,impact_2004to2013)

not_ath_matrix<-not_ath_selected %>%
  group_by(combo) %>%
  summarize(total_request = n())

plot(jitter(not_ath_selected$Latitude,5)~not_ath_selected$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)



### PERVIG
# Katie went through all in-bounds PERVIG, CARVEL, ICTPUN, and PIMVIG and found all valid lots (i.e., present, big enough)
# Now I need to randomly select the lots to be dissected from among these valid lots.

# Here is the code to generate a list of all of the "in-play" lots.

katie_valid<-read.csv("lot_selection/old/lots_suggested_for_dissection_ROUND3 - lots_suggested_for_dissection_ROUND3.csv")

katie_valid_PERVIG <- katie_valid %>%
  filter(ScientificName == "Percina vigil")

# We are including maybes because there are way too few jars in the control without them.

katie_valid_PERVIG <- katie_valid_PERVIG %>%
  filter(Valid.Invalid == "VALID" | Valid.Invalid == "maybe")

per_vig_matrix<-katie_valid_PERVIG %>%
  group_by(CI, decade) %>%
  summarize(total_available = n())

katie_valid_PERVIG$combo<-paste(katie_valid_PERVIG$CI,katie_valid_PERVIG$decade,sep="_")

#control_1954to1963<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "control_1954-1963"), 3, replace = F), ]#not enough for 5
control_1964to1973<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "control_1964-1973"), 3, replace = F), ]
control_1974to1983<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "control_1974-1983"), 3, replace = F), ]
#control_1984to1993<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "control_1984-1993"), 1, replace = F), ]#not enough for 5
control_1994to2003<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "control_1994-2003"), 2, replace = F), ]#not enough for 5
#control_2004to2013<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "control_2004-2013"), 2, replace = F), ]#not enough for 5
impact_1954to1963<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "impact_1954-1963"), 5, replace = F), ]
impact_1964to1973<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "impact_1964-1973"), 5, replace = F), ]
impact_1974to1983<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "impact_1974-1983"), 5, replace = F), ]
impact_1984to1993<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "impact_1984-1993"), 5, replace = F), ]
impact_1994to2003<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "impact_1994-2003"), 5, replace = F), ]
#impact_2004to2013<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "impact_2004-2013"), 4, replace = F), ]#don't bother because we don't have this decade for control


per_vig_selected<-rbind.data.frame(control_1964to1973,control_1974to1983,
                                   control_1994to2003,
                                   impact_1954to1963,impact_1964to1973,impact_1974to1983,impact_1984to1993,
                                   impact_1994to2003)

per_vig_matrix<-per_vig_selected %>%
  group_by(combo) %>%
  summarize(total_request = n())

plot(jitter(per_vig_selected$Latitude,5)~per_vig_selected$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)


### CARVEL
# Katie went through all in-bounds PERVIG, CARVEL, ICTPUN, and PIMVIG and found all valid lots (i.e., present, big enough)
# Now I need to randomly select the lots to be dissected from among these valid lots.

# Here is the code to generate a list of all of the "in-play" lots.

katie_valid<-read.csv("lot_selection/old/lots_suggested_for_dissection_ROUND3 - lots_suggested_for_dissection_ROUND3.csv")

katie_valid_CARVEL <- katie_valid %>%
  filter(ScientificName == "Carpiodes velifer")

# We are including maybes because there are way too few jars in the control without them.

katie_valid_CARVEL <- katie_valid_CARVEL %>%
  filter(Valid.Invalid == "VALID" | Valid.Invalid == "maybe")

per_vig_matrix<-katie_valid_CARVEL %>%
  group_by(CI, decade) %>%
  summarize(total_available = n())

katie_valid_CARVEL$combo<-paste(katie_valid_CARVEL$CI,katie_valid_CARVEL$decade,sep="_")

control_1954to1963<-katie_valid_CARVEL[ sample( which( katie_valid_CARVEL$combo == "control_1954-1963"), 5, replace = F), ]
control_1964to1973<-katie_valid_CARVEL[ sample( which( katie_valid_CARVEL$combo == "control_1964-1973"), 5, replace = F), ]
control_1974to1983<-katie_valid_CARVEL[ which( katie_valid_CARVEL$combo == "control_1974-1983"), ]
control_1984to1993<-katie_valid_CARVEL[ sample( which( katie_valid_CARVEL$combo == "control_1984-1993"), 5, replace = F), ]
control_1994to2003<-katie_valid_CARVEL[ which( katie_valid_CARVEL$combo == "control_1994-2003"), ]
control_2004to2013<-katie_valid_CARVEL[ sample( which( katie_valid_CARVEL$combo == "control_2004-2013"), 2, replace = F), ]
#impact_1954to1963<-katie_valid_CARVEL[ sample( which( katie_valid_CARVEL$combo == "impact_1954-1963"), 5, replace = F), ]
impact_1964to1973<-katie_valid_CARVEL[ sample( which( katie_valid_CARVEL$combo == "impact_1964-1973"), 5, replace = F), ]
impact_1974to1983<-katie_valid_CARVEL[ which( katie_valid_CARVEL$combo == "impact_1974-1983"), ]
impact_1984to1993<-katie_valid_CARVEL[ sample( which( katie_valid_CARVEL$combo == "impact_1984-1993"), 2, replace = F), ]
impact_1994to2003<-katie_valid_CARVEL[ sample( which( katie_valid_CARVEL$combo == "impact_1994-2003"), 3, replace = F), ]
impact_2004to2013<-katie_valid_CARVEL[ which( katie_valid_CARVEL$combo == "impact_2004-2013"), ]


car_vel_selected<-rbind.data.frame(control_1954to1963,control_1964to1973,control_1974to1983,
                                   control_1984to1993,control_1994to2003,control_2004to2013,
                                   impact_1964to1973,impact_1974to1983,impact_1984to1993,
                                   impact_1994to2003,impact_2004to2013)

car_vel_matrix<-car_vel_selected %>%
  group_by(combo) %>%
  summarize(total_request = n())

write.csv(car_vel_selected, file="lot_selection/final_lots/CAR_VEL_final_lots_2024.02.07.csv")

plot(jitter(car_vel_selected$Latitude,5)~car_vel_selected$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)


### PIMVIG
# Katie went through all in-bounds PERVIG, CARVEL, ICTPUN, and PIMVIG and found all valid lots (i.e., present, big enough)
# Now I need to randomly select the lots to be dissected from among these valid lots.

# Here is the code to generate a list of all of the "in-play" lots.

katie_valid<-read.csv("lot_selection/old/lots_suggested_for_dissection_ROUND3 - lots_suggested_for_dissection_ROUND3.csv")

katie_valid_PIMVIG <- katie_valid %>%
  filter(ScientificName == "Pimephales vigilax")

# We are including maybes because there are way too few jars in the control without them.

katie_valid_PIMVIG <- katie_valid_PIMVIG %>%
  filter(Valid.Invalid == "VALID" | Valid.Invalid == "maybe")

pim_vig_matrix<-katie_valid_PIMVIG %>%
  group_by(CI, decade) %>%
  summarize(total_available = n())

katie_valid_PIMVIG$combo<-paste(katie_valid_PIMVIG$CI,katie_valid_PIMVIG$decade,sep="_")

control_1954to1963<-katie_valid_PIMVIG[ sample( which( katie_valid_PIMVIG$combo == "control_1954-1963"), 3, replace = F), ]
control_1964to1973<-katie_valid_PIMVIG[ sample( which( katie_valid_PIMVIG$combo == "control_1964-1973"), 5, replace = F), ]
control_1974to1983<-katie_valid_PIMVIG[ sample( which( katie_valid_PIMVIG$combo == "control_1974-1983"), 5, replace = F), ]
control_1984to1993<-katie_valid_PIMVIG[ sample( which( katie_valid_PIMVIG$combo == "control_1984-1993"), 5, replace = F), ]
control_1994to2003<-katie_valid_PIMVIG[ sample( which( katie_valid_PIMVIG$combo == "control_1994-2003"), 5, replace = F), ]
control_2004to2013<-katie_valid_PIMVIG[ which( katie_valid_PIMVIG$combo == "control_2004-2013"), ]
#impact_1954to1963<-katie_valid_PIMVIG[ sample( which( katie_valid_PIMVIG$combo == "impact_1954-1963"), 5, replace = F), ]
impact_1964to1973<-katie_valid_PIMVIG[ sample( which( katie_valid_PIMVIG$combo == "impact_1964-1973"), 5, replace = F), ]
impact_1974to1983<-katie_valid_PIMVIG[ sample( which( katie_valid_PIMVIG$combo == "impact_1974-1983"), 5, replace = F), ]
impact_1984to1993<-katie_valid_PIMVIG[ sample( which( katie_valid_PIMVIG$combo == "impact_1984-1993"), 5, replace = F), ]
impact_1994to2003<-katie_valid_PIMVIG[ sample( which( katie_valid_PIMVIG$combo == "impact_1994-2003"), 5, replace = F), ]
impact_2004to2013<-katie_valid_PIMVIG[ sample( which( katie_valid_PIMVIG$combo == "impact_2004-2013"), 2, replace = F), ]


pim_vig_selected<-rbind.data.frame(control_1954to1963,control_1964to1973,control_1974to1983,
                                   control_1984to1993,control_1994to2003,control_2004to2013,
                                   impact_1964to1973,impact_1974to1983,impact_1984to1993,
                                   impact_1994to2003,impact_2004to2013)

pim_vig_matrix<-pim_vig_selected %>%
  group_by(combo) %>%
  summarize(total_request = n())

plot(jitter(pim_vig_selected$Latitude,5)~pim_vig_selected$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)




### Now bring everything together into the main dataset - the newly selected PERVIG, CARVEL, and PIMVIG plus
### everything from the main dataset. Remember that you already replaced the bad ICTPUN above, so that's sorted.


car_vel_final <- car_vel_selected %>%
  mutate(number_individuals_requested = 4)

gam_aff_final <- gam_aff_selected %>%
  mutate(number_individuals_requested = 4)

hyb_nuc_final <- hyb_nuc_selected %>%
  mutate(number_individuals_requested = 4)

ict_pun_final <- ict_pun_selected %>%
  mutate(number_individuals_requested = 2)

not_ath_final <- not_ath_selected %>%
  mutate(number_individuals_requested = 4)

per_vig_final <- per_vig_selected %>%
  mutate(number_individuals_requested = 2)

pim_vig_final <- pim_vig_selected %>%
  mutate(number_individuals_requested = 4)


### Trim so that the datasets have the same number of columns

ncol(car_vel_final)
ncol(gam_aff_final)
ncol(hyb_nuc_final)
ncol(ict_pun_final)
ncol(not_ath_final)
ncol(per_vig_final)
ncol(pim_vig_final)

per_vig_final<-per_vig_final[ , -which(names(per_vig_final) %in% c("X.1","X","Valid.Invalid","notes"))]
car_vel_final<-car_vel_final[ , -which(names(car_vel_final) %in% c("X.1","X","Valid.Invalid","notes"))]
pim_vig_final<-pim_vig_final[ , -which(names(pim_vig_final) %in% c("X.1","X","Valid.Invalid","notes"))]

final_dataset <- rbind(car_vel_final, gam_aff_final, hyb_nuc_final, ict_pun_final, not_ath_final, per_vig_final, pim_vig_final)


### Now write the final file.

write.csv(final_lots, file="lot_selection/final_lots/final_lots_2024.02.07.csv")


# Prior to Katie's in-person winnowing
tally_up<-lots_MS_number_available_to_dissect %>%
  group_by(ScientificName) %>%
  summarize(total_lots = n(), total_individuals_requested = sum(number_individuals_requested))
tally_up

sum(tally_up$total_individuals_requested)
sum(tally_up$total_lots)

# After Katie's in-person winnowing
tally_up<-final_lots %>%
  group_by(ScientificName) %>%
  summarize(total_lots = n(), total_individuals_requested = sum(number_individuals_requested))
tally_up

sum(tally_up$total_individuals_requested)
sum(tally_up$total_lots)

