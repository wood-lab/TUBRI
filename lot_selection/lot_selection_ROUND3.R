# Selecting specimens for TUBRI project
# Chelsea Wood
# July 2023

library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales) #for transparency
library(tidyverse)
library(ggmap)
library(cowplot)


# We tried to dig up all of the lots for Round 1 (see lot_selection_ROUND1.R and 
# lots_suggested_for_dissection_ROUND1.csv), but many are missing or the fish are too small.

# Katie is currently in the collection and we're going to try another approach. I will generate
# a list of all the "in-play" lots (i.e., all lots that are suitable).

# She will then find all of those lots and pull all of the ones that are present and big enough.
# I will then take her list and randomly select specimens to maximize sampling across the four 
# quadrants.

# Here is the code to generate a list of all of the "in-play" lots.

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

bounds<-c(left=-89.87, bottom=30.67, right=-89.79, top=30.81)
map_check<-get_stamenmap(bounds, zoom=11, maptype = "terrain-background") %>% ggmap()+
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

map_check<-get_stamenmap(bounds, zoom=12, maptype = "terrain-background") %>% ggmap()+
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


# Cool, now you're ready to select lots.
# If fewer than 8 individuals are in the lot, throw it out (unless it's a big fish and there are few
# individuals per lot - then you can change the criterion to 4 fish per lot).  

# CAR_VEL
car_vel<-lots_MS_lower %>%
  filter(ScientificName =="Carpiodes velifer")

car_vel<-car_vel[which(car_vel$IndividualCount >= 4),]

car_vel_plot<-plot(jitter(car_vel$Latitude,5)~car_vel$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)
car_vel_plot

car_vel_plot<-ggplot(data=car_vel,aes(x=YearCollected,y=jitter(Latitude,5)))+
  annotate(geom = "rect", xmin = -Inf, xmax = 1973, ymin = -Inf, ymax = Inf,
           fill = "royalblue3", colour = NA, alpha = 0.7)+
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 30.76,
           fill = "red2", colour = NA, alpha = 0.5)+
  geom_point(color="black",fill="darkgray",shape=21,size=4)+
  geom_hline(yintercept = 30.76, lty = 2)+
  geom_vline(xintercept = 1973, lty = 2)+
  xlab("")+
  ylab("Latitude (°N)")+
  ggtitle("Carpiodes velifer")+
  #ylim()+
  #xlim()+
  theme_bw()+
  theme(plot.title=element_text(face="italic",hjust = 0.5,size=20),
        plot.margin = unit(c(0,0,0,0), "cm"), text=element_text(family='sans',size=20),
        axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),legend.text = element_text(size=14),
        legend.title = element_text(size=14))
car_vel_plot


# GAM_AFF
gam_aff<-lots_MS_lower %>%
  filter(ScientificName =="Gambusia affinis")

gam_aff<-gam_aff[which(gam_aff$IndividualCount >= 8),]

plot(jitter(gam_aff$Latitude,5)~gam_aff$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

gam_aff_plot<-ggplot(data=gam_aff,aes(x=YearCollected,y=jitter(Latitude,5)))+
  annotate(geom = "rect", xmin = -Inf, xmax = 1973, ymin = -Inf, ymax = Inf,
           fill = "royalblue3", colour = NA, alpha = 0.7)+
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 30.76,
           fill = "red2", colour = NA, alpha = 0.5)+
  geom_point(color="black",fill="gray",shape=21,size=4)+
  geom_hline(yintercept = 30.76, lty = 2)+
  geom_vline(xintercept = 1973, lty = 2)+
  xlab("")+
  ylab("")+
  ggtitle("Gambusia affinis")+
  #ylim()+
  #xlim()+
  theme_bw()+
  theme(plot.title=element_text(face="italic",hjust = 0.5,size=20),
        plot.margin = unit(c(0,0,0,0), "cm"), text=element_text(family='sans',size=20),
        axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),legend.text = element_text(size=14),
        legend.title = element_text(size=14))
gam_aff_plot


#HYB_NUC
hyb_nuc<-lots_MS_lower %>%
  filter(ScientificName =="Hybognathus nuchalis")

hyb_nuc<-hyb_nuc[which(hyb_nuc$IndividualCount >= 8),]

plot(jitter(hyb_nuc$Latitude,5)~hyb_nuc$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

hyb_nuc_plot<-ggplot(data=hyb_nuc,aes(x=YearCollected,y=jitter(Latitude,5)))+
  annotate(geom = "rect", xmin = -Inf, xmax = 1973, ymin = -Inf, ymax = Inf,
           fill = "royalblue3", colour = NA, alpha = 0.7)+
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 30.76,
           fill = "red2", colour = NA, alpha = 0.5)+
  geom_point(color="black",fill="gray",shape=21,size=4)+
  geom_hline(yintercept = 30.76, lty = 2)+
  geom_vline(xintercept = 1973, lty = 2)+
  xlab("")+
  ylab("")+
  ggtitle("Hybognathus nuchalis")+
  #ylim()+
  #xlim()+
  theme_bw()+
  theme(plot.title=element_text(face="italic",hjust = 0.5,size=20),
        plot.margin = unit(c(0,0,0,0), "cm"), text=element_text(family='sans',size=20),
        axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),legend.text = element_text(size=14),
        legend.title = element_text(size=14))
hyb_nuc_plot


#ICT_PUN
ict_pun<-lots_MS_lower %>%
  filter(ScientificName =="Ictalurus punctatus")

# This is a big fish and there are fewer specimens per lot than others, so lower the threshold
ict_pun<-ict_pun[which(ict_pun$IndividualCount >= 4),]

plot(jitter(ict_pun$Latitude,5)~ict_pun$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

ict_pun_plot<-ggplot(data=ict_pun,aes(x=YearCollected,y=jitter(Latitude,5)))+
  annotate(geom = "rect", xmin = -Inf, xmax = 1973, ymin = -Inf, ymax = Inf,
           fill = "royalblue3", colour = NA, alpha = 0.7)+
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 30.76,
           fill = "red2", colour = NA, alpha = 0.5)+
  geom_point(color="black",fill="gray",shape=21,size=4)+
  geom_hline(yintercept = 30.76, lty = 2)+
  geom_vline(xintercept = 1973, lty = 2)+
  xlab("")+
  ylab("")+
  ggtitle("Ictalurus punctatus")+
  #ylim()+
  #xlim()+
  theme_bw()+
  theme(plot.title=element_text(face="italic",hjust = 0.5,size=20),
        plot.margin = unit(c(0,0,0,0), "cm"), text=element_text(family='sans',size=20),
        axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),legend.text = element_text(size=14),
        legend.title = element_text(size=14))
ict_pun_plot


#NOT_ATH
not_ath<-lots_MS_lower %>%
  filter(ScientificName =="Notropis atherinoides")

not_ath<-not_ath[which(not_ath$IndividualCount >= 8),]

plot(jitter(not_ath$Latitude,5)~not_ath$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

not_ath_plot<-ggplot(data=not_ath,aes(x=YearCollected,y=jitter(Latitude,5)))+
  annotate(geom = "rect", xmin = -Inf, xmax = 1973, ymin = -Inf, ymax = Inf,
           fill = "royalblue3", colour = NA, alpha = 0.7)+
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 30.76,
           fill = "red2", colour = NA, alpha = 0.5)+
  geom_point(color="black",fill="gray",shape=21,size=4)+
  geom_hline(yintercept = 30.76, lty = 2)+
  geom_vline(xintercept = 1973, lty = 2)+
  xlab("")+
  ylab("Latitude (°N)")+
  ggtitle("Notropis atherinoides")+
  #ylim()+
  #xlim()+
  theme_bw()+
  theme(plot.title=element_text(face="italic",hjust = 0.5,size=20),
        plot.margin = unit(c(0,0,0,0), "cm"), text=element_text(family='sans',size=20),
        axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),legend.text = element_text(size=14),
        legend.title = element_text(size=14))
not_ath_plot


#PER_VIG
per_vig<-lots_MS_lower %>%
  filter(ScientificName =="Percina vigil")

#this is a big fish and there are fewer specimens per lot than others, so lower the threshold

per_vig<-per_vig[which(per_vig$IndividualCount >= 4),]

plot(jitter(per_vig$Latitude,5)~per_vig$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

per_vig_plot<-ggplot(data=per_vig,aes(x=YearCollected,y=jitter(Latitude,5)))+
  annotate(geom = "rect", xmin = -Inf, xmax = 1973, ymin = -Inf, ymax = Inf,
           fill = "royalblue3", colour = NA, alpha = 0.7)+
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 30.76,
           fill = "red2", colour = NA, alpha = 0.5)+
  geom_point(color="black",fill="gray",shape=21,size=4)+
  geom_hline(yintercept = 30.76, lty = 2)+
  geom_vline(xintercept = 1973, lty = 2)+
  xlab("")+
  ylab("")+
  ggtitle("Percina vigil")+
  #ylim()+
  #xlim()+
  theme_bw()+
  theme(plot.title=element_text(face="italic",hjust = 0.5,size=20),
        plot.margin = unit(c(0,0,0,0), "cm"), text=element_text(family='sans',size=20),
        axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),legend.text = element_text(size=14),
        legend.title = element_text(size=14))
per_vig_plot


#PIM_VIG
pim_vig<-lots_MS_lower %>%
  filter(ScientificName =="Pimephales vigilax")

pim_vig<-pim_vig[which(pim_vig$IndividualCount >= 8),]

plot(jitter(pim_vig$Latitude,5)~pim_vig$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)

pim_vig_plot<-ggplot(data=pim_vig,aes(x=YearCollected,y=jitter(Latitude,5)))+
  annotate(geom = "rect", xmin = -Inf, xmax = 1973, ymin = -Inf, ymax = Inf,
           fill = "royalblue3", colour = NA, alpha = 0.7)+
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 30.76,
           fill = "red2", colour = NA, alpha = 0.5)+
  geom_point(color="black",fill="gray",shape=21,size=4)+
  geom_hline(yintercept = 30.76, lty = 2)+
  geom_vline(xintercept = 1973, lty = 2)+
  xlab("")+
  ylab("")+
  ggtitle("Pimephales vigilax")+
  #ylim()+
  #xlim()+
  theme_bw()+
  theme(plot.title=element_text(face="italic",hjust = 0.5,size=20),
        plot.margin = unit(c(0,0,0,0), "cm"), text=element_text(family='sans',size=20),
        axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),legend.text = element_text(size=14),
        legend.title = element_text(size=14))
pim_vig_plot


# Compile request for MS_lower

lots_MS_number_available_to_dissect<-rbind.data.frame(car_vel,gam_aff,hyb_nuc,ict_pun,not_ath,per_vig,pim_vig)


# Plot request

final_figure <- ggdraw(plot=NULL,xlim=c(0,40),ylim=c(0,20))+
  draw_plot(car_vel_plot,x=0,y=10,width=10,height=10)+
  draw_plot(gam_aff_plot,x=10,y=10,width=10,height=10)+
  draw_plot(hyb_nuc_plot,x=20,y=10,width=10,height=10)+
  draw_plot(ict_pun_plot,x=30,y=10,width=10,height=10)+
  draw_plot(not_ath_plot,x=0,y=0,width=10,height=10)+
  draw_plot(per_vig_plot,x=10,y=0,width=10,height=10)+
  draw_plot(pim_vig_plot,x=20,y=0,width=10,height=10)+
  draw_label("Year",x=20,y=0.45,size=20)
final_figure

dev.off()


# Katie needs replacements for some ICTPUN specimens

View(ict_pun)

ict_pun_34760 <- lots_MS_number_available_to_dissect %>%
  filter(ScientificName == "Ictalurus punctatus" & CatalogNumber == "34760")

ict_pun_replace_34760 <- lots_MS_number_available_to_dissect %>%
  filter(ScientificName == "Ictalurus punctatus" & decade == "1964-1973" & CI == "control")

write.csv(ict_pun_replace_34760, file="lot_selection/replacements/ICTPUN/ict_pun_replace_34760.csv")

ict_pun_37883 <- lots_MS_number_available_to_dissect %>%
  filter(ScientificName == "Ictalurus punctatus" & CatalogNumber == "37883")

write.csv(ict_pun_replace_34760, file="lot_selection/replacements/ICTPUN/ict_pun_replace_34760_37883.csv")

ict_pun_115142 <- lots_MS_number_available_to_dissect %>%
  filter(ScientificName == "Ictalurus punctatus" & CatalogNumber == "115142")

ict_pun_replace_115142 <- lots_MS_number_available_to_dissect %>%
  filter(ScientificName == "Ictalurus punctatus" & decade == "1974-1983" & CI == "control")

write.csv(ict_pun_replace_34760, file="lot_selection/replacements/ICTPUN/ict_pun_replace_115142.csv")

ict_pun_147855 <- lots_MS_number_available_to_dissect %>%
  filter(ScientificName == "Ictalurus punctatus" & CatalogNumber == "147855")

ict_pun_replace_147855 <- lots_MS_number_available_to_dissect %>%
  filter(ScientificName == "Ictalurus punctatus" & decade == "1984-1993" & CI == "control")

write.csv(ict_pun_replace_147855, file="lot_selection/replacements/ICTPUN/ict_pun_replace_147855.csv")

ict_pun_147873 <- lots_MS_number_available_to_dissect %>%
  filter(ScientificName == "Ictalurus punctatus" & CatalogNumber == "147873")

ict_pun_replace_147873 <- lots_MS_number_available_to_dissect %>%
  filter(ScientificName == "Ictalurus punctatus" & decade == "1984-1993" & CI == "impact")

write.csv(ict_pun_replace_147873, file="lot_selection/replacements/ICTPUN/ict_pun_replace_147873.csv")


# Almost all the ICTPUN to be dissected have been pulled.  We just needed to replace a handful of them.
# That's what the code above is doing.  We haven't yet pulled the replacements, but once we do, we will have
# ICTPUN fully ready to go!

# All that you need to do is compile the list of replacements above, put it in numerical order, pull them all,
# and then randomly select amongst them.



# Katie went through all in-bounds PERVIG and found all valid lots (i.e., present, big enough)
# Now I need to randomly select the lots to be dissected from among these valid lots.

# Here is the code to generate a list of all of the "in-play" lots.

katie_valid_PERVIG<-read.csv("lot_selection/lots_suggested_for_dissection_LOTS_KATIE_IDS_AS_VALID.csv")

katie_valid_PERVIG <- katie_valid_PERVIG %>%
  filter(ScientificName == "Percina vigil")

# We are including maybes because there are way too few jars in the control without them.

katie_valid_PERVIG <- katie_valid_PERVIG %>%
  filter(Valid_Invalid == "VALID" | Valid_Invalid == "maybe")

per_vig_matrix<-katie_valid_PERVIG %>%
  group_by(CI, decade) %>%
  summarize(total_available = n())

katie_valid_PERVIG$combo<-paste(katie_valid_PERVIG$CI,katie_valid_PERVIG$decade,sep="_")

#control_1954to1963<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "control_1954-1963"), 3), ]#not enough for 5
control_1964to1973<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "control_1964-1973"), 3), ]
control_1974to1983<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "control_1974-1983"), 3), ]
#control_1984to1993<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "control_1984-1993"), 1), ]#not enough for 5
control_1994to2003<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "control_1994-2003"), 2), ]#not enough for 5
#control_2004to2013<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "control_2004-2013"), 2), ]#not enough for 5
impact_1954to1963<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "impact_1954-1963"), 5), ]
impact_1964to1973<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "impact_1964-1973"), 5), ]
impact_1974to1983<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "impact_1974-1983"), 5), ]
impact_1984to1993<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "impact_1984-1993"), 5), ]
impact_1994to2003<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "impact_1994-2003"), 5), ]
#impact_2004to2013<-katie_valid_PERVIG[ sample( which( katie_valid_PERVIG$combo == "impact_2004-2013"), 4), ]#don't bother because we don't have this decade for control


per_vig_selected<-rbind.data.frame(control_1964to1973,control_1974to1983,
                                   control_1994to2003,
                                   impact_1954to1963,impact_1964to1973,impact_1974to1983,impact_1984to1993,
                                   impact_1994to2003)

per_vig_matrix<-per_vig_selected %>%
  group_by(combo) %>%
  summarize(total_request = n())

plot(jitter(per_vig_selected$Latitude,5)~per_vig_selected$YearCollected)+abline(a = 30.76, b = 0, lty = 2)+abline(v = 1973, lty = 2)





# Export the sheet

write.csv(per_vig_selected, file="lot_selection/final_lots/per_vig_inbound_valid_selected.csv")

# The sheet above is ready to go - these are the FINAL lots that will be targeted for dissection.
# The maybes and valids are all pulled so they should be easy to find when you get to TUBRI!


####### HERE'S WHERE WE'RE AT (2 Aug 2023) ####### 

# PER_VIG = just a few to pull and then ready to go!
# ICT_PUN = just a few to pull and then ready to go! (going through the jars for replacements is going to be
# time consuming because there are sooooo many jars)
# CAR_VEL = held in a separate building without AC, now moved to the main building. 100% done and stored in the main builing!

# GAM_AFF = We need to send to Justin a list of the NAs to see if there are bigger jars floating about that
# match those lot numbers /// All in-bounds (maybes/VALIDs) have been pulled, we just need to select from 
# amongst those.

# HYB_NUC = Looked through all the big jars (not the 4oz jars, which are likely to be way too small).
# All in-bounds (maybes/VALIDs) have been pulled, we just need to select from amongst those.

# PIM_VIG = Looked through all the big jars (not the 4oz jars, which are likely to be way too small).
# All in-bounds (maybes/VALIDs) have been pulled, we just need to select from amongst those.

# NOT_ATH = Looked through all the big jars (not the 4oz jars, which are likely to be way too small).
# All in-bounds (maybes/VALIDs) have been pulled, we just need to select from amongst those.

# Have the entire team go through the lot selection process together with the entire group, 
# maybe in week 1 as we're getting our feet under us, trading off with parasite dissection training.


# Practice dissections should not be counted since they weren't always on randomly selected lots

