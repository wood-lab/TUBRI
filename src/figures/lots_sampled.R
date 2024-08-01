# SuPeR Parasites map and plot
# Creating a map and plot just like the one in the CAREER proposal, but updated to reflect the lat, long, 
# and year of all the samples we ACTUALLY processed (not just the ones that are available at TUBRI).

library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales) #for transparency
library(rgdal)
library(maps)
library(mapdata)
library(maptools) #for shapefiles
library(scales) #for transparency
library(ggmap)
library(ggsn)
library(tidyverse)
library(tidyr)
library(dplyr)
library(ggmap)
library(cowplot)
library(sp)
library(rgeos)
library(raster)
library(maptools)
library(rgdal)
library(magick)
library(pdftools)


# Load data about TUBRI lots

lots<-read.csv("data/processed/Full_dataset_with_psite_life_history_info_2024.08.01.csv", header = T, sep = ",")
str(lots)

# MS lower map
survey_sites_MS_lower<-read.csv("TUBRI/MS/survey_sites_MS_lower.csv")
bounds<-c(left=-89.855, bottom=30.69, right=-89.81, top=30.81)
map_MS_lower<-get_stamenmap(bounds, zoom=12, maptype = "terrain-background") %>% ggmap()+
  geom_point(y=30.76558198250576, x=-89.83398463056083, size=4)+
  geom_point(data = survey_sites_MS_lower, aes(x=Longitude,y=Latitude,fill=CI),shape=21,size=4)+
  scale_fill_manual(values=c("white","#bdbdbd"))+
  xlab("")+
  ylab("")+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=4),axis.text.y=element_text(size=8))
#geom_rect(xmin=-106.89,ymin=34.8,xmax=-106.45,ymax=35.2102673,fill="grey",alpha=0.2,linetype="blank")+
#geom_hline(yintercept=35.2102673,linetype="dashed")+
#geom_map(data=water_AL_df,map=water_AL_df,aes(x=long,y=lat,map_id=id),color="lightsteelblue3",fill="lightsteelblue3")
#geom_map(data=polygon_df,map=polygon_df,aes(x=long,y=lat,map_id=id),color="black",fill=NA)+
#annotate("text",x=-106.6,y=35.097,label="City of Albuquerque",size=5)+
#annotate("text",x=-106.58,y=35.33,label="Rio Grande River",size=5,angle=55)
map_MS_lower

# create legend
survey_sites_legend<-read.csv("TUBRI/AL/survey_sites_legend.csv")
bounds<-c(left=-87.6, bottom=31.5, right=-87.35, top=32.05)
legend<-get_stamenmap(bounds, zoom=11, maptype = "terrain-background") %>% ggmap()+
  geom_point(data = survey_sites_legend, aes(x=Longitude,y=Latitude,fill=CI),shape=21)+
  scale_fill_manual(values=c("white","black","#bdbdbd"),limits=c("control","mill","impact"))+
  xlab("")+
  ylab("")+
  theme(legend.position = "right", legend.title = element_blank(), legend.key=element_blank(), legend.text = element_text(size=20))+
  guides(fill = guide_legend(override.aes = list(size=7)))+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=4),axis.text.y=element_text(size=8))
#geom_rect(xmin=-106.89,ymin=34.8,xmax=-106.45,ymax=35.2102673,fill="grey",alpha=0.2,linetype="blank")+
#geom_hline(yintercept=35.2102673,linetype="dashed")+
#geom_map(data=water_AL_df,map=water_AL_df,aes(x=long,y=lat,map_id=id),color="lightsteelblue3",fill="lightsteelblue3")
#geom_map(data=polygon_df,map=polygon_df,aes(x=long,y=lat,map_id=id),color="black",fill=NA)+
#annotate("text",x=-106.6,y=35.097,label="City of Albuquerque",size=5)+
#annotate("text",x=-106.58,y=35.33,label="Rio Grande River",size=5,angle=55)
legend


#MS lower plot

lots_MS_lower<-read.csv("TUBRI/MS/MS_lower_search_results.csv")


#Filter out just the five species you want to focus on

lots_MS_lower <- lots_MS_lower %>%
  filter(ScientificName == "Pimephales vigilax" | ScientificName == "Notropis atherinoides" | ScientificName == "Gambusia affinis" 
         | ScientificName == "Ictalurus punctatus" | ScientificName == "Lepomis macrochirus" | ScientificName == "Carpiodes velifer" 
         | ScientificName == "Hybognathus nuchalis" | ScientificName == "Micropterus salmoides" 
         | ScientificName == "Lepomis megalotis" | ScientificName == "Percina vigil")

min(lots_MS_lower$YearCollected)
max(lots_MS_lower$YearCollected)
thing<-sort(lots_MS_lower$YearCollected)
2005-1963

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

library(wesanderson)
wes_palette("Zissou1")
pal<-wes_palette(name = "Zissou1", 10, type = "continuous")

plot_MS_lower<-ggplot(data=lots_MS_lower,aes(x=YearCollected,y=CI))+
  scale_y_discrete(limits=rev,labels=c("P","C"))+
  annotate(geom = "rect", xmin = -Inf, xmax = 1972.5, ymin = -Inf, ymax = Inf,
           fill = "grey", colour = NA, alpha = 0.5)+
  annotate(geom = "rect", xmin = 1966.75, xmax = 1969.25, ymin = -Inf, ymax = Inf,
           fill = "darkgrey", colour = NA, alpha = 0.75)+
  annotate(geom = "rect", xmin = 1975.75, xmax = 1977.25, ymin = -Inf, ymax = Inf,
           fill = "darkgrey", colour = NA, alpha = 0.75)+
  annotate(geom = "rect", xmin = 1980.75, xmax = 1981.25, ymin = -Inf, ymax = Inf,
           fill = "darkgrey", colour = NA, alpha = 0.75)+
  annotate(geom = "rect", xmin = 1984.75, xmax = 1986.25, ymin = -Inf, ymax = Inf,
           fill = "darkgrey", colour = NA, alpha = 0.75)+
  annotate(geom = "rect", xmin = 1999.75, xmax = 2001.25, ymin = -Inf, ymax = Inf,
           fill = "darkgrey", colour = NA, alpha = 0.75)+
  geom_point(data=lots_MS_lower,aes(x=YearCollected,fill=ScientificName),position=position_jitter(width=0.15),shape=21)+
  scale_fill_manual(values=pal[1:10],name='')+
  xlab("")+
  ylab("")+
  geom_vline(xintercept=1972.5,linetype="dashed")+
  #Katrina: geom_vline(xintercept=2005,linetype="dotted")+
  xlim(1962,2006)+
  theme_bw()+
  theme(legend.position = "none")+
  guides(fill = guide_legend(override.aes = list(size=7)))+
  theme(plot.margin = unit(c(0,0,0,-0.6), "cm"), text=element_text(family='sans',size=20),
        axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),legend.text = element_text(size=18))
plot_MS_lower


#Reduce to two sites

final_figure <- ggdraw(plot=NULL,xlim=c(0,15),ylim=c(0,5))+
  draw_image(magick::image_read_pdf("TUBRI/TUBRI_map.pdf"),x=0.1,y=0,width=5,height=4.75)+
  #draw_plot(legend,x=9.3,y=0,width=2,height=5)+
  draw_plot(map_MS_lower,x=5,y=0,width=2,height=5)+
  draw_plot(plot_MS_lower,x=7,y=0,width=3,height=5)+
  draw_plot(map_AL,x=10,y=0,width=2,height=5)+
  draw_plot(plot_AL,x=11.9,y=0,width=3,height=5)+
  draw_label("(a)",x=0.25,y=4.75,size=30)+
  draw_label("(b)",x=5,y=4.75,size=30)+
  draw_label("(c)",x=10.4,y=4.75,size=30)
final_figure