# Testing BACI response of parasites to pulp mill / Clean Water Act
# Written by Chelsea Wood (chelwood@uw.edu)
# August 2024

# Read in the dataset

full_dataset_with_LH<-read.csv("data/processed/Full_dataset_with_psite_life_history_info_2024.08.01.csv")


# Take out myxos, since they were counted differently (you'll run a parallel analysis for them)
minus_myxos <- full_dataset_with_LH %>%
  filter(Life_History!="Myxozoa")


# Test how things changed through time

library(MASS)
time_model_draft<-glm.nb(psite_count~YearCollected,data=minus_myxos)
summary(time_model_draft)

str(minus_myxos)

model_draft<-glmer.nb(psite_count~CI*before_after+scaled_TL_mm+(1|Fish_sp.x/psite_spp.x)+(1|CatalogNumber),
                      data=minus_myxos,family="nbinom")
summary(model_draft)


library(ggeffects)
big_predictions<-ggeffect(model_draft,c("before_after","CI"))

big_predictions$group

big_plot<-ggplot(big_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=big_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("before or after CWA")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=(rev(levels(big_predictions$x))))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
big_plot


