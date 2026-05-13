# Testing BACI response of parasites to pulp mill / Clean Water Act
# Written by Chelsea Wood (chelwood@uw.edu)
# August 2024

# Load libraries
library(tidyverse)
library(MASS)
library(lme4)
library(emmeans)
library(ggeffects)
 

# Read in the dataset

full_dataset_with_LH<-read.csv("data/processed/Full_dataset_with_psite_life_history_info_2024.11.21.csv")


# Take out myxos, since they were counted differently (you'll run a parallel analysis for them).
minus_myxos <- full_dataset_with_LH %>%
  filter(Parasite_taxonomic_group!="Myxozoa")


# You should also remove any ciliates or nematomorphs - ciliates because we're not sure we always detected
# them, nematomorph because I'm not convinced it was a nematomorph and because it isn't a fish parasite
# anyway, merely gut contents.

minus_myxos <- minus_myxos %>%
  filter(Parasite_taxonomic_group!="Nematomorpha")

minus_myxos <- minus_myxos %>%
  filter(Parasite_taxonomic_group!="Protozoa")

levels(as.factor(minus_myxos$Parasite_taxonomic_group))

str(minus_myxos)


# Let's do some initial tallies

total_n_fish_lines <- minus_myxos %>%
  group_by(IndividualFishID,Fish_sp.x) %>%
  summarize(count = n())

total_n_fish <- total_n_fish_lines %>%
  group_by(Fish_sp.x) %>%
  summarize(count = n())

sum(total_n_fish$count)

sum(minus_myxos$psite_count,na.rm=TRUE)

psite_species <- levels(as.factor(minus_myxos$fish_psite_combo))
length(psite_species)


# Now lets remove any parasites that occur at <5% prevalence in their host species.

minus_myxos$positive <- ifelse(minus_myxos$psite_count > 0 , 1, 0)

#### what about NAs?

thing <- minus_myxos %>%
  group_by(fish_psite_combo) %>%
  summarise(pos = sum(as.numeric(positive),na.rm=T), count = length(positive))

thing$prev <- thing$pos/thing$count
view(thing)
length(thing$fish_psite_combo)

high_prev_psites <- thing %>%
  filter(prev>=0.05)
view(high_prev_psites)
length(high_prev_psites$fish_psite_combo)

high_prev_psites$fish_psite_combo

final_dataset <- minus_myxos %>%
  filter(fish_psite_combo == "Carpiodes velifer_CEST.MEGAN" | fish_psite_combo == "Carpiodes velifer_TREM.META.UNK" |
           fish_psite_combo == "Gambusia affinis_MONO.SASE" | fish_psite_combo == "Gambusia affinis_TREM.POS" |
           fish_psite_combo == "Hybognathus nuchalis_MONO.DACT" | fish_psite_combo == "Hybognathus nuchalis_TREM.ASS" |
           fish_psite_combo == "Hybognathus nuchalis_TREM.BC" | fish_psite_combo == "Hybognathus nuchalis_TREM.CM" |
           fish_psite_combo == "Hybognathus nuchalis_TREM.META.GO" | fish_psite_combo == "Hybognathus nuchalis_TREM.META.HET" |
           fish_psite_combo == "Hybognathus nuchalis_TREM.META.SP"| fish_psite_combo == "Ictalurus punctatus_MONO.IP" |
           fish_psite_combo == "Ictalurus punctatus_NEM.NMTS" | fish_psite_combo == "Ictalurus punctatus_NEM.PHAR" |
           fish_psite_combo == "Ictalurus punctatus_NEM.SP" | fish_psite_combo == "Ictalurus punctatus_TREM.ALLO" |
           fish_psite_combo == "Ictalurus punctatus_TREM.LG" | fish_psite_combo == "Ictalurus punctatus_TREM.MEGICT" |
           fish_psite_combo == "Notropis atherinoides_MONO.ALL" | fish_psite_combo == "Notropis atherinoides_TREM.CYST" |
           fish_psite_combo == "Notropis atherinoides_TREM.LARV" | fish_psite_combo == "Notropis atherinoides_TREM.META" |
           fish_psite_combo == "Percina vigil_ACANTH.CYA" | fish_psite_combo == "Percina vigil_NEM.LARV" |
           fish_psite_combo == "Pimephales vigilax_META.UNK" | fish_psite_combo == "Pimephales vigilax_MONO.DACT" |
           fish_psite_combo == "Pimephales vigilax_TREM.BUC" | fish_psite_combo == "Pimephales vigilax_TREM.MET" |
           fish_psite_combo == "Pimephales vigilax_TREM.METAF" | fish_psite_combo == "Pimephales vigilax_TREM.METAGS" |
           fish_psite_combo == "Pimephales vigilax_TREM.SSS")

levels(as.factor(final_dataset$fish_psite_combo))

sum(final_dataset$psite_count,na.rm=TRUE)


# Test how things changed through time

time_model_draft<-glm.nb(psite_count~YearCollected,data=final_dataset)
summary(time_model_draft)

model_draft_glm_CI<-glm.nb(psite_count~YearCollected*CI+Latitude+Fish_sp.x*TotalLength_mm,data=final_dataset)
summary(model_draft_glm_CI)

model_draft_glm_LH<-glm.nb(psite_count~YearCollected*Life_History+Latitude+Fish_sp.x*TotalLength_mm,data=final_dataset)
summary(model_draft_glm_LH)

model_draft<-glmer.nb(psite_count~CI*before_after*Life_History+Fish_sp.x*scale(TotalLength_mm)+
                        (1|fish_psite_combo)+(1|CatalogNumber),
                      data=final_dataset,family="nbinom")

summary(model_draft)


model_draft<-glmer.nb(psite_count~CI*before_after*Life_History+Fish_sp.x*scale(TotalLength_mm)+
                        (1|fish_psite_combo)+(1|CatalogNumber),
                      data=final_dataset,family="nbinom")
summary(model_draft)

big_predictions<-ggeffect(model_draft,c("before_after", "Life_History"))

big_predictions$x

big_plot<-ggplot(big_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=big_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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




model_draft_2<-glmer.nb(psite_count~CI*before_after+Fish_sp.x*scale(TotalLength_mm)+
                        (1|fish_psite_combo)+(1|CatalogNumber),
                      data=final_dataset,family="nbinom")
summary(model_draft_2)

big_predictions<-ggeffect(model_draft_2,c("before_after", "CI"))

big_predictions$x

big_plot<-ggplot(big_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=big_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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



model_draft_3<-glmer.nb(psite_count~CI*before_after+Fish_sp.x*scale(TotalLength_mm)+
                          (before_after|fish_psite_combo)+(1|CatalogNumber),
                        data=final_dataset,family="nbinom")
summary(model_draft_3)

big_predictions<-ggeffect(model_draft_3,c("before_after", "CI"))

big_predictions$x

big_plot<-ggplot(big_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=big_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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


# Extract the fixed and random effects

a<-fixef(model_draft_3)
b<-ranef(model_draft_3,condVar=TRUE)


# Extract the variances of the random effects

qq<-attr(b[[2]],"postVar")
e<-sqrt(qq)
e<-e[2,2,]


# Calculate the CIs

liminf=(b[[2]][2]+a[2])-(e*2)
mean_=(b[[2]][2]+a[2])
limsup=(b[[2]][2]+a[2])+(e*2)

dotchart(mean_[,1],labels=rownames(mean_),cex=0.5)

# add CIs

for(i in 1:nrow(mean_)){
  lines(x=c(liminf[i,1],limsup[i,1]),y=
          c(i,i))
}

raneff_data<-cbind.data.frame(rownames(mean_),mean_,e,(e^2),mean_-(e*2),mean_+(e*2),row.names=NULL)
names(raneff_data)<-c("psite_taxon","interaction","sd","var","min","max")

# Visualize

ranef_plot<-ggplot(raneff_data,aes(psite_taxon,interaction))+
  geom_point(size=3)+
  geom_errorbar(data=raneff_data,mapping=aes(ymin=min,ymax=max),width=0.5)+
  geom_hline(yintercept = -0.3334432, linetype = "dotted")+
  geom_hline(yintercept = 0, linetype = "solid")+
  coord_flip()+
  xlab("parasite taxon")+
  ylab("change over time")+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=14),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_y_reverse()
                   #labels=c("CHONAR"=expression(paste(italic("Chondracanthus narium")))))
#theme(legend.position="none")
ranef_plot


# check the data cleaning script to ensure that all psite groupings are valid
# pick out the sig interactions for the individual parasite taxa and plot them
# maybe the narrative is "no overall change, big change for some psites, variable directions"
# add an analysis of inflection point?


model_draft_4<-glmer.nb(psite_count~CI*before_after+Fish_sp.x*scale(TotalLength_mm)+
                          (CI*before_after|fish_psite_combo)+(1|CatalogNumber),
                        data=final_dataset,family="nbinom")
summary(model_draft_4)

model_draft_5<-glmer.nb(psite_count~CI*before_after+Fish_sp.x*scale(TotalLength_mm)+
                          (CI*before_after|fish_psite_combo)+(1|YearCollected)+(1|CatalogNumber),
                        data=final_dataset,family="nbinom")
summary(model_draft_4)

big_predictions<-ggeffect(model_draft_4,c("before_after", "CI"))

big_predictions$x

big_plot<-ggplot(big_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=big_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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


# Extract the fixed and random effects

a<-fixef(model_draft_4)
b<-ranef(model_draft_4,condVar=TRUE)


# Extract the variances of the random effects

qq<-attr(b[[2]],"postVar")
e<-sqrt(qq)
e<-e[4,2,]


# Calculate the CIs

liminf=(b[[2]][4]+a[11])-(e*2)
mean_=(b[[2]][4]+a[11])
limsup=(b[[2]][4]+a[11])+(e*2)

dotchart(mean_[,1],labels=rownames(mean_),cex=0.5)

# add CIs

for(i in 1:nrow(mean_)){
  lines(x=c(liminf[i,1],limsup[i,1]),y=
          c(i,i))
}

raneff_data<-cbind.data.frame(rownames(mean_),mean_,e,(e^2),mean_-(e*2),mean_+(e*2),row.names=NULL)
names(raneff_data)<-c("psite_taxon","interaction","sd","var","min","max")

# Visualize

ranef_plot<-ggplot(raneff_data,aes(psite_taxon,interaction))+
  geom_point(size=3)+
  geom_errorbar(data=raneff_data,mapping=aes(ymin=min,ymax=max),width=0.5)+
  geom_hline(yintercept = -0.52998176, linetype = "dotted")+
  geom_hline(yintercept = 0, linetype = "solid")+
  coord_flip()+
  xlab("parasite taxon")+
  ylab("change over time")+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=14),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_y_reverse()
#labels=c("CHONAR"=expression(paste(italic("Chondracanthus narium")))))
#theme(legend.position="none")
ranef_plot


### Make individual prediction plots for each fish_psite_combo

raneff_predictions<-ggpredict(model_draft_4,c("before_after", "CI", "fish_psite_combo"))

str(raneff_predictions)
raneff_predictions$facet

raneff_plots<-ggplot(raneff_predictions,aes(x,predicted),group=group,color=group)+
  facet_wrap(vars(facet),nrow=9,ncol=4)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=raneff_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=(rev(levels(raneff_predictions$x))))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
raneff_plots





# Break it up into individual host species

# Check to see whether you need to drop low-prevalence parasites.
# Individual host species analyses may be under-powered.

levels(as.factor(minus_myxos$Fish_sp.x))

car_vel_data <- minus_myxos %>%
  filter(Fish_sp.x=="Carpiodes velifer")
levels(as.factor(car_vel_data$Fish_sp.x))

gam_aff_data <- minus_myxos %>%
  filter(Fish_sp.x=="Gambusia affinis")
levels(as.factor(gam_aff_data$Fish_sp.x))

hyb_nuc_data <- minus_myxos %>%
  filter(Fish_sp.x=="Hybognathus nuchalis")
levels(as.factor(hyb_nuc_data$Fish_sp.x))

ict_pun_data <- minus_myxos %>%
  filter(Fish_sp.x=="Ictalurus punctatus")
levels(as.factor(ict_pun_data$Fish_sp.x))

not_ath_data <- minus_myxos %>%
  filter(Fish_sp.x=="Notropis atherinoides")
levels(as.factor(not_ath_data$Fish_sp.x))

per_vig_data <- minus_myxos %>%
  filter(Fish_sp.x=="Percina vigil")
levels(as.factor(per_vig_data$Fish_sp.x))

pim_vig_data <- minus_myxos %>%
  filter(Fish_sp.x=="Pimephales vigilax")
levels(as.factor(pim_vig_data$Fish_sp.x))

levels(as.factor(car_vel_data$fish_psite_combo))


# Need to account for temporal autocorrelation
# Probably need to scale to get models to behave
# Then you might need to back-transform to plot

car_vel_data$YearCollected

car_vel_model<-glmer.nb(psite_count~CI*before_after*Life_History+
                          (TotalLength_mm|fish_psite_combo)+(1|CatalogNumber)+(1|YearCollected),
                        data=car_vel_data,family="nbinom")

car_vel_model_simplified<-glmer.nb(psite_count~CI*before_after+
                                     (scale(TotalLength_mm)|fish_psite_combo)+(1|CatalogNumber)+(1|YearCollected),
                                   data=car_vel_data,family="nbinom")



summary(car_vel_model_simplified)

car_vel_predictions<-ggeffect(car_vel_model,c("before_after", "CI", "Life_History"))

car_vel_predictions_simplified<-ggeffect(car_vel_model_simplified,c("before_after", "CI"))

car_vel_plot<-ggplot(car_vel_predictions_simplified,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=car_vel_predictions_simplified,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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
car_vel_plot




gam_aff_model<-glmer.nb(psite_count~CI*before_after*Life_History+
                          (TotalLength_mm|fish_psite_combo)+(1|CatalogNumber),
                        data=gam_aff_data,family="nbinom")

gam_aff_model_simplified<-glmer.nb(psite_count~CI*before_after+
                                     (scale(TotalLength_mm)|fish_psite_combo)+(1|CatalogNumber)+(1|YearCollected),
                                   data=gam_aff_data,family="nbinom")

summary(gam_aff_model_simplified)

gam_aff_predictions<-ggeffect(gam_aff_model_simplified,c("before_after", "CI"))

car_vel_predictions_simplified<-ggeffect(car_vel_model_simplified,c("before_after", "CI"))

gam_aff_plot<-ggplot(gam_aff_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=gam_aff_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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
gam_aff_plot



hyb_nuc_model_simplified<-glmer.nb(psite_count~CI*before_after+
                                     (scale(TotalLength_mm)|fish_psite_combo)+(1|CatalogNumber)+(1|YearCollected),
                                   data=hyb_nuc_data,family="nbinom")

summary(hyb_nuc_model_simplified)

hyb_nuc_predictions<-ggeffect(hyb_nuc_model_simplified,c("before_after", "CI"))

hyb_nuc_plot<-ggplot(hyb_nuc_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=hyb_nuc_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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
hyb_nuc_plot


ict_pun_model_simplified<-glmer.nb(psite_count~CI*before_after+
                                     (scale(TotalLength_mm)|fish_psite_combo)+(1|CatalogNumber)+(1|YearCollected),
                                   data=ict_pun_data,family="nbinom")

summary(ict_pun_model_simplified)

ict_pun_predictions<-ggeffect(ict_pun_model_simplified,c("before_after", "CI"))

ict_pun_plot<-ggplot(ict_pun_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=ict_pun_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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
ict_pun_plot



not_ath_model_simplified<-glmer.nb(psite_count~CI*before_after+
                                     (scale(TotalLength_mm)|fish_psite_combo)+(1|CatalogNumber)+(1|YearCollected),
                                   data=not_ath_data,family="nbinom")

summary(not_ath_model_simplified)

not_ath_predictions<-ggeffect(not_ath_model_simplified,c("before_after", "CI"))

not_ath_plot<-ggplot(not_ath_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=not_ath_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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
not_ath_plot


per_vig_model_simplified<-glmer.nb(psite_count~CI*before_after+
                                     (scale(TotalLength_mm)|fish_psite_combo)+(1|CatalogNumber)+(1|YearCollected),
                                   data=per_vig_data,family="nbinom")

summary(per_vig_model_simplified)

per_vig_predictions<-ggeffect(per_vig_model_simplified,c("before_after", "CI"))

per_vig_plot<-ggplot(per_vig_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=per_vig_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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
per_vig_plot



pim_vig_model_simplified<-glmer.nb(psite_count~CI*before_after+
                                     (scale(TotalLength_mm)|fish_psite_combo)+(1|CatalogNumber)+(1|YearCollected),
                                   data=pim_vig_data,family="nbinom")

summary(pim_vig_model_simplified)

pim_vig_predictions<-ggeffect(pim_vig_model_simplified,c("before_after", "CI"))

pim_vig_plot<-ggplot(pim_vig_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=per_vig_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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
pim_vig_plot



### What if we just focused on time?  I checked life history, and it seems to have almost no influence
### on temporal pattern

car_vel_model_LH<-glmer.nb(psite_count~YearCollected+scale(TotalLength_mm)+
                             (1|fish_psite_combo),
                           data=car_vel_data,family="nbinom")

summary(car_vel_model_LH)

car_vel_LH_predictions<-ggeffect(car_vel_model_LH,c("YearCollected", "Life_History"))

car_vel_LH_plot<-ggplot(car_vel_LH_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=car_vel_LH_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  #scale_x_discrete(limits=(rev(levels(big_predictions$x))))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
car_vel_LH_plot



gam_aff_model_LH<-glmer.nb(psite_count~YearCollected+scale(TotalLength_mm)+(1|fish_psite_combo),
                           data=gam_aff_data,family="nbinom")

gam_aff_model_LH<-glm.nb(psite_count~YearCollected+scale(TotalLength_mm),
                         data=gam_aff_data)

summary(gam_aff_model_LH)

gam_aff_LH_predictions<-ggeffect(gam_aff_model_LH,c("before_after", "Life_History"))

gam_aff_LH_plot<-ggplot(gam_aff_LH_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=gam_aff_LH_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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
gam_aff_LH_plot



hyb_nuc_model_LH<-glmer.nb(psite_count~YearCollected+
                             (scale(TotalLength_mm)|fish_psite_combo)+(1|CatalogNumber),
                           data=hyb_nuc_data,family="nbinom")

summary(hyb_nuc_model_LH)

hyb_nuc_LH_predictions<-ggeffect(hyb_nuc_model_LH,c("YearCollected", "Life_History"))

hyb_nuc_LH_plot<-ggplot(hyb_nuc_LH_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=hyb_nuc_LH_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  #scale_x_discrete(limits=(rev(levels(big_predictions$x))))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
hyb_nuc_LH_plot



ict_pun_model_LH<-glm.nb(psite_count~YearCollected+scale(TotalLength_mm),
                         data=ict_pun_data)

summary(ict_pun_model_LH)

ict_pun_LH_predictions<-ggeffect(ict_pun_model_LH,c("before_after", "Life_History"))

ict_pun_LH_plot<-ggplot(ict_pun_LH_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=ict_pun_LH_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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
ict_pun_LH_plot



not_ath_model_LH<-glmer.nb(psite_count~YearCollected+scale(TotalLength_mm)+
                             (1|fish_psite_combo)+(1|CatalogNumber),
                           data=not_ath_data,family="nbinom")

summary(not_ath_model_LH)

not_ath_LH_predictions<-ggeffect(not_ath_model_LH,c("before_after", "Life_History"))

not_ath_LH_plot<-ggplot(not_ath_LH_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=not_ath_LH_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
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
not_ath_LH_plot



per_vig_model_LH<-glmer.nb(psite_count~YearCollected+scale(TotalLength_mm)+
                             (1|fish_psite_combo),
                           data=per_vig_data,family="nbinom")

summary(per_vig_model_LH)

per_vig_LH_predictions<-ggeffect(per_vig_model_LH,c("YearCollected", "Life_History"))

per_vig_LH_plot<-ggplot(per_vig_LH_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=per_vig_LH_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  #scale_x_discrete(limits=(rev(levels(big_predictions$x))))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
per_vig_LH_plot



pim_vig_model_LH<-glm.nb(psite_count~YearCollected+scale(TotalLength_mm),
                         data=pim_vig_data)

summary(pim_vig_model_LH)

pim_vig_LH_predictions<-ggeffect(pim_vig_model_LH,c("YearCollected", "Life_History"))

pim_vig_LH_plot<-ggplot(pim_vig_LH_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=pim_vig_LH_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  #scale_x_discrete(limits=(rev(levels(big_predictions$x))))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
pim_vig_LH_plot
