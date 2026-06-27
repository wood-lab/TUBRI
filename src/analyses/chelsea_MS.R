#### Testing BACI response of parasites to BA*CI (Clean Water Act / pulp mill)
#### Written by Chelsea Wood (chelwood@uw.edu)
#### begun August 2024, completed May 2026

# Load libraries
library(tidyverse)
library(MASS)
library(lme4)
library(emmeans)
library(ggeffects)
library(ggplot2)
 

# Read in the dataset

full_dataset_with_LH<-read.csv("data/processed/Full_dataset_with_psite_life_history_info_2024.11.21.csv")


# Take out myxos, since they were counted differently (you'll run a parallel analysis for them).
minus_myxos <- full_dataset_with_LH %>%
  filter(Parasite_taxonomic_group!="Myxozoa")


# You should also remove any ciliates or nematomorphs - ciliates because we're not sure we always detected
# them, nematomorphs because I'm not convinced it was a nematomorph and because it isn't a fish parasite
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

body_size_fish <- minus_myxos %>%
  group_by(Fish_sp.x) %>%
  summarize(TL = mean(TotalLength_mm), min_TL = min(TotalLength_mm), max_TL = max(TotalLength_mm))

sum(total_n_fish$count)

sum(minus_myxos$psite_count,na.rm=TRUE)

psite_species <- levels(as.factor(minus_myxos$fish_psite_combo))
length(psite_species)




# Now lets remove any parasites that occur at <5% prevalence in their host species.

minus_myxos$positive<-vector("numeric",length(minus_myxos$fish_psite_combo))

for(i in 1:length(minus_myxos$fish_psite_combo)) {
  
  if(is.na(minus_myxos$psite_count[i])){
    minus_myxos$positive[i] <- NA
  } else {
    
    if(minus_myxos$psite_count[i] > 0){
      minus_myxos$positive[i] <- 1
      
    } else {
      
      if(minus_myxos$psite_count[i] == 0){
        minus_myxos$positive[i] <- 0
      }
    }
  }
}

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


# Investigate why there are some repeated prevalences

check <- minus_myxos %>%
  filter(fish_psite_combo == "Hybognathus nuchalis_TREM.ASS" |
           fish_psite_combo == "Hybognathus nuchalis_TREM.BC")
View(check)

check <- minus_myxos %>%
  filter(fish_psite_combo == "Hybognathus nuchalis_TREM.META.GO" |
           fish_psite_combo == "Hybognathus nuchalis_TREM.META.HET")
View(check)

check <- minus_myxos %>%
  filter(fish_psite_combo == "Ictalurus punctatus_NEM.PHAR" |
           fish_psite_combo == "Ictalurus punctatus_NEM.SP")
View(check)

check <- minus_myxos %>%
  filter(fish_psite_combo == "Pimephales vigilax_TREM.BUC" |
           fish_psite_combo == "Pimephales vigilax_TREM.SSS")
View(check)

# All looks fine! These repetitions are not errors. 
# They arise because, within a given fish species, the denominator is the same across all parasite taxa 
# found in that fish species (= to the number of fish individuals examined), so any parasite taxa found 
# in the same number of fish individuals (i.e., with the same numerator) will exhibit the same prevalence. 
# For example, both Neascus sp. and Trematode sp. 1 were found in 18 of 221 Hybognathus nuchalis examined 
# (= 0.0814). Table caption updated accordingly.



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


### Test out various models, starting with the three-way interaction between BA*CI*life_history

model_draft<-glmer.nb(psite_count~CI*before_after*Life_History+Fish_sp.x*scale(TotalLength_mm)+
                        (1|fish_psite_combo)+(1|CatalogNumber),
                      data=final_dataset,family="nbinom")
summary(model_draft)

# Three-way interaction is not significant, but the two-way interaction between BA*life_history is marginal. 
# Let's have a look at what's up.

big_predictions<-ggeffect(model_draft,c("before_after", "Life_History"))

big_predictions$x

big_plot<-ggplot(big_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=big_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("")+
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


# Since three-way interaction isn't significant, drop it and stick with the core hypothesis BA*CI.

model_draft_2<-glmer.nb(psite_count~CI*before_after+Fish_sp.x*scale(TotalLength_mm)+
                        (1|fish_psite_combo)+(1|CatalogNumber),
                      data=final_dataset,family="nbinom")
summary(model_draft_2)

# Two-way interaction BA*CI not significant.

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


# That's fine, but without a random slope for each fish_psite_combo, we only know the net effect - not the
# specific effect on each psite taxon.

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


# Okay, but (before_after|fish_psite_combo) is just the slope of BA for each fish_psite combo. We want the
# slope of BA*CI for each fish_psite combo.

model_draft_4<-glmer.nb(psite_count~CI*before_after+Fish_sp.x*scale(TotalLength_mm)+
                          (CI*before_after|fish_psite_combo)+(1|CatalogNumber),
                        data=final_dataset,family="nbinom")
summary(model_draft_4)


# Sure, but you forgot to account for temporal autocorrelation.

model_draft_5<-glmer.nb(psite_count~CI*before_after+Fish_sp.x*scale(TotalLength_mm)+
                          (CI*before_after|fish_psite_combo)+(1|YearCollected)+(1|CatalogNumber),
                        data=final_dataset,family="nbinom")
summary(model_draft_5)


# Finally: you don't just want the slope of the BA*CI for each fish_psite_combo - you also want to calculate
# a random intercept for each fish_psite_combo, because then each fish_psite_combo can take a different
# baseline abundance.

model_draft_5b<-glmer.nb(psite_count~CI*before_after+Fish_sp.x*scale(TotalLength_mm)+
                          (1+CI*before_after|fish_psite_combo)+(1|YearCollected)+(1|CatalogNumber),
                        data=final_dataset,family="nbinom")
summary(model_draft_5b)


# Still getting some warning messages - try scaling all numerical values, including year. And make sure that
# CatalogNumber is getting treated as a categorical variable.

final_dataset$YearCollected_sc <- scale(final_dataset$YearCollected)
final_dataset$CatalogNumber_chr <- as.character(final_dataset$CatalogNumber)

model_draft_5c<-glmer.nb(psite_count~CI*before_after+Fish_sp.x*scale(TotalLength_mm)+
                           (1+CI*before_after|fish_psite_combo)+(1|YearCollected_sc)+
                           (1|CatalogNumber_chr),
                         data=final_dataset,family="nbinom")
summary(model_draft_5c)

big_predictions<-ggeffect(model_draft_5c,c("before_after", "CI"))

big_plot<-ggplot(big_predictions,aes(x,predicted),group=group,color=group)+
  geom_point(aes(group=group,color=group),size=4,pch=19,position = position_dodge(width = 0.5))+
  geom_errorbar(data=big_predictions,mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),
                width=0.04,position = position_dodge(width = 0.5))+
  geom_line(aes(group=group,color=group),position = position_dodge(width = 0.5))+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=(rev(levels(big_predictions$x))))+
  theme(legend.position="top",legend.title = element_blank(),
        legend.text = element_text(size=12))
big_plot


# Extract the fixed and random effects

a<-fixef(model_draft_5c)
b<-ranef(model_draft_5c,condVar=TRUE)


# Extract the variances of the random effects

qq<-attr(b[[3]],"postVar")
e<-sqrt(qq)
e<-e[4,2,]


# Calculate the CIs

liminf=(b[[3]][4]+a[11])-(e*2)
mean_=(b[[3]][4]+a[11])
limsup=(b[[3]][4]+a[11])+(e*2)

dotchart(mean_[,1],labels=rownames(mean_),cex=0.5)

# add CIs

for(i in 1:nrow(mean_)){
  lines(x=c(liminf[i,1],limsup[i,1]),y=
          c(i,i))
}

raneff_data<-cbind.data.frame(rownames(mean_),mean_,e,(e^2),mean_-(e*2),mean_+(e*2),row.names=NULL)
names(raneff_data)<-c("psite_taxon","interaction","sd","var","min","max")
psite_code_break<-read.csv("data/translating_psite_codes.csv", sep=",", header=TRUE)
raneff_data_named<-merge(raneff_data, psite_code_break, by.x = "psite_taxon", by.y = "data_combo")

# Visualize

library(forcats)

ranef_plot<-ggplot(raneff_data_named,aes(x = fct_reorder(psite_uniform_code, interaction),interaction))+
  geom_point(size=3)+
  geom_errorbar(data=raneff_data_named,mapping=aes(ymin=min,ymax=max),width=0.5)+
  geom_hline(yintercept = -0.52998176, linetype = "dotted")+
  geom_hline(yintercept = 0, linetype = "solid")+
  coord_flip()+
  xlab("parasite taxon")+
  ylab("interaction effect BA*CI")+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),
        axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=14),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
#labels=c("CHONAR"=expression(paste(italic("Chondracanthus narium")))))
#theme(legend.position="none")
ranef_plot


### Make individual prediction plots for each fish_psite_combo

raneff_predictions<-ggpredict(model_draft_5c,c("before_after", "CI", "fish_psite_combo"))

str(raneff_predictions)
raneff_predictions$facet

filtered_predictions <- raneff_predictions %>%
  filter(facet=="Pimephales vigilax_TREM.SSS" | facet == "Pimephales vigilax_TREM.METAGS" |
           facet == "Pimephales vigilax_TREM.BUC" | facet == "Hybognathus nuchalis_TREM.BC" |
           facet == "Pimephales vigilax_META.UNK" | facet == "Percina vigil_NEM.LARV" |
           facet == "Carpiodes velifer_TREM.META.UNK" | facet == "Ictalurus punctatus_TREM.LG" |
           facet == "Hybognathus nuchalis_TREM.META.GO" | facet == "Pimephales vigilax_MONO.DACT" |
           facet == "Hybognathus nuchalis_TREM.META.HET" | facet == "Notropis atherinoides_TREM.META" |
           facet == "Notropis atherinoides_TREM.LARV" | facet == "Ictalurus punctatus_MONO.IP" |
           facet == "Gambusia affinis_TREM.POS" | facet == "Notropis atherinoides_MONO.ALL")

raneff_plots<-ggplot(filtered_predictions,aes(x,predicted),group=group,color=group)+
  facet_wrap(vars(facet),nrow=4,ncol=4)+
  geom_point(aes(group=group,color=group),size=4,pch=19)+
  geom_errorbar(data=filtered_predictions,
                mapping=aes(x=x,ymin=conf.low,ymax=conf.high,group=group,color=group),width=0.03)+
  geom_line(aes(group=group,color=group))+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  xlab("year")+
  ylab("predicted parasite abundance\n per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),
        axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=(rev(levels(filtered_predictions$x))))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
raneff_plots


# Okay, so model_draft_5c is the final version of the BA*CI model. But it would be cool if we could turn
# BA into a continuous variable... see if the models will run?

model_draft_6<-glmer.nb(psite_count~CI*YearCollected_sc+Fish_sp.x*scale(TotalLength_mm)+
                          (1+CI*YearCollected_sc|fish_psite_combo)+(1|CatalogNumber_chr),
                        data=final_dataset,family="nbinom")
summary(model_draft_6)

big_predictions<-ggeffect(model_draft_6,c("YearCollected_sc", "CI"))

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

a<-fixef(model_draft_6)
b<-ranef(model_draft_6,condVar=TRUE)


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

# Whole lotta nuthin'. Not very informative. Would probably be better just to plot the raw data across time, 
# or to do the changepoint analysis originally planned if you want to explore the granular temporal trend.



### Okay now time to look at the raw data - first in BA form, then in granular temporal gradient, with
### changepoint analysis.

str(final_dataset)

raw_plot<-ggplot(final_dataset,aes(before_after,psite_count),group=CI,color=CI)+
  #facet_wrap(vars(fish_psite_combo),nrow=9,ncol=4)+
  geom_point(aes(group=CI,color=CI),position=position_jitter(width=0.15),size=3,pch=19)+
  geom_violin(aes(fill=CI))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  #xlab("year")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  ylim(0,1)+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=(rev(levels(raneff_predictions$x))))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
raw_plot

# Not very useful. Maybe we can aggregate this a bit to make the pattern easier to see, 
# averaging within individual fish?

within_fish <- final_dataset %>%
  group_by(CI,before_after,IndividualFishID) %>%
  summarize(mean_psite_count = mean(psite_count),summed_psite_count = sum(psite_count))

raw_plot<-ggplot(within_fish,aes(before_after,mean_psite_count,group=str_c(CI,before_after)))+
  #facet_wrap(vars(fish_psite_combo),nrow=9,ncol=4)+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("mean parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  #ylim(0,1)+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=(rev(levels(as.factor(within_fish$before_after)))))+
  theme(legend.position="top",legend.title = element_blank(),
        legend.text = element_text(size=12))
raw_plot

truncated_raw_plot<-ggplot(within_fish,aes(before_after,mean_psite_count,group=str_c(CI,before_after)))+
  #facet_wrap(vars(fish_psite_combo),nrow=9,ncol=4)+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("mean parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=(rev(levels(as.factor(within_fish$before_after)))))+
  theme(legend.position="top",legend.title = element_blank(),
        legend.text = element_text(size=12))+
  coord_cartesian(ylim = c(0, 10))
truncated_raw_plot


# Gorgeous. Now put it all together into one figure.
library(cowplot)
main_results_figure <- ggdraw(plot=NULL,xlim=c(0,15),ylim=c(0,5))+
  draw_plot(raw_plot,x=0,y=0.25,width=4.75,height=4.75)+
  draw_plot(truncated_raw_plot,x=5,y=0.25,width=4.75,height=4.75)+
  draw_plot(big_plot,x=10,y=0.25,width=4.75,height=4.75)+
  draw_label("time period (relative to 1973 Clean Water Act)",x=7.75,y=0.25,size=18)+
  draw_label("(a)",x=0.9,y=4.75,size=25)+
  draw_label("(b)",x=5.9,y=4.75,size=25)+
  draw_label("(c)",x=10.9,y=4.75,size=25)
main_results_figure


# Now do the same thing for individual parasite taxa

raw_plot_per_psite<-ggplot(final_dataset,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  facet_wrap(vars(fish_psite_combo),nrow=8,ncol=4)+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("parasite abundance\n per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=(rev(levels(as.factor(within_fish$before_after)))))+
  theme(legend.position="top",legend.title = element_blank(),
        legend.text = element_text(size=12))+
  coord_cartesian(ylim = c(0, 1))
raw_plot_per_psite


# Okay - a lot going on, every parasite has its own story. I think we need to focus just on those parasites
# with significant interaction effects, and show what they're up to. As in the other main figure, this
# needs to be the overall raw data and the cartesian version of the plot (zoomed in on where the action is).
# All we need to do is isolate each parasite and then plot it.

psite_code_break<-read.csv("data/translating_psite_codes.csv", sep=",", header=TRUE)
final_data_named<-merge(final_dataset, psite_code_break, by.x = "fish_psite_combo", by.y = "data_combo")


PIMVIGTREMSP3_data <- final_data_named %>%
  filter(psite_uniform_code == "PIMVIG-TREM-SP3")

PIMVIGTREMSP3_plot<-ggplot(PIMVIGTREMSP3_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("PIMVIG-TREM-SP3")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
PIMVIGTREMSP3_plot

PIMVIGTREMSP3_plot_trunc<-ggplot(PIMVIGTREMSP3_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("PIMVIG-TREM-SP3")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 1))
PIMVIGTREMSP3_plot_trunc


PIMVIGTREMBU2_data <- final_data_named %>%
  filter(psite_uniform_code == "PIMVIG-TREM-BU2")

PIMVIGTREMBU2_plot<-ggplot(PIMVIGTREMBU2_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("PIMVIG-TREM-BU2")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 1))+
  theme(legend.position = "none")
PIMVIGTREMBU2_plot

PIMVIGTREMBU2_plot_trunc<-ggplot(PIMVIGTREMBU2_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("PIMVIG-TREM-BU2")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 20))
PIMVIGTREMBU2_plot_trunc


PIMVIGTREMBU1_data <- final_data_named %>%
  filter(psite_uniform_code == "PIMVIG-TREM-BU1")

PIMVIGTREMBU1_plot<-ggplot(PIMVIGTREMBU1_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("PIMVIG-TREM-BU1")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
PIMVIGTREMBU1_plot

PIMVIGTREMBU1_plot_trunc<-ggplot(PIMVIGTREMBU1_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("PIMVIG-TREM-BU1")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 1))
PIMVIGTREMBU1_plot_trunc


HYBNUCTREMNEA_data <- final_data_named %>%
  filter(psite_uniform_code == "HYBNUC-TREM-NEA")

HYBNUCTREMNEA_plot<-ggplot(HYBNUCTREMNEA_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("HYBNUC-TREM-NEA")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
HYBNUCTREMNEA_plot

HYBNUCTREMNEA_plot_trunc<-ggplot(HYBNUCTREMNEA_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("HYBNUC-TREM-NEA")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 2))
HYBNUCTREMNEA_plot_trunc


PIMVIGTREMSP1_data <- final_data_named %>%
  filter(psite_uniform_code == "PIMVIG-TREM-SP1")

PIMVIGTREMSP1_plot<-ggplot(PIMVIGTREMSP1_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("PIMVIG-TREM-SP1")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
PIMVIGTREMSP1_plot

PIMVIGTREMSP1_plot_trunc<-ggplot(PIMVIGTREMSP1_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("PIMVIG-TREM-SP1")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 2))
PIMVIGTREMSP1_plot_trunc


PERVIGNEMAASC_data <- final_data_named %>%
  filter(psite_uniform_code == "PERVIG-NEMA-ASC")

PERVIGNEMAASC_plot<-ggplot(PERVIGNEMAASC_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("PERVIG-NEMA-ASC")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
PERVIGNEMAASC_plot

PERVIGNEMAASC_plot_trunc<-ggplot(PERVIGNEMAASC_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("PERVIG-NEMA-ASC")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 0.5))
PERVIGNEMAASC_plot_trunc


CARVELTREMSP1_data <- final_data_named %>%
  filter(psite_uniform_code == "CARVEL-TREM-SP1")

CARVELTREMSP1_plot<-ggplot(CARVELTREMSP1_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("CARVEL-TREM-SP1")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
CARVELTREMSP1_plot

CARVELTREMSP1_plot_trunc<-ggplot(CARVELTREMSP1_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("CARVEL-TREM-SP1")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 0.5))
CARVELTREMSP1_plot_trunc


ICTPUNTREMSP1_data <- final_data_named %>%
  filter(psite_uniform_code == "ICTPUN-TREM-SP1")

ICTPUNTREMSP1_plot<-ggplot(ICTPUNTREMSP1_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("ICTPUN-TREM-SP1")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
ICTPUNTREMSP1_plot

ICTPUNTREMSP1_plot_trunc<-ggplot(ICTPUNTREMSP1_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("ICTPUN-TREM-SP1")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 5))
ICTPUNTREMSP1_plot_trunc


HYBNUCTREMSP2_data <- final_data_named %>%
  filter(psite_uniform_code == "HYBNUC-TREM-SP2")

HYBNUCTREMSP2_plot<-ggplot(HYBNUCTREMSP2_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("HYBNUC-TREM-SP2")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
HYBNUCTREMSP2_plot

HYBNUCTREMSP2_plot_trunc<-ggplot(HYBNUCTREMSP2_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("HYBNUC-TREM-SP2")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 1))
HYBNUCTREMSP2_plot_trunc


PIMVIGMONODAC_data <- final_data_named %>%
  filter(psite_uniform_code == "PIMVIG-MONO-DAC")

PIMVIGMONODAC_plot<-ggplot(PIMVIGMONODAC_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("PIMVIG-MONO-DAC")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
PIMVIGMONODAC_plot

PIMVIGMONODAC_plot_trunc<-ggplot(PIMVIGMONODAC_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("PIMVIG-MONO-DAC")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 4))
PIMVIGMONODAC_plot_trunc


HYBNUCTREMSP3_data <- final_data_named %>%
  filter(psite_uniform_code == "HYBNUC-TREM-SP3")

HYBNUCTREMSP3_plot<-ggplot(HYBNUCTREMSP3_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("HYBNUC-TREM-SP3")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
HYBNUCTREMSP3_plot

HYBNUCTREMSP3_plot_trunc<-ggplot(HYBNUCTREMSP3_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("HYBNUC-TREM-SP3")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 0.5))
HYBNUCTREMSP3_plot_trunc


NOTATHTREMSP2_data <- final_data_named %>%
  filter(psite_uniform_code == "NOTATH-TREM-SP2")

NOTATHTREMSP2_plot<-ggplot(NOTATHTREMSP2_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("NOTATH-TREM-SP2")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
NOTATHTREMSP2_plot

NOTATHTREMSP2_plot_trunc<-ggplot(NOTATHTREMSP2_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("NOTATH-TREM-SP2")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 5))
NOTATHTREMSP2_plot_trunc


NOTATHTREMHET_data <- final_data_named %>%
  filter(psite_uniform_code == "NOTATH-TREM-HET")

NOTATHTREMHET_plot<-ggplot(NOTATHTREMHET_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("NOTATH-TREM-HET")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
NOTATHTREMHET_plot

NOTATHTREMHET_plot_trunc<-ggplot(NOTATHTREMHET_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("NOTATH-TREM-HET")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 2))
NOTATHTREMHET_plot_trunc


ICTPUNMONOLIG_data <- final_data_named %>%
  filter(psite_uniform_code == "ICTPUN-MONO-LIG")

ICTPUNMONOLIG_plot<-ggplot(ICTPUNMONOLIG_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("ICTPUN-MONO-LIG")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
ICTPUNMONOLIG_plot

ICTPUNMONOLIG_plot_trunc<-ggplot(ICTPUNMONOLIG_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("ICTPUN-MONO-LIG")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 10))
ICTPUNMONOLIG_plot_trunc


GAMAFFTREMPOS_data <- final_data_named %>%
  filter(psite_uniform_code == "GAMAFF-TREM-POS")

GAMAFFTREMPOS_plot<-ggplot(GAMAFFTREMPOS_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("GAMAFF-TREM-POS")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
GAMAFFTREMPOS_plot

GAMAFFTREMPOS_plot_trunc<-ggplot(GAMAFFTREMPOS_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("GAMAFF-TREM-POS")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 8))
GAMAFFTREMPOS_plot_trunc


NOTATHMONODAC_data <- final_data_named %>%
  filter(psite_uniform_code == "NOTATH-MONO-DAC")

NOTATHMONODAC_plot<-ggplot(NOTATHMONODAC_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("NOTATH-MONO-DAC")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")
NOTATHMONODAC_plot

NOTATHMONODAC_plot_trunc<-ggplot(NOTATHMONODAC_data,aes(before_after,psite_count,group=str_c(CI,before_after)))+
  geom_point(aes(color=str_c(CI)),position = position_jitterdodge(seed = 1, dodge.width = 0.9, 
                                                                  jitter.width = 0.8),
             size=1,pch=19)+
  geom_violin(aes(fill=str_c(CI)),position = position_dodge(width = 0.9),alpha=0.5)+
  stat_summary(
    fun = mean, 
    geom = "point", 
    aes(color = CI),
    size=5, 
    position = position_dodge(width = 0.9),
    alpha = 0.75
  )+
  stat_summary(
    fun = "mean", 
    geom = "line",
    aes(group = CI, color = CI),
    position = position_dodge(width = 0.9),
    linetype = "solid", 
    linewidth = 0.5)+
  scale_color_manual(values=c("cadetblue","burlywood4"))+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  ggtitle("NOTATH-MONO-DAC")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=12,hjust=0.5,face="plain"),axis.text.y=element_text(size=10),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(labels=NULL,limits=(rev(levels(as.factor(within_fish$before_after)))))+
  scale_y_continuous(labels = label_number(accuracy = 0.1))+
  theme(legend.position = "none")+
  coord_cartesian(ylim = c(0, 8))
NOTATHMONODAC_plot_trunc


# Now put it all together

psite_taxon_raw_together <- ggdraw(plot=NULL,xlim=c(0,8),ylim=c(0,16))+
  draw_plot(PIMVIGTREMSP3_plot_trunc,x=0,y=14,width=4,height=2)+
  draw_plot(PIMVIGTREMBU2_plot_trunc,x=0,y=12,width=4,height=2)+
  draw_plot(PIMVIGTREMBU1_plot_trunc,x=0,y=10,width=4,height=2)+
  draw_plot(HYBNUCTREMNEA_plot_trunc,x=0,y=8,width=4,height=2)+
  draw_plot(PIMVIGTREMSP1_plot_trunc,x=0,y=6,width=4,height=2)+
  draw_plot(PERVIGNEMAASC_plot_trunc,x=0,y=4,width=4,height=2)+
  draw_plot(CARVELTREMSP1_plot_trunc,x=0,y=2,width=4,height=2)+
  draw_plot(ICTPUNTREMSP1_plot_trunc,x=0,y=0,width=4,height=2)+
  draw_plot(HYBNUCTREMSP2_plot_trunc,x=4,y=14,width=4,height=2)+
  draw_plot(PIMVIGMONODAC_plot_trunc,x=4,y=12,width=4,height=2)+
  draw_plot(HYBNUCTREMSP3_plot_trunc,x=4,y=10,width=4,height=2)+
  draw_plot(NOTATHTREMSP2_plot_trunc,x=4,y=8,width=4,height=2)+
  draw_plot(NOTATHTREMHET_plot_trunc,x=4,y=6,width=4,height=2)+
  draw_plot(ICTPUNMONOLIG_plot_trunc,x=4,y=4,width=4,height=2)+
  draw_plot(GAMAFFTREMPOS_plot_trunc,x=4,y=2,width=4,height=2)+
  draw_plot(NOTATHMONODAC_plot_trunc,x=4,y=0,width=4,height=2)
#draw_label("(a)",x=0.25,y=4.75,size=30)+
#draw_label("(b)",x=5,y=4.75,size=30)+
#draw_label("(c)",x=10.4,y=4.75,size=30)
psite_taxon_raw_together

# Make titles smaller
# Add annotations for x- and y-axis labels
# Minimize white space
# Create an identical set with full axis for supplement
# For big raneff plot, reverse the order of the y-axis to match
# Update plots with NEMA (and caption)




raw_plots_per_psite<-ggplot(final_dataset,aes(before_after,psite_count),group=CI,color=CI)+
  facet_wrap(vars(fish_psite_combo),nrow=9,ncol=4)+
  #geom_point(aes(group=CI,color=CI),position=position_jitter(width=0.15),size=4,pch=19)+
  geom_violin(aes(fill=CI))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  ylim(0,10)+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=(rev(levels(raneff_predictions$x))))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
raw_plots_per_psite


# Now over time

raw_plot_time<-ggplot(final_dataset,aes(YearCollected,psite_count),group=CI,color=CI)+
  #facet_wrap(vars(fish_psite_combo),nrow=9,ncol=4)+
  geom_point(aes(group=CI,color=CI),position=position_jitter(width=0.15),size=3,pch=19)+
  geom_smooth(aes(color=CI))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  ylim(0,50)+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  #scale_x_discrete(limits=(rev(levels(raneff_predictions$x))))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
raw_plot_time


raw_plot_time_per_psite<-ggplot(final_dataset,aes(YearCollected,psite_count),group=CI,color=CI)+
  facet_wrap(vars(fish_psite_combo),nrow=9,ncol=4)+
  geom_point(aes(group=CI,color=CI),position=position_jitter(width=0.15),size=3,pch=19)+
  #geom_smooth(aes(color=CI))+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  ylim(0,10)+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  #scale_x_discrete(limits=(rev(levels(raneff_predictions$x))))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
raw_plot_time_per_psite


raw_plot_body_size<-ggplot(final_dataset,aes(YearCollected,TotalLength_mm),group=CI,color=CI)+
  facet_wrap(vars(Fish_sp.x),nrow=2,ncol=4)+
  geom_point(aes(group=CI,color=CI),position=position_jitter(width=0.15),size=3,pch=19)+
  geom_smooth(aes(color=CI),method="lm")+
  #scale_color_manual(name = c(""),values=plasma_pal)+
  xlab("year")+
  ylab("total length of host fish (mm)")+
  theme_minimal()+
  #ylim(0,10)+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=12,color="black"),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),
        panel.grid.minor=element_line(color=NA),panel.spacing = unit(1, "lines"))+
  #scale_x_discrete(limits=(rev(levels(raneff_predictions$x))))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
raw_plot_body_size




### Make a nice plot of all lots

library(maps)
library(mapdata)
#library(maptools) #for shapefiles
library(scales) #for transparency
#library(rgdal)
library(maps)
library(mapdata)
library(ggmap)
#library(ggsn)
library(tidyverse)
library(tidyr)
library(dplyr)
library(cowplot)
library(sp)
#library(rgeos)
library(raster)
library(magick)
library(pdftools)


# First you need to reduce the dataset to make each fish a row

fish_sampled <- final_dataset %>%
  group_by(Fish_sp.x,IndividualFishID,CI,Latitude,Longitude,YearCollected) %>%
  summarize(count = n())


# Make the map

bounds<-c(left=-89.855, bottom=30.69, right=-89.81, top=30.81)
map<-get_stadiamap(bounds, zoom=12, maptype = "stamen_terrain_background") %>% ggmap()+
  geom_point(data = fish_sampled, aes(x=Longitude,y=Latitude,fill=CI),shape=21,size=4)+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  geom_point(y=30.76558198250576, x=-89.83398463056083, size=4)+
  xlab("")+
  ylab("")+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=15,angle=45,hjust=1),axis.text.y=element_text(size=20))
#geom_rect(xmin=-106.89,ymin=34.8,xmax=-106.45,ymax=35.2102673,fill="grey",alpha=0.2,linetype="blank")+
#geom_hline(yintercept=35.2102673,linetype="dashed")+
#geom_map(data=water_AL_df,map=water_AL_df,aes(x=long,y=lat,map_id=id),color="lightsteelblue3",fill="lightsteelblue3")
#geom_map(data=polygon_df,map=polygon_df,aes(x=long,y=lat,map_id=id),color="black",fill=NA)+
#annotate("text",x=-106.6,y=35.097,label="City of Albuquerque",size=5)+
#annotate("text",x=-106.58,y=35.33,label="Rio Grande River",size=5,angle=55)
map

# Make the inset map

bounds<-c(left=-124.77, bottom=24.52, right=-66.95, top=49.38)
inset_map_us<-get_stadiamap(bounds, zoom=5, maptype = "stamen_terrain_background") %>% ggmap()+
  xlab("")+
  ylab("")+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=15,angle=45,hjust=1),axis.text.y=element_text(size=20))+
  geom_rect(xmin=-92,ymin=28.922129,xmax=-88.821197,ymax=31,fill="black",alpha=0.1,linetype="blank")
#geom_hline(yintercept=35.2102673,linetype="dashed")+
#geom_map(data=water_AL_df,map=water_AL_df,aes(x=long,y=lat,map_id=id),color="lightsteelblue3",fill="lightsteelblue3")
#geom_map(data=polygon_df,map=polygon_df,aes(x=long,y=lat,map_id=id),color="black",fill=NA)+
#annotate("text",x=-106.6,y=35.097,label="City of Albuquerque",size=5)+
#annotate("text",x=-106.58,y=35.33,label="Rio Grande River",size=5,angle=55)
inset_map_us

bounds<-c(left=-92, bottom=28.922129, right=-88.821197, top=31)
inset_map_la<-get_stadiamap(bounds, zoom=7, maptype = "stamen_terrain") %>% ggmap()+
  xlab("")+
  ylab("")+
  theme(legend.position = "none")+
  theme(plot.margin = unit(c(0,0,0,0), "cm"),axis.title.x=element_text(size=20),axis.title.y=element_text(size=20),
        axis.text.x=element_text(size=15,angle=45,hjust=1),axis.text.y=element_text(size=20))+
  geom_rect(xmin=-89.855,ymin=30.69,xmax=-89.81,ymax=30.81,fill="black",alpha=0.2,linetype="blank")
#geom_hline(yintercept=35.2102673,linetype="dashed")+
#geom_map(data=water_AL_df,map=water_AL_df,aes(x=long,y=lat,map_id=id),color="lightsteelblue3",fill="lightsteelblue3")
#geom_map(data=polygon_df,map=polygon_df,aes(x=long,y=lat,map_id=id),color="black",fill=NA)+
#annotate("text",x=-106.6,y=35.097,label="City of Albuquerque",size=5)+
#annotate("text",x=-106.58,y=35.33,label="Rio Grande River",size=5,angle=55)
inset_map_la


# Plot the fish sampled

fish_plot<-ggplot(data=fish_sampled,aes(x=YearCollected,y=CI))+
  scale_y_discrete(limits=rev,labels=c("I","C"))+
  facet_wrap(vars(Fish_sp.x),nrow=4,ncol=2)+
  geom_point(data=fish_sampled,aes(x=YearCollected,fill=CI),
             position=position_jitter(width=0.15),shape=21,size=3)+
  scale_fill_manual(values=c("cadetblue","burlywood4"))+
  xlab("")+
  ylab("")+
  annotate("rect", xmin = -Inf, xmax = 1972+(292/366), ymin = -Inf, ymax = Inf,
           fill = "black", alpha = 0.2)+
  geom_vline(xintercept=1972+(292/366),linetype="dashed")+
  geom_hline(yintercept=1.5,linetype="solid")+
  #Katrina: geom_vline(xintercept=2005,linetype="dotted")+
  xlim(1962,2006)+
  theme_bw()+
  theme(legend.position = "none")+
  guides(fill = guide_legend(override.aes = list(size=7)))+
  scale_x_continuous(breaks=seq(1960,2010,10))+
  theme(plot.margin = unit(c(0,0,0,-0.6), "cm"), text=element_text(family='sans',size=20),
        axis.text.y = element_text(size=14),axis.text.x = element_text(size=14),
        legend.text = element_text(size=18),strip.text = element_text(face="italic"),
        strip.background = element_blank())
fish_plot


# Set the map and the plot side-by-side
library(cowplot)
final_figure <- ggdraw(plot=NULL,xlim=c(0,10),ylim=c(0,10))+
  draw_plot(map,x=0.5,y=-0.1,width=6,height=9.75)+
  draw_plot(inset_map_us,x=-0.8,y=5.5,width=4,height=3)+
  draw_plot(inset_map_la,x=-0.8,y=1.5,width=4,height=3.5)+
  draw_plot(fish_plot,x=4.8,y=0,width=5,height=10)
#draw_label("(a)",x=0.25,y=4.75,size=30)+
#draw_label("(b)",x=5,y=4.75,size=30)+
#draw_label("(c)",x=10.4,y=4.75,size=30)
final_figure




