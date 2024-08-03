# Shyanne's data exploration, plotting, and analysis
# by Shyanne Christner and Chelsea Wood
# christner.shyanne03@gmail.com and chelwood@uw.edu
# Started 4 July 2024


# Load up the packagaes you need

library(ggplot2)
library(tidyverse)
library(lme4)
library(car)
library(ggeffects)
library(reshape2)


# R is really useful for doing quick tallies and plots.
# Let's do some together!

library(tidyverse)

ict_pun_data<-read.csv("data/processed/Ictalurus_punctatus_processed_human_readable_UPDATED_2024.07.10.csv")
pim_vig_data <-read.csv("data/processed/Pimephales_vigilax_processed_human_readable_UPDATED_2024.07.10.csv")
not_ath_data <-read.csv("data/processed/Notropis_atherinoides_processed_human_readable.csv")

# Let's do some quick data tallies

tally_up <- ict_pun_data %>%
  group_by(combo) %>%
  summarize(mean.MONO.IP = mean(MONO.IP))

tally_up <- ict_pun_data %>%
  group_by(YearCollected) %>%
  summarize(mean.MONO.IP = mean(MONO.IP))
View(tally_up)


# Let's do some quick plots

plot(ict_pun_data$MONO.IP~ict_pun_data$Latitude)+abline(v=30.76)


# How about some quick stats?

summary(lm(pim_vig_data$MONO.IP~pim_vig_data$Latitude))




### PRELIMINARY ANALYSIS - CHELSEA, 26 JULY 2024

full_data<-read.csv("data/processed/Full_dataset_with_psite_life_history_info_2024.08.01.csv", header = T, sep = ",")

colnames(full_data)[20]<-"scaled_TL_mm"

View(full_data)


# First create a model where you are looking at Life_History (direct versus complex)
# Note that this model may take a few minutes to run.

trimmed_data <- full_data %>%
  filter(Life_History=="Direct" | Life_History=="Complex")


# Remember that you need to remove myxos from this analysis, because they are quantified differently than the
# other parasites (i.e., by number of cysts or presence/absence, not by number of individuals).

trimmed_data <- trimmed_data %>%
  filter(Parasite_taxonomic_group!="Myxozoa")


# The model won't run with all the data, and I narrowed the culprit down to Hybognathus nuchalis - it doesn't have 
# enough directly transmitted parasites to make a fair comparison between direct versus complex!  So we will have to 
# drop it from this particular analysis.

trimmed_data <- trimmed_data %>%
  filter(Fish_sp.x!="Hybognathus nuchalis")

model_1<-glmer.nb(psite_count~CI*Life_History+(1|Fish_sp.x/psite_spp.x)+
                    offset(log(scaled_TL_mm)),data=trimmed_data,family="nbinom")

summary(model_1)

levels(as.factor(trimmed_data$Fish_sp.x))


# At some point, you should probably control for season
# (1|MonthCollected)


# Let's save the model output down somewhere so that we can grab it again when we need it without having to wait
# for the model to run and converge

library(stargazer)
model_1_output<-stargazer::stargazer(model_1, type = "text")
write.table(model_1_output,"REU/shyanne/model_1_output_2024.08.01")

# Check VIFs to make sure that there is no collinearity
library(car)
vif(model_1)

# Now plot it so that you can visualize the results

library(ggeffects)
color_predict<-ggeffect(model_1,c("Life_History","CI"))

library(ggplot2)
color_plot<-ggplot(color_predict,aes(c(1.95,0.95,2.05,1.05),predicted),group=Life_History,linetype=Life_History)+
  geom_errorbar(data=color_predict,mapping=aes(x=c(1.95,0.95,2.05,1.05),ymin=conf.low,ymax=conf.high),width=0.03)+
  geom_line(aes(group=group,linetype=group))+
  geom_point(size=6,pch=21,fill=c("burlywood4","cadetblue","burlywood4","cadetblue"))+
  xlab("treatment")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  scale_linetype_discrete(name = c("life history"),labels=c("direct","complex"))+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color=c("cadetblue","burlywood4")),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=rev(levels(color_predict$group)),labels=c("control","impact"))+
  annotate("text",label="effect of treatment:\np < 0.0001",x = 1.9, y = 0.65, size = 6)+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=14))
color_plot




# Now create a model where you are looking at Parasite_taxonomic_group (Nematoda, Trematoda, etc)
# Note that this model will take a WHILE to run. It took about 37 minutes on my computer.

# To do this, we need to trim out a lot of things from the dataset, like leeches, nematomorphs, and myxos.

stuff <- trimmed_data %>%
  group_by(Parasite_taxonomic_group) %>%
  summarize(count = n())

# Take out the copepods

trimmed_data <- trimmed_data %>%
  filter(Parasite_taxonomic_group!="Copepoda")

library(lme4)
model_2<-glmer.nb(psite_count~CI*Parasite_taxonomic_group+(1|Fish_sp.x/psite_spp.x)+
                    offset(log(scaled_TL_mm)),data=trimmed_data,family="nbinom")


# At some point, you should probably control for season
# (1|MonthCollected)

summary(model_2)

# Check VIFs to make sure that there is no collinearity
library(car)
vif(model_2)


# Let's save the model output down somewhere so that we can grab it again when we need it without having to wait
# for the model to run and converge

model_2_output<-stargazer::stargazer(model_2, type = "text")
write.table(model_2_output,"REU/shyanne/model_2_output_2024.08.02")


# Now let's plot it

library(ggeffects)
color_predict_2<-ggeffect(model_2,c("Parasite_taxonomic_group","CI"))

# Need to exclude Hirudinea because error bars are wacky (probably because we found so few). For now, I've
# just taken the error bars off.

# Set up a nice color palette first.
library(viridis)
plasma_pal <- c("red", viridis::plasma(n = 5))
pal<-viridis(n=6)

pal.bands(coolwarm)

par(op)
color_predict_2$group

color_plot_2<-ggplot(color_predict_2,aes(c(0.85,0.9,0.95,1,1.05,1.85,1.9,1.95,2,2.05),predicted),
                     group=Parasite_taxonomic_group,color=Parasite_taxonomic_group)+
  geom_point(aes(group=x,color=x),size=4,pch=19)+
  geom_errorbar(data=color_predict_2,mapping=aes(x=c(0.85,0.9,0.95,1,1.05,1.85,1.9,1.95,2,2.05),
                                                 ymin=conf.low,ymax=conf.high,group=x,color=x),width=0.03)+
  geom_line(aes(group=x,color=x))+
  scale_color_manual(name = c("parasite taxonomic group"),values=plasma_pal)+
  xlab("treatment")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color=c("cadetblue","burlywood4")),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=(levels(color_predict_2$group)))+
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=12))
color_plot_2



# Now make a dedicated analysis and plot for Myxozoa

myxo_data <- full_data %>%
  filter(Parasite_taxonomic_group=="Myxozoa")

# For some myxos, we counted cysts, while others were presence (1) / absence (0). We need to turn the whole
# myxo dataset into presence/absence, and the loop below does that.

myxo_data$presence<-vector("character",length(myxo_data$fish_psite_combo))

for(i in 1:length(myxo_data$fish_psite_combo)) {
  
  if(is.na(myxo_data$psite_count[i])){
    myxo_data$presence[i] <- NA
  } else {
    
    if(myxo_data$psite_count[i]>0){
      myxo_data$presence[i] <- 1
      
    } else {
      
      if(myxo_data$psite_count[i]==0){
        myxo_data$presence[i] <- 0
      }
    }
  }
}
myxo_data$presence
myxo_data$presence<-as.numeric(myxo_data$presence)

model_myxo<-glmer(presence~CI+(1|Fish_sp.x/psite_spp.x)+
                    offset(log(scaled_TL_mm)),data=myxo_data,family="binomial")
summary(model_myxo)


# Now let's plot it

library(ggeffects)
myxo_predict<-ggemmeans(model_myxo,"CI")

myxo_predict$x
myxo_predict$predicted

myxo_plot<-ggplot(myxo_predict,aes(x,predicted))+
  geom_errorbar(data=myxo_predict,mapping=aes(x=x,y=predicted,ymin=conf.low,ymax=conf.high),width=0.03)+
  geom_line(aes(group=group,linetype=group))+
  geom_point(size=6,pch=21,fill=c("cadetblue","burlywood4"))+
  xlab("treatment")+
  ylab("predicted likelihood of parasite presence\n per myxozoan taxon per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  scale_linetype_discrete(name = c("life history"),labels=c("direct","complex"))+
  theme(legend.position="none")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color=c("cadetblue","burlywood4")),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  annotate("text",label="effect of treatment:\np = 0.058",x = 1.9, y = 0.045, size = 6)+
  scale_x_discrete(limits=(levels(myxo_predict$x)),labels=c("control","impact"))
  
myxo_plot

