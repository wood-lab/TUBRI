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




### PRELIMINARY ANALYSIS - CHELSEA, 22 JULY 2024

full_data<-read.csv("data/processed/Full_dataset_with_psite_life_history_info_2024.07.19.csv", header = T, sep = ",")
colnames(full_data)[20]<-"scaled_TL_mm"

View(full_data)

# First create a model where you are looking at Life_History (direct versus complex)
# Note that this model may take a few minutes to run.

library(lme4)
model_1<-glmer.nb(psite_count~CI*Life_History+(1|Fish_sp.x/psite_spp.x)+
                                offset(log(scaled_TL_mm)),data=full_data,family="nbinom")

# At some point, you should probably control for season
# (1|MonthCollected)

summary(model_1)

# Let's save the model output down somewhere so that we can grab it again when we need it without having to wait
# for the model to run and converge

library(stargazer)
model_1_output<-stargazer::stargazer(model_1, type = "text")
write.table(model_1_output,"REU/shyanne/model_1_output_2024.07.15")

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
  theme(legend.position="top",legend.title = element_text(size = 18),
        legend.text = element_text(size=14))
color_plot




# Now create a model where you are looking at Parasite_taxonomic_group (Nematoda, Trematoda, etc)
# Note that this model will take a WHILE to run. It took about 37 minutes on my computer.

library(lme4)
model_2<-glmer.nb(psite_count~CI*Parasite_taxonomic_group+(1|Fish_sp.x/psite_spp.x)+
                    offset(log(scaled_TL_mm)),data=full_data,family="nbinom")


# At some point, you should probably control for season
# (1|MonthCollected)

summary(model_2)

# Check VIFs to make sure that there is no collinearity
library(car)
vif(model_2)


# Let's save the model output down somewhere so that we can grab it again when we need it without having to wait
# for the model to run and converge

model_2_output<-stargazer::stargazer(model_2, type = "text")
write.table(model_2_output,"REU/shyanne/model_2_output_2024.07.15")

read.table("REU/shyanne/model_2_output_2024.07.15")

# Now let's plot it

library(ggeffects)
color_predict_2<-ggeffect(model_2,c("Parasite_taxonomic_group","CI"))

# Need to exclude Hirudinea because error bars are wacky (probably because we found so few). For now, I've
# just taken the error bars off.

# Set up a nice color palette first.
library(viridis)
pal<-viridis(n=9)

color_plot_2<-ggplot(color_predict_2,aes(group,predicted),group=Parasite_taxonomic_group,
                   color=Parasite_taxonomic_group)+
  geom_point(aes(group=x,color=x),size=4,pch=19)+
  #geom_errorbar(data=color_predict_2,mapping=aes(x="CI",ymin=conf.low,ymax=conf.high),width=0.03)+
  geom_line(aes(group=x,color=x,linetype=x))+
  scale_color_manual(name = c("parasite taxonomic group"),values=pal)+
  xlab("treatment")+
  ylab("predicted parasite abundance\n per parasite taxon per host individual")+
  theme_minimal()+
  #labs(linetype="parasite life history strategy")+
  theme(plot.title=element_text(size=18,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=16),
        axis.text.x=element_text(size=18,color=c("cadetblue","burlywood4")),axis.title.x=element_text(size=16),
        panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color=NA),
        panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))+
  scale_x_discrete(limits=(levels(color_predict_2$group)))+
  theme(legend.position="left",legend.title = element_text(size = 18),
        legend.text = element_text(size=14))
color_plot_2

