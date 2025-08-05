#SuPeR PARASITES SUMMER 2024
#"TITLE OF PROJECT"
#Code written by Gabby Commisso and Desmond Boyd


setwd("C:/Users/Test 1/OneDrive/Documents/WOODLAB/TUBRI")
library(googlesheets4)
library(ggplot2)
library(lme4)
library(lmerTest)
library(ggeffects)
library(readr)
#################################################################SIZE DATA FROM DISSECTED FISHES
#############################ictalurus
#loading in dada
ictpunraw<-Ictalurus_punctatus_Datasheet_2024_06_28
View(ictpunraw)

#plotting ALL
plot(ictpunraw$StandardLength_mm~ictpunraw$YearCollected)

#which rows contain fish where selection was influenced by size?
#13, 17, 35, 39, 44, 46, 56, 57, 65, 67, 69, 70, 74, 79, 80, 87
ictpuntrim <- ictpunraw[-c(13, 17, 35, 39, 44, 46, 56, 57, 65, 67, 69, 70, 74, 79, 80, 87), ]
View(ictpuntrim)

#plot only randomly selected fishes
ggplot(data=ictpuntrim, aes(YearCollected, StandardLength_mm))+
  geom_point(size=2)+
  geom_smooth(method=lm, color="blue")+
  xlab("Year collected")+
  ylab("Standard Length in mm")+
  annotate("text", x = -Inf, y = Inf, 
           label = paste("p-value: 0.637", "\n",
                         "slope: 0.067"),
           hjust = -0.5, vjust = 1.5, 
           size = 4, 
           color = "black", 
           fontface = "italic") +
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))

#simple linear model summary
summary(lm(ictpuntrim$StandardLength_mm~ictpuntrim$YearCollected))

################################notropis
#loading in data
notropraw<-Notropis_atherinoides_Datasheet_2024_01_10
View(notropraw)

#plotting ALL
plot(notropraw$StandardLength_mm~notropraw$YearCollected)

#which rows contain size-selected fish?
#35, 36, 37, 54, 75, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 96, 99, 100, 101, 102, 106, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 121, 123, 125, 126, 127, 128, 129, 130, 131, 132, 133, 136, 137, 138, 139, 141, 142, 144, 145, 146, 147, 148, 149, 151, 152, 154, 155, 156, 158, 159, 160, 161, 162, 164, 165, 168, 169, 171, 172, 173, 174, 175, 176, 177, 178, 179, 181, 182, 183, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208
notroptrim <- notropraw[-c(35, 36, 37, 54, 75, 81, 82, 83, 84, 85, 86, 87, 88, 90, 91, 92, 93, 94, 96, 99, 100, 101, 102, 106, 110, 111, 112, 113, 114, 115, 116, 117, 118, 119, 120, 121, 123, 125, 126, 127, 128, 129, 130, 131, 132, 133, 136, 137, 138, 139, 141, 142, 144, 145, 146, 147, 148, 149, 151, 152, 154, 155, 156, 158, 159, 160, 161, 162, 164, 165, 168, 169, 171, 172, 173, 174, 175, 176, 177, 178, 179, 181, 182, 183, 185, 186, 187, 188, 189, 190, 191, 192, 193, 194, 195, 196, 197, 198, 199, 200, 201, 202, 203, 204, 205, 206, 207, 208), ]
View(notroptrim)

#plot only randomly selected fishes
ggplot(data=notroptrim, aes(YearCollected, StandardLength_mm))+
  geom_point(size=2)+
  geom_smooth(method=lm, color="red")+
  xlab("Year collected")+
  ylab("Standard Length in mm")+
  annotate("text", x = -Inf, y = Inf, 
           label = paste("p-value: 0.001", "\n",
                         "slope: -0.274"),
           hjust = -0.5, vjust = 1.5, 
           size = 4, 
           color = "black", 
           fontface = "italic") +
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))

#model summary
summary(lm(notroptrim$StandardLength_mm~notroptrim$YearCollected)) 

#############################hybog
#loading in data
hybograw <-read_csv("data/raw/Hybognathus_nuchalis_Datasheet_2024.07.25.csv")
View(hybograw)

#plotting ALL
ggplot(data=hybograw, aes(YearCollected, StandardLength_mm), color="#006008")+
  geom_point(size=2)+
  geom_smooth(method=lm, color="#006008")+
  xlab("Year collected")+
  ylab("Standard Length in mm")+
  annotate("text", x = -Inf, y = Inf, 
           label = paste("p-value: 0.012", "\n",
                         "slope: -0.172"),
           hjust = -0.5, vjust = 1.5, 
           size = 4, 
           color = "black", 
           fontface = "italic") +
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))

#create model and add regression line
summary(lm(hybograw$StandardLength_mm~hybograw$YearCollected)) 
#nearly every fish was selected for size; even still there is a significant size decrease



#######################################################SIZE DATA FROM FULLY RESELECTED FISHES

desmond_data <- read_sheet("https://docs.google.com/spreadsheets/d/1Ix7ZkoTA7AZDOnA3Rqmxt8cB9hPGSwEiwbmCek9uUOQ/edit?gid=0#gid=0") #paste the URL of the google sheet, then follow the prompts
gs4_auth(scopes="spreadsheets.readonly") #says that people with access to this code can only read the spreadsheet! follow prompts and ALLOW tidyverse to see your spreadsheets

#main research question: how did body length change over time?
View(desmond_data)


###############################all species
SL_time<-ggplot(desmond_data,aes(Year,SL, color=Species))+
  scale_color_manual(values=c("#008000", "blue","red"))+
  geom_point(size=2)+
  xlab("Year collected")+
  ylab("Standard Length in mm")+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
 SL_time

###########################ictalurus
#subset
ict_desmond <- subset(desmond_data, Species == "Ictalurus_punctatus")
View(ict_desmond)
#model
ict_model <- lmerTest::lmer(SL~Year + (1|Cat_Num), data= ict_desmond)
ict_sum <- summary(ict_model)
plot(ict_model)
#predicting line
ict_pred <- ggpredict(ict_model, c("Year [n=100]"))
#plotting
plot(ict_pred, show_data = TRUE, dot_alpha = 1, color = "blue")+ 
  labs(x = 'Year Collected', y = 'Standard Length in mm',title=NULL)+
  annotate("text", x = -Inf, y = Inf, 
           label = paste("p-value: 0.855", "\n",
             "slope: 0.098"),
            hjust = -0.5, vjust = 1.5, 
            size = 4, 
            color = "black", 
            fontface = "italic") +
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))

#############################notropis
#subset
notrop_desmond <- subset(desmond_data, Species == "Notropis_atherinoides")
#model
notrop_model <- lmerTest::lmer(SL~Year + (1|Cat_Num), data= notrop_desmond)
summary(notrop_model)
plot(notrop_model)
#predicting line
notrop_pred <- ggpredict(notrop_model, c("Year [n=100]"))
#plot
plot(notrop_pred, show_data = TRUE, dot_alpha = 1, color = "red")+ 
  labs(x = 'Year Collected', y = 'Standard Length in mm',title=NULL)+
  annotate("text", x = -Inf, y = Inf, 
           label = paste("p-value: 0.009", "\n",
                         "slope: -0.513"),
           hjust = -0.5, vjust = 1.8, 
           size = 4, 
           color = "black", 
           fontface = "italic") +
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))

##############################hybog
#subset
hybog_desmond <- subset(desmond_data, Species == "Hybognathus_nuchalis")
View(hybog_desmond)
#model
hybog_model <- lmerTest::lmer(SL~Year + (1|Cat_Num), data= hybog_desmond)
summary(hybog_model)
plot(hybog_model)
#predicting line
hybog_pred <- ggpredict(hybog_model, c("Year [n=100]"))
#plot
plot(hybog_pred, show_data = TRUE, dot_alpha = 1, color = "#008000")+ 
  labs(x = 'Year Collected', y = 'Standard Length in mm',title=NULL)+
  annotate("text", x = -Inf, y = Inf, 
           label = paste("p-value: 0.0004", "\n",
                         "slope: -0.724"),
           hjust = -0.5, vjust = 1.75, 
           size = 4, 
           color = "black", 
           fontface = "italic") +
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))








##################################################WHAT ABOUT WEIGHT:LENGTH RATIO?
# this data won't be used for august's poster session, but by plotting the graphs you can see we show a pretty intense decrease in weight to body length ratio! this means that not only are fish getting shorter, they're getting thinner for their size. this suggests that the issue may not be exclusively age class shifting/suttkus fishing out the sampling sites, but something else
###################subset data
#ictalurus
ict_desmond$weight_to_sl <- ict_desmond$Weight/ict_desmond$SL
View(ict_desmond)
#notropis
notrop_desmond$weight_to_sl <- notrop_desmond$Weight/notrop_desmond$SL
View(notrop_desmond)
#hybog
hybog_desmond$weight_to_sl <- hybog_desmond$Weight/hybog_desmond$SL
View(hybog_desmond)

###################plots
#ictalurus
#model
ict_model_sl <- lmerTest::lmer(weight_to_sl~Year + (1|Cat_Num), data= ict_desmond)
ict_sum_sl <- summary(ict_model_sl)
plot(ict_model_sl)
#predicting line
ict_pred_sl <- ggpredict(ict_model_sl, c("Year [n=100]"))
#plotting
plot(ict_pred_sl, show_data = TRUE, dot_alpha = 1, color = "blue")+ 
  labs(x = 'Year Collected', y = 'Weight-to-Length Ratio (mg/mm)',title=NULL)+
  annotate("text", x = -Inf, y = Inf, 
           label = paste("p-value: 0.455", "\n",
                         "slope: 1.143"),
           hjust = -0.5, vjust = 1.5, 
           size = 4, 
           color = "black", 
           fontface = "italic") +
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))


#notropis
#model
notrop_model_sl <- lmerTest::lmer(weight_to_sl~Year + (1|Cat_Num), data= notrop_desmond)
summary(notrop_model_sl)
plot(notrop_model)
#predicting line
notrop_pred_sl <- ggpredict(notrop_model_sl, c("Year [n=100]"))
#plot
plot(notrop_pred_sl, show_data = TRUE, dot_alpha = 1, color = "red")+ 
  labs(x = 'Year Collected', y = 'Weight-to-Length Ratio (mg/mm)',title=NULL)+
  annotate("text", x = -Inf, y = Inf, 
           label = paste("p-value: 0.0266", "\n",
                         "slope: -0.422"),
           hjust = -0.5, vjust = 1.8, 
           size = 4, 
           color = "black", 
           fontface = "italic") +
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))


#hybog
#model
hybog_model_sl <- lmerTest::lmer(weight_to_sl~Year + (1|Cat_Num), data= hybog_desmond)
summary(hybog_model_sl)
plot(hybog_model_sl)
#predicting line
hybog_pred_sl <- ggpredict(hybog_model_sl, c("Year [n=100]"))
#plot
plot(hybog_pred_sl, show_data = TRUE, dot_alpha = 1, color = "#008000")+ 
  labs(x = 'Year Collected', y = 'Standard Length in mm',title=NULL)+
  annotate("text", x = -Inf, y = Inf, 
           label = paste("p-value: 0.0028", "\n",
                         "slope: -1.159"),
           hjust = -0.5, vjust = 1.75, 
           size = 4, 
           color = "black", 
           fontface = "italic") +
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))


#Possible weight/SL~Year Summaries
ict_model1 <- lmerTest::lmer((Weight/SL)~Year + (1|Cat_Num), data= ict_desmond)
summary(ict_model1)

notrop_model1 <- lmerTest::lmer((Weight/SL)~Year + (1|Cat_Num), data= notrop_desmond)
summary(notrop_model1)

hybog_model1 <- lmerTest::lmer((Weight/SL)~Year + (1|Cat_Num), data= hybog_desmond)
summary(hybog_model1)

install.packages("usethis")
library(usethis)
use_git_config(user.name= "Fisherguy18", user.email="djboyd@email.sc.edu")






library(glmmTMB)


#ELIANA'S PROJECT 2025
hybog_polished <- read_csv("data/processed/Hybognathus_nuchalis_processed_human_readable_UPDATED_2024.08.01.csv")
hybog_polished <- hybog_polished %>% 
  janitor::clean_names()
View(hybog_polished)

summary(lm(hybograw$TotalLength_mm~hybograw$YearCollected))

hybog_model <- glmmTMB(total_length_mm ~ year_collected + ci + (1|catalog_number), data= hybog_polished)
summary(hybog_model)

hybog_model1 <- glmmTMB(total_length_mm ~ year_collected + ci, data= hybog_polished)
summary(hybog_model1)

AIC(hybog_model)
AIC(hybog_model1)


hybog_pred <- ggpredict(hybog_model, c("year_collected", "ci"))

plot(hybog_pred, show_data = TRUE, dot_alpha = 1) + 
  labs(x = 'Year Collected', y = 'Standard Length in mm',title=NULL) +
  scale_color_manual(values = c("control" = "#76A973", 
                                "impact" = "#825F87")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())


gamaff_polished <- read_csv("data/processed/Gambusia_affinis_processed_human_readable_UPDATED_2024.08.17.csv")
gamaff_polished <- gamaff_polished %>%
  janitor::clean_names()
View(gamaff_polished)

gamaff_model <- glmmTMB(total_length_mm ~ year_collected + ci + (1|catalog_number), data= gamaff_polished)
summary(gamaff_model)

gamaff_model1 <- glmmTMB(total_length_mm ~ year_collected + ci, data= gamaff_polished)
summary(gamaff_model1)

AIC(gamaff_model)
AIC(gamaff_model1)

gamaff_pred <- ggpredict(gamaff_model, c("year_collected", "ci"))

plot(gamaff_pred, show_data = TRUE, dot_alpha = 1) + 
  labs(x = 'Year Collected', y = 'Standard Length in mm',title=NULL) +
  scale_color_manual(values = c("control" = "#76A973", 
                                "impact" = "#825F87")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line())


######using desmond's data
big_model <- glmmTMB(TL~Year + Species + (1|Cat_Num), data=desmond_data)
summary(big_model)

big_pred <- ggpredict(big_model, c("Year", "Species"))

plot(big_pred, show_data = TRUE, show_ci=TRUE, dot_alpha = 1) + 
  labs(x = 'Year Collected', y = 'Standard Length in mm',title=NULL) +
  scale_color_manual(values = c("Ictalurus_punctatus" = "#deebf7", 
                                "Notropis_atherinoides" = "#9ecae1", 
                                "Hybognathus_nuchalis" = "#3182bd")) +
  scale_fill_manual(values = c("Ictalurus_punctatus" = "#deebf7", 
                               "Notropis_atherinoides" = "#9ecae1", 
                               "Hybognathus_nuchalis" = "#3182bd")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        legend.direction = "horizontal", 
        legend.position = "top")

plot(big_pred, show_data = TRUE, show_ci=TRUE, dot_alpha = 1) + 
  labs(x = 'Year Collected', y = 'Total Length in mm',title=NULL) +
  scale_color_manual(values = c("Ictalurus_punctatus" = "#08519c", 
                                "Notropis_atherinoides" = "#9ecae1", 
                                "Hybognathus_nuchalis" = "#3182bd")) +
  scale_fill_manual(values = c("Ictalurus_punctatus" = "#08519c", 
                               "Notropis_atherinoides" = "#9ecae1", 
                               "Hybognathus_nuchalis" = "#3182bd")) +
  theme_bw() +
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border=element_blank(),
        axis.line=element_line(),
        legend.direction = "horizontal", 
        legend.position = "top")
