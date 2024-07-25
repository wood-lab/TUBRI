#SuPeR PARASITES SUMMER 2024
#"TITLE OF PROJECT"
#Code written by Gabby Commisso and Desmond Boyd


setwd("C:/Users/Test 1/OneDrive/Documents/WOODLAB/TUBRI")
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
plot(ictpuntrim$StandardLength_mm~ictpuntrim$YearCollected)

#create model and add regression line
summary(lm(ictpuntrim$StandardLength_mm~ictpuntrim$YearCollected))
abline(lm(ictpuntrim$StandardLength_mm~ictpuntrim$YearCollected))

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
plot(notroptrim$StandardLength_mm~notroptrim$YearCollected)

#create model and add regression line
summary(lm(notroptrim$StandardLength_mm~notroptrim$YearCollected)) 
abline(lm(notroptrim$StandardLength_mm~notroptrim$YearCollected))

#############################hybog
#loading in data
hybograw <-Hybognathus_nuchalis_Datasheet_2024_07_25
View(hybograw)

#plotting ALL
plot(hybograw$StandardLength_mm~hybograw$YearCollected)

#create model and add regression line
summary(lm(hybograw$StandardLength_mm~hybograw$YearCollected)) 
abline(lm(hybograw$StandardLength_mm~hybograw$YearCollected))
#nearly every fish was selected for size; even still there is a significant size decrease



#######################################################SIZE DATA FROM FULLY RESELECTED FISHES

install.packages("googlesheets4")
library(googlesheets4)
desmond_data <- read_sheet("https://docs.google.com/spreadsheets/d/1Ix7ZkoTA7AZDOnA3Rqmxt8cB9hPGSwEiwbmCek9uUOQ/edit?gid=0#gid=0") #paste the URL of the google sheet, then follow the prompts
gs4_auth(scopes="spreadsheets.readonly") #says that people with access to this code can only read the spreadsheet! follow prompts and ALLOW tidyverse to see your spreadsheets

#main research question: how did body length change over time?
library(ggplot2)

View(desmond_data)
#plotting both species together
SL_time<-ggplot(desmond_data,aes(Year,SL, color=Species))+
  scale_color_manual(values=c("#FCD12A", "#0571b0","#ca0020"))+
  geom_point(size=4)+
  xlab("Year collected")+
  ylab("Standard Length in mm")+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
 SL_time

#plotting only ictalurus 
ict_only <- ggplot(subset(desmond_data, Species %in% "Ictalurus_punctatus"), aes(Year,SL))+
  geom_point(size=4)+
  xlab("Year collected")+
  ylab("Standard Length in mm")+
  geom_smooth(method = "lm", formula = y ~ x)+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
ict_only

#plotting only notropis
notrop_only <- ggplot(subset(desmond_data, Species %in% "Notropis_atherinoides"), aes(Year,SL))+
  geom_point(size=4)+
  xlab("Year collected")+
  ylab("Standard Length in mm")+
  geom_smooth(method = "lm", formula = y ~ x)+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
notrop_only

#plotting only hybog
hybog_only <- ggplot(subset(desmond_data, Species %in% "Hybognathus_nuchalis"), aes(Year,SL))+
  geom_point(size=4)+
  xlab("Year collected")+
  ylab("Standard Length in mm")+
  geom_smooth(method = "lm", formula = y ~ x)+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
hybog_only


##################################################WHAT ABOUT WEIGHT:LENGTH RATIO?
##########subset data
#ictalurus
ict_desmond <- subset(desmond_data, Species == "Ictalurus_punctatus")
ict_desmond$weight_to_sl <- ict_desmond$Weight/ict_desmond$SL
View(ict_desmond)
#notropis
notrop_desmond <- subset(desmond_data, Species == "Notropis_atherinoides")
notrop_desmond$weight_to_sl <- notrop_desmond$Weight/notrop_desmond$SL
View(notrop_desmond)
#hybog
hybog_desmond <- subset(desmond_data, Species == "Hybognathus_nuchalis")
hybog_desmond$weight_to_sl <- hybog_desmond$Weight/hybog_desmond$SL
View(hybog_desmond)

#########plots
#ictalurus
ict_w2SL_only <- ggplot(ict_desmond, aes(Year,weight_to_sl))+
  geom_point(size=4)+
  xlab("Year collected")+
  ylab("Weight/Length Ratio")+
  geom_smooth(method = "lm", formula = y ~ x)+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
ict_w2SL_only

#notropis
notrop_w2SL_only <- ggplot(notrop_desmond, aes(Year,weight_to_sl))+
  geom_point(size=4)+
  xlab("Year collected")+
  ylab("Weight/Length Ratio")+
  geom_smooth(method = "lm", formula = y ~ x)+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
notrop_w2SL_only

#hybog
hybog_w2SL_only <- ggplot(hybog_desmond, aes(Year,weight_to_sl))+
  geom_point(size=4)+
  xlab("Year collected")+
  ylab("Weight/Length Ratio")+
  geom_smooth(method = "lm", formula = y ~ x)+
  theme_minimal()+
  theme(plot.title=element_text(size=14,hjust=0.5,face="plain"),axis.text.y=element_text(size=14),axis.title.y=element_text(size=14),axis.text.x=element_text(size=14),axis.title.x=element_text(size=10),panel.background=element_rect(fill="white",color="black"),panel.grid.major=element_line(color="grey95"),panel.grid.minor=element_line(color=NA),plot.margin=unit(c(0,0,0,0),"cm"))
hybog_w2SL_only

##################################################MODELS FOR RESELECTED FISHES
#ictalurus
ict_model <- lm(SL~Year, data= ict_desmond)
summary(ict_model)
#the model summary shows that there is not a very good fit of this model to the data (low R-sq value), and that there is no significant slope (pvalue, or Pr>t, is much greater than 0.05). 

#notropis
notrop_model <- lm(SL~Year, data= notrop_desmond)
summary(notrop_model)

#hybog
hybog_model <- lm(SL~Year, data= hybog_desmond)
summary(hybog_model)