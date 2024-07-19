setwd("C:/Users/Test 1/OneDrive/Documents/WOODLAB/TUBRI")


#SIZE DATA FROM DISSECTED FISHES
#ICTALURUS PUNCTATUS

#processed data
library(readr)
ictpun <- read_csv("data/processed/Ictalurus_punctatus_processed_machine_readable_UPDATED_2024.07.10.csv")
View(ictpun)
ictpunraw<-Ictalurus_punctatus_Datasheet_2024_06_28


#rawdata
View(ictpunraw)



#plotting ALL
plot(ictpunraw$StandardLength_mm~ictpunraw$YearCollected)

#which rows contain fish where selection was influenced by size?
#13, 17, 35, 39, 44, 46, 56, 57, 65, 67, 69, 70, 74, 79, 80, 87
ictpuntrim <- ictpunraw[-c(13, 17, 35, 39, 44, 46, 56, 57, 65, 67, 69, 70, 74, 79, 80, 87), ]
View(ictpuntrim)

#plot only randomly selected fishes
plot(ictpuntrim$StandardLength_mm~ictpuntrim$YearCollected)
summary(lm(ictpuntrim$StandardLength_mm~ictpuntrim$YearCollected))
abline(lm(ictpuntrim$StandardLength_mm~ictpuntrim$YearCollected))



#NOTROPIS ATHERNOIDES
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
summary(lm(notroptrim$StandardLength_mm~notroptrim$YearCollected)) 
abline(lm(notroptrim$StandardLength_mm~notroptrim$YearCollected))
#yay! it looks like we have a strong trend in decreasing body size for notropis fishes




#SIZE DATA FROM FULLY RESELECTED FISHES

install.packages("googlesheets4")
library(googlesheets4)
desmond_data <- read_sheet("https://docs.google.com/spreadsheets/d/1Ix7ZkoTA7AZDOnA3Rqmxt8cB9hPGSwEiwbmCek9uUOQ/edit?gid=0#gid=0") #paste the URL of the google sheet, then follow the prompts
gs4_auth(scopes="spreadsheets.readonly") #says that people with access to this code can only read the spreadsheet! follow prompts and ALLOW tidyverse to see your spreadsheets
View(desmond_data)
#main research question: how did body length change over time?
library(ggplot2)

#plotting both species together
SL_time<-ggplot(desmond_data,aes(Year,SL, color=Species))+
  scale_color_manual(values=c("#0571b0","#ca0020"))+
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

#above, I added a line based on a linear model that represents the data! In this case, it is super simple. We can see there isn't a strong trend in any direction in this data for ictalurus, so let's take a closer look at a linear model for this data. 

#MODELS FOR RESELECTED FISHES

#Ictalurus punctatus
ict_model <- lm(SL~Year, data= desmond_data)
summary(ict_model)
#the model summary shows that there is not a very good fit of this model to the data (low R-sq value), and that there is no significant slope (pvalue, or Pr>t, is much greater than 0.05). 

#we'll make the same type of model for notropis, but since we don't have a ton of data yet, we'll wait!