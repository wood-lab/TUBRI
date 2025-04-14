### Myxo measurements

library(readxl)
library(ggplot2)

myxo_measures <- read_excel("Manuscripts/Myxozoan_measurements.xlsx")
View(myxo_measures)

myxo_measuresb <-subset(myxo_measures,IndividualFishID != "Myxobolus ovalis" )
myxo_measuresc <-subset(myxo_measures,IndividualFishID != "Myxobolus ovalis" &
                          IndividualFishID != "Myxobolus bellus" &
                          IndividualFishID != "Myxobolus discrepans" &
                          IndividualFishID != "Myxobolus obliquus" &
                          IndividualFishID != "Myxobolus meglitschi" )

# Spore lengths
ggplot(myxo_measuresb, aes(x= IndividualFishID,
                         y=Spore_length))+
  geom_point()+apatheme+ylab("Spore length (µm)")+
  ggtitle("Spore length")

ggsave(file="Manuscripts/Figures/Measurements/sporelength.png", width=240, height=180, dpi=1000, units = "mm")

# Spore width
ggplot(myxo_measuresb, aes(x= IndividualFishID,
                          y=Spore_width))+
  geom_point()+apatheme+ylab("Spore width (µm)")+
  ggtitle("Spore width")

ggsave(file="Manuscripts/Figures/Measurements/sporewidth.png", width=240, height=180, dpi=1000, units = "mm")

# Polar capsule ratio
myxo_measuresc$length_ratio <- (myxo_measuresc$L_pc_length)/(myxo_measuresc$R_pc_length)
myxo_measuresc$width_ratio <- (myxo_measuresc$L_pc_width)/(myxo_measuresc$R_pc_width)

ggplot(myxo_measuresc, aes(x= IndividualFishID,
                          y=length_ratio))+
  geom_point()+apatheme+ylab("LPC/RPC length")+
  ggtitle("Left and right polar capsule (LPC/RPC) length ratio ")

ggsave(file="Manuscripts/Figures/Measurements/pc_length_ratio.png", width=240, height=180, dpi=1000, units = "mm")

ggplot(myxo_measuresc, aes(x= IndividualFishID,
                          y=width_ratio))+
  geom_point()+apatheme+ylab("LPC/RPC width")+
  ggtitle("Left and right polar capsule (LPC/RPC) width ratio ")

ggsave(file="Manuscripts/Figures/Measurements/pc_width_ratio.png", width=240, height=180, dpi=1000, units = "mm")
