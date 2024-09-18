## This script curates chemical data from USGS including nutrients and element concentrations
# Daki Diaz-Morales
# dmdiaz@uw.edu; diazdakeishla@gmail.com


### Nutrients----
## Import data for nutrients

nutrients <- read_csv("data/Physicochemical/nutrients_resultphyschem.csv")

# Read USGS-gage site metadata

siteinfo <-  read_csv("data/Physicochemical/station_info.csv")

# Add info on sampling sites (latitude, longitude, CI, etc.)

nutrients_withmetadata <- merge(nutrients, siteinfo, by.x = "MonitoringLocationIdentifier", by.y = "MonitoringLocationIdentifier", all.x = TRUE)

# Summarize by taking the mean per year

nutrients_means <- nutrients_withmetadata %>% group_by(YearCollected,CI,Measure) %>% 
  summarise(Result=mean(Result,na.rm=TRUE),
            .groups = 'drop')

# Nitrogen species: nitrogen available to plants occurs in the form of 
# ammonia (NH3), ammonium (NH4+), nitrate (NO3−), nitrite (NO2-), and 
# organic nitrogen. However, the data that we have most complete is that of inorganic N
# (nitrate and nitrite). Nitrate and nitrite have data for most years except between 1977-1984. 
# However, we have data for that window of time as "inorganic nitrogen (nitrat and nitrite).
# Therefore, we use that to fill the gap. 


# Take only nitrate, nitrite, and their sum (inorganic nitrogen (nitrate and nitrite))

nitrate_nitrite <- subset(nutrients_means, Measure == "Nitrite" |
                            Measure == "Nitrate"|
                            Measure == "Inorganic nitrogen (nitrate and nitrite)")


# From the following plot you can see that nitrate and nitrite have the longest time span
# Inorganic nitrogen was measure in most recent years. So we will the gap in nitrate and nitrite with inorganic N.

# Evaluate dynamic of nutrients across years

ggplot(nutrients_means, aes(x= as.factor(YearCollected),
                            y=Result, color=CI))+
  geom_point()+apatheme+
  ggtitle("Nutrients per year")+
  facet_wrap("Measure")


# We remove data before 1977 and after 1984 to not be reduntant on nitrate and nitrite concentatrions individually measured for those years.
inorgN_a <- nitrate_nitrite %>% filter(!(Measure == "Inorganic nitrogen (nitrate and nitrite)" &
                                               YearCollected < 1977))

inorganicN <- inorgN_a %>% filter(!(Measure == "Inorganic nitrogen (nitrate and nitrite)" &
                                           YearCollected > 1984))

# OK, now we have to still sum individual measurements of nitrate and nitrite to have total inorganic (nitrate and nitrite)

inorganicN$Measure <- as.factor(inorganicN$Measure)

inorganicN <- inorganicN %>% 
  mutate(Measure_type = case_when(Measure == "Nitrate" ~ "nitrate_nitrite",
                                Measure == "Nitrite" ~ "nitrate_nitrite",
                              Measure == "Inorganic nitrogen (nitrate and nitrite)"~ 
                              "inorg_N"))


inorganicN_sum <- inorganicN %>% group_by(YearCollected,CI,Measure_type) %>% 
  summarise(Result=sum(Result,na.rm=TRUE),
            .groups = 'drop')


# Re-evaluate dynamic of nitrates and nitrites across years

ggplot(inorganicN_sum, aes(x= as.factor(YearCollected),
                           y=Result, color=CI))+
  geom_point()+apatheme+
  ggtitle("Total nitrate and nitrite concentration per year")

### Phosphates did not change much over time so we keep it out to reduce dimensionality