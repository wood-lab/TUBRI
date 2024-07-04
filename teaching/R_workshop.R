# How to use R
# A workshop
# Written by Chelsea Wood (chelwood@uw.edu)
# 4 July 2024
# Hi SuPeR Parasites!!! HELLO
# Hello my name is IMANI and I am working with DAKI. <3

# To run a thing in R, you put your cursor on it and hit COMMAND+RETURN (on Mac) or CONTROL+RETURN (on PC).
# Anything with a # in front of it will not run, because it is annotation - words designed to be human-readable,
# but meant to be ignored by R.

# This is how you find out what directory R is reading from - i.e., what folder on your computer is R going to
# when you ask it to retrieve a file?  getwd = get working directory.

getwd()

# If you need to change your working directory, paste the file structure between the quotation marks below.
# To get the file structure, open the folder and then "Get info" about that folder. Copy the file structure
# and paste it below.

setwd("C:/Users/imani/OneDrive/Desktop/TUBRI_Monogenea_Project")

# Once your working directory is set, you're ready to read in the data!  If there are sub-folders inside your 
# working directory, you'll need to specify them as I have below.
# The <- command tells R what a thing is called.  So you can read the line below as,
# Look at this csv file, and name it pim_vig_data.  When I call pim_vig_data, I want you to give me the csv file.

#pim_vig_data<-read.csv("data/processed/Ictalurus_punctatus_processed_human_readable.csv")
pim_vig_data<-read.csv("Ictalurus_punctatus_processed_human_readable.csv")

# To see your data as a spreadsheet in a new tab, use the View command.

View(pim_vig_data)

# You can also see your data in the console below by running the name of the dataset.

pim_vig_data

# With your dataset ready to go, you can now call individual columns in that dataset using $ plus the variable
# name, exactly as it is spelled and capitalized in the dataset.

pim_vig_data$YearCollected

pim_vig_data$MONO.IP

# You can also call individual rows, columns, or cells using brackets.

pim_vig_data[33,2]

# You can now perform simple data manipulations on the variables.  For example, if I want the total number of 
# MONO.IP worms we counted in Pimephales vigilax, I can write,

sum(pim_vig_data$MONO.IP)

# Or maybe I want the mean number of MONO.IP per fish. It's simple as this:

mean(pim_vig_data$MONO.IP)

# I can even call up a quick histogram of the number of MONO.IP in each fish:

hist(pim_vig_data$MONO.IP)

# Or I can plot the number of MONO.IP over time - plot(y~x).

plot(pim_vig_data$MONO.IP~pim_vig_data$YearCollected)

plot(pim_vig_data$MONO.IP~jitter(pim_vig_data$YearCollected,10))

plot(pim_vig_data$TREM.ALLO~pim_vig_data$YearCollected)

plot(pim_vig_data$TREM.ALLO~jitter(pim_vig_data$YearCollected,10))

# Looks like there's a pattern there!  Let's do a simple statistical test to see if it is significant?
# A model is an equation that represents the data. Here, our model is y = mx + b.

summary(lm(pim_vig_data$MONO.IP~pim_vig_data$YearCollected))

summary(lm(pim_vig_data$TREM.ALLO~pim_vig_data$YearCollected))
