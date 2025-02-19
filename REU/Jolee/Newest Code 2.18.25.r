### Jolee Code REU Project ###
### created 7/11/24        ###
### by Jolee Thirtyacre    ###
### and Connor Whalen
## completely updated and moved from JT_Markdown_EXT.Rmd on 8/16/24


# Welcome
#For helpful instructions on how to start using R Markdown, visit [this link](https://fish497.github.io/website/lectures/week_07/intro_rmarkdown.html).
#To insert code chunks, press command+option+i (mac), ctrl+alt+i (windows)
#Anything inside a code chunk will run as normal code, anything outside a code chunk will appear in the document as plain text. 
#data used for this project are located on the [Wood Lab GitHub](https://github.com/wood-lab/TUBRI)
#this project is part of a larger project focused on reconstructing change in parasite burden in river fishes over time. If you want to learn more about the larger project, you can do so at the [Wood Lab Website](https://chelsealwood.com/research/historical-parasite-ecology/_)

# Questions we are asking: 
#Is increased parasitism associated with decreased body condition?
#This is complicated because we are proposing an observational study and the causality that we are attempting to describe can go in both directions. Fish that have high parasitism could be skinner, but also skinner fish could be more susceptible to increased parasitism. We won’t be able to parse this out using our dataset, but it is still an interesting question worth answering. We just need to be careful to not make any absolute statements based on this (which is why we use the language ‘associated’ instead of ‘affected by’). 
#Is the type (taxonomic distinction and/or life stage) of parasite important in determining whether parasitism is associated with decreased body condition?
#Parasites exist in different parts of the body based on what they are: monogenes will always be on the gills, metacercariae can be most anywhere on the body, etc. 
#To answer this question, we would need to go through the data and identify what parasites are in their adult stage versus their larval stage. This would involve a bit of research into what the different (larval or adult) stages are present for each of the ones we have found. 

# Methods
#see code below for full explainations of how parasites were grouped and analyzed
### Rules for discriminating between adult and larval: 
#by looking at our compiled parasite ID guides, we will designate all found parasites as either adult or larval on a case-by-case basis
#For any parasite we have identified to a specific taxonomic level (family, genus, etc.), we will look at that life cycle in the literature and ID what stage we have. 


#Analyses

#### Install packages
install.packages("tidyverse")
install.packages("janitor")
install.packages("ggeffects")
install.packages("rmarkdown")
install.packages("caret")
install.packages("emmeans")
install.packages("tibble")
install.packages("lmerTest")
install.packages("DHARMa")
install.packages("ggpmisc")
install.packages("smatr")
install.packages("scales")
install.packages("lme4")
install.packages("Matrix")

install.packages("car")
install.packages("ggpubr")
install.packages("glmmTMB")

#### Load packages
library(Matrix)
library(ggeffects)
library(tidyverse)
library(lmerTest)
library(emmeans)
library(broom)
library(ggpubr)
library(smatr)
library(stringr)
library(caret)
library(tibble)
library(car)
library(tidyr)
library(readr)
library(dbplyr)
library(lubridate)
library(ggplot2)
library(glmmTMB)
library(DHARMa)
library(scales)
library(ggpmisc)
library(lme4)


setwd("E:/SENIOR Senior Year/Autumn Quarter/Parasite Research/Correct Data")

### Load datasets 
P.vigilax.data <- read.csv("Pimephales_vigilax_processed_human_readable_UPDATED_2024.08.16.csv") %>% 
  janitor::clean_names() %>%
  mutate(weight_mg = as.numeric(weight_mg)) %>% 
  mutate(weight_g = (weight_mg/1000))

I.punctatus.data <- read.csv("Ictalurus_punctatus_processed_human_readable.csv_UPDATED_2024.08.25.csv")  %>%
  janitor::clean_names() %>%
  mutate(day_collected = as.integer(day_collected)) %>%
  mutate(weight_mg = as.numeric(weight_mg))%>% 
  mutate(weight_g = (weight_mg/1000))

N.atherinoides.data <- read.csv("Notropis_atherinoides_processed_human_readable_UPDATED_2024.08.16.csv") %>% 
  janitor::clean_names() %>%
  mutate(weight_mg = as.numeric(weight_mg)) %>% 
  mutate(weight_g = (weight_mg/1000))

H.nuchalis.data <- read.csv("Hybognathus_nuchalis_processed_human_readable_UPDATED_2024.08.01.csv") %>%
  janitor::clean_names() %>%
  mutate(weight_mg = as.numeric(weight_mg)) %>% 
  mutate(weight_g = (weight_mg/1000))

P.vigil.data <- read.csv("Percina_vigil_processed_human_readable_UPDATED_2024.08.01.csv") %>%
  janitor::clean_names() %>%
  mutate(day_collected = as.integer(day_collected)) %>%
  mutate(weight_mg = as.numeric(weight_mg)) %>% 
  mutate(weight_g = (weight_mg/1000))

C.velifer.data <- read.csv("Carpiodes_velifer_processed_human_readable_UPDATED_2024.08.16.csv") %>%
  janitor::clean_names() %>%
  mutate(weight_mg = as.numeric(weight_mg)) %>% 
  mutate(weight_g = (weight_mg/1000))

G.affinis.data <- read.csv("Gambusia_affinis_processed_human_readable_UPDATED_2024.08.17.csv") %>%
  janitor::clean_names() %>%
  mutate(weight_mg = as.numeric(weight_mg)) %>% 
  mutate(weight_g = (weight_mg/1000))

# Question 1
## Giving adult and larval IDs so we can look at life stage compared to total parasite burden (for answering question 1).

P.vigilax.data <- P.vigilax.data %>%
  mutate(adult = rowSums(select(., 
                                "acanth_bigb", 
                                "acanth_bkr",
                                "nem_cona",
                                "nem_cap",
                                "nem_dich",
                                "nem_larv",
                                "nem_trbk",
                                "nem_unk",
                                "cest_godz",
                                "cest_ntg",
                                "cope_cl",
                                "cope_grab",
                                "cest_pea",
                                "cest_l",
                                "cest_vit",
                                "cest_unk",
                                "mono_gyro",
                                "mono_dact",
                                "mono_lg",
                                "trem_rnda",
                                "trem_unk"))) %>%
  mutate(larval = rowSums(select(.,
                                 "meta_hs",
                                 "meta_unk", 
                                 "myx_go",
                                 "cest_comp",
                                 "myxo_bc",
                                 "myx_thel",
                                 "myxo_sbad",
                                 "myxo_unk",
                                 "myxo_up",
                                 "trem_met",
                                 "trem_sss",
                                 "trem_buc",
                                 "trem_metaf",
                                 "trem_metags",
                                 "trem_ms",
                                 "trem_nea"))) %>%
  mutate(adult = as.integer(adult)) %>% 
  mutate(larval = as.integer(larval)) %>%
  mutate(psite_count = rowSums(select(.,
                                      "adult", 
                                      "larval"))) %>% 
  mutate(weight_g = (weight_mg/1000)) %>%
  janitor::clean_names()

####
I.punctatus.data <- I.punctatus.data %>%
  mutate(larval = rowSums(select(., 
                                 "nem_cyst",
                                 "trem_meta_gg",
                                 "trem_meta_unk",
                                 "trem_lg",
                                 "nem_cyst",
                                 "cest_comp",
                                 "myx_tail"))) %>%
  mutate(adult = rowSums(select(.,
                                "acan_ostr",
                                "cest_meg",
                                # "leech_kut", removed entirely so it's not counted
                                "nem_nmts",
                                "mono_ip",
                                "mono_unk",
                                "nem_phar",
                                "nem_pty",
                                "nem_sp",
                                "nem_taf",
                                "nem_unk",
                                # "nmorph", removed entirely so it's not counted
                                "trem_2",
                                "trem_allo",
                                "trem_ble",
                                "trem_crep",
                                "trem_f",
                                "trem_lip",
                                "trem_megict",
                                "trem_smile",
                                "trem_phyls",
                                "trem_unk"))) %>%
  # select(-nmorph, -leech_kut) %>% # this removes the column but still counts the parasite in the adult/larval/psite_count
  mutate(adult = as.integer(adult)) %>% 
  mutate(larval = as.integer(larval)) %>%
  mutate(psite_count = rowSums(select(.,
                                      "adult", 
                                      "larval"))) %>% 
  mutate(weight_g = (weight_mg/1000)) %>%
  janitor::clean_names()

####
N.atherinoides.data <- N.atherinoides.data %>%
  mutate(larval = rowSums(select(., 
                                 "nem_cyst",
                                 "trem_larv",
                                 "trem_meta",
                                 "myx_sp", 
                                 "cest_comp",
                                 "cest_lvs",
                                 "myx_round",
                                 "trem_cyst"))) %>%#
  mutate(adult = rowSums(select(.,
                                "cest_hyd",
                                "cest_tri",
                                "cope_a",
                                "mono_all",
                                "mono_unk",
                                "nem_bran",
                                "nem_cap",
                                "nem_spky",
                                "nem_tfs",
                                "nem_twas",
                                "nem_unkn",
                                "trem_gb",
                                "trem_l"))) %>%
  mutate(adult = as.integer(adult)) %>% 
  mutate(larval = as.integer(larval)) %>%
  mutate(psite_count = rowSums(select(.,
                                      "adult", 
                                      "larval"))) %>% 
  mutate(weight_g = (weight_mg/1000)) %>%
  janitor::clean_names()

####
H.nuchalis.data <- H.nuchalis.data %>%
  mutate(larval = rowSums(select(., 
                                 "nem_daft",
                                 "trem_cm",
                                 "trem_meta_het",
                                 "trem_meta_sp",
                                 "trem_bc",
                                 "trem_meta_go",
                                 "trem_meta_es",
                                 "meta_unk",
                                 "trem_meta",
                                 "trem_meta_unk",
                                 "myx_ak",
                                 "myx_fi",
                                 "myx_tin",
                                 "myx_gl",
                                 "myx_ot",
                                 "myx_eye",
                                 "myx_gc"))) %>%#
  mutate(adult = rowSums(select(.,
                                "cest_bb",
                                "nem_unk",
                                "trem_ass",
                                "trem_ns",
                                "trem_diplo",
                                "trem_unk",
                                "cest_unk",
                                "mono_dact",
                                "mono_gdac",
                                "acan_ad",
                                "cope_pood"))) %>%
  mutate(adult = as.integer(adult)) %>% 
  mutate(larval = as.integer(larval)) %>%
  mutate(psite_count = rowSums(select(.,
                                      "adult", 
                                      "larval"))) %>% 
  mutate(weight_g = (weight_mg/1000)) %>%
  janitor::clean_names()

#### 
P.vigil.data <- P.vigil.data %>%
  mutate(adult = rowSums(select(., 
                                "nem_heath",
                                "nem_squig",
                                "acanth_bb",
                                "trem_phlb",
                                "cest_bot",
                                "mono_gyro"))) %>%#
  mutate(larval = rowSums(select(.,
                                 "nem_larv",
                                 "acanth_cya",
                                 "acanth_thin",
                                 "trem_sk",
                                 "trem_chonk",
                                 "trem_endous",
                                 "trem_cl",
                                 "trem_ispot"))) %>%
  mutate(adult = as.integer(adult)) %>% 
  mutate(larval = as.integer(larval)) %>%
  mutate(psite_count = rowSums(select(.,
                                      "adult", 
                                      "larval"))) %>% 
  mutate(weight_g = (weight_mg/1000)) %>%
  janitor::clean_names()

#### 
C.velifer.data <- C.velifer.data %>%
  mutate(adult = rowSums(select(., 
                                "cest_ws",
                                "cest_megan",
                                "cest_tri",
                                "cest_arch",
                                "cope_hny",
                                "nem_cerl",
                                "trem_sep",
                                "trem_ors",
                                "trem_cv"))) %>%#
  mutate(larval = rowSums(select(.,
                                 "trem_larv",
                                 "trem_meta_unk", # added, is this correct or are we excluding it in the study?
                                 "myx_s",
                                 "myx_f",
                                 "myx_e",
                                 "myx_g",
                                 "myx_cm",
                                 "myx_kl"))) %>%
  mutate(adult = as.integer(adult)) %>% 
  mutate(larval = as.integer(larval)) %>%
  mutate(psite_count = rowSums(select(.,
                                      "adult", 
                                      "larval"))) %>% 
  mutate(weight_g = (weight_mg/1000)) %>%
  janitor::clean_names()

####
G.affinis.data <- G.affinis.data %>%
  mutate(larval = rowSums(select(., 
                                 "acanth_nemor",
                                 "acanth_blve",
                                 "cest_mipi",
                                 "cest_unk",
                                 "nem_larv", 
                                 "nem_unk",
                                 "meta_bolism",
                                 "meta_z",
                                 "mono_sase",
                                 "mono_g",
                                 "myx_stwy", 
                                 "myx_thin",
                                 "myx_upc"))) %>%
  mutate(adult = rowSums(select(.,
                                "cest_woad",
                                "cest_y",
                                "cest_triga",
                                "nem_oodle",
                                "nem_w",
                                "nem_pretz",
                                "nem_esis",
                                "trem_shy",
                                "trem_scorp",
                                "trem_ga",
                                "trem_pos"))) %>%
  mutate(adult = as.integer(adult)) %>% 
  mutate(larval = as.integer(larval)) %>%
  mutate(psite_count = rowSums(select(.,
                                      "adult", 
                                      "larval"))) %>% 
  mutate(weight_g = (weight_mg/1000)) %>%
  janitor::clean_names()

## Create new columns for calculating condition factor (using smi - scaled mass index)

P.vigilax.data.smi <- P.vigilax.data %>%
  mutate(log_weight = log(weight_g)) %>%
  mutate(log_standard_length = log(standard_length_mm)) %>%
  janitor::clean_names()

I.punctatus.data.smi <- I.punctatus.data %>%
  mutate(log_weight = log(weight_g)) %>%
  mutate(log_standard_length = log(standard_length_mm)) %>% 
  mutate(b_smi = (log_weight - log_standard_length)) %>%
  janitor::clean_names()

N.atherinoides.data.smi <- N.atherinoides.data %>%
  mutate(log_weight = log(weight_g)) %>%
  mutate(log_standard_length = log(standard_length_mm)) %>% 
  mutate(b_smi = (log_weight - log_standard_length)) %>%
  janitor::clean_names()

H.nuchalis.data.smi <- H.nuchalis.data %>%
  mutate(log_weight = log(weight_g)) %>%
  mutate(log_standard_length = log(standard_length_mm)) %>% 
  mutate(b_smi = (log_weight - log_standard_length)) %>%
  janitor::clean_names()

P.vigil.data.smi <- P.vigil.data %>%
  mutate(log_weight = log(weight_g)) %>%
  mutate(log_standard_length = log(standard_length_mm)) %>% 
  mutate(b_smi = (log_weight - log_standard_length)) %>%
  janitor::clean_names()

C.velifer.data.smi <- C.velifer.data %>%
  mutate(log_weight = log(weight_g)) %>%
  mutate(log_standard_length = log(standard_length_mm)) %>% 
  mutate(b_smi = (log_weight - log_standard_length)) %>%
  janitor::clean_names()

G.affinis.data.smi <- G.affinis.data %>%
  mutate(log_weight = log(weight_g)) %>%
  mutate(log_standard_length = log(standard_length_mm)) %>% 
  mutate(b_smi = (log_weight - log_standard_length)) %>%
  janitor::clean_names()


## get scaling exponent
pimvig_m <- lm(log_weight ~ log_standard_length, data = P.vigilax.data.smi)
pimvig_c <- tidy(pimvig_m)
pimvig_s <- pimvig_c %>%
  filter(term == "log_standard_length") %>%
  pull(estimate) # 3.199591

ictpun_m <- lm(log_weight ~ log_standard_length, data = I.punctatus.data.smi)
ictpun_c <- tidy(ictpun_m)
ictpun_s <- ictpun_c %>%
  filter(term == "log_standard_length") %>%
  pull(estimate) # 2.89255

notath_m <- lm(log_weight ~ log_standard_length, data = N.atherinoides.data.smi)
notath_c <- tidy(notath_m)
notath_s <- notath_c %>%
  filter(term == "log_standard_length") %>%
  pull(estimate) # 3.13197476

hybnuc_m <- lm(log_weight ~ log_standard_length, data = H.nuchalis.data.smi)
hybnuc_c <- tidy(hybnuc_m)
hybnuc_s <- hybnuc_c %>%
  filter(term == "log_standard_length") %>%
  pull(estimate) # 2.765758

pervig_m <- lm(log_weight ~ log_standard_length, data = P.vigil.data.smi)
pervig_c <- tidy(pervig_m)
pervig_s <- pervig_c %>%
  filter(term == "log_standard_length") %>%
  pull(estimate) # 3.1723261

carvel_m <- lm(log_weight ~ log_standard_length, data = C.velifer.data.smi)
carvel_c <- tidy(carvel_m)
carvel_s <- carvel_c %>%
  filter(term == "log_standard_length") %>%
  pull(estimate) # 3.26232447

gamaff_m <- lm(log_weight ~ log_standard_length, data = G.affinis.data.smi)
gamaff_c <- tidy(gamaff_m)
gamaff_s <- gamaff_c %>%
  filter(term == "log_standard_length") %>%
  pull(estimate) # 3.6403417


## Find average length for all pooled samples (per species)

pimvig_avg_l <- mean(P.vigilax.data$standard_length_mm, na.rm = TRUE)
ictpun_avg_l <- mean(I.punctatus.data$standard_length_mm, na.rm = TRUE)
notath_avg_l <- mean(N.atherinoides.data$standard_length_mm, na.rm = TRUE)
hybnuc_avg_l <- mean(H.nuchalis.data$standard_length_mm, na.rm = TRUE)
pervig_avg_l <- mean(P.vigil.data$standard_length_mm, na.rm = TRUE)
carvel_avg_l <- mean(C.velifer.data$standard_length_mm, na.rm = TRUE)
gamaff_avg_l <- mean(G.affinis.data$standard_length_mm, na.rm = TRUE)


## Calculate new SMI Condition Factor  (uses average standard length for it's numerator)
P.vigilax.data.smi <- P.vigilax.data.smi %>%
  mutate(smi_pimvig_2 = weight_g*(pimvig_avg_l / standard_length_mm) ^ pimvig_s) %>%
  janitor::clean_names()

I.punctatus.data.smi <- I.punctatus.data.smi %>%
  mutate(smi_ictpun_2 = weight_g*(ictpun_avg_l / standard_length_mm) ^ ictpun_s) %>%
  janitor::clean_names()

N.atherinoides.data.smi <- N.atherinoides.data.smi %>%
  mutate(smi_notath_2 = weight_g*(notath_avg_l / standard_length_mm) ^ notath_s) %>%
  janitor::clean_names()

H.nuchalis.data.smi <- H.nuchalis.data.smi %>%
  mutate(smi_hybnuc_2 = weight_g*(hybnuc_avg_l / standard_length_mm) ^ hybnuc_s) %>%
  janitor::clean_names()

P.vigil.data.smi <- P.vigil.data.smi %>%
  mutate(smi_pervig_2 = weight_g*(pervig_avg_l / standard_length_mm) ^ pervig_s) %>%
  janitor::clean_names()

C.velifer.data.smi <- C.velifer.data.smi %>%
  mutate(smi_carvel_2 = weight_g*(carvel_avg_l / standard_length_mm) ^ carvel_s) %>%
  janitor::clean_names()

G.affinis.data.smi <- G.affinis.data.smi %>%
  mutate(smi_gamaff_2 = weight_g*(gamaff_avg_l / standard_length_mm) ^ gamaff_s) %>%
  janitor::clean_names()



### pull out the specific fish with 0 parasites and fish with only 1 parasite group > 0

## functions filtering fish w/o parasites & w/ one species

# function filtering fish w/ 1 parasite species
infected.smi <- function(data, species_columns) {
  data %>%
    filter(rowSums(select(., all_of(species_columns)) > 0) == 1) %>%
    mutate(state = ("infected")) %>%
    select(everything())
}

# function filtering fish w/o parasites
uninfected.smi <- function(data, species_columns) {
  data %>%
    filter(rowSums(select(., all_of(species_columns)) > 0) == 0) %>%
    mutate(state = ("uninfected")) %>%
    select(everything())
}

# function filtering fish w/o parasites & w/ one species
fish.filter.smi <- function(data, species_columns) {
  data %>%
    filter(rowSums(select(., all_of(species_columns)) > 0) <= 1) %>%
    select(everything())
}


# P. vigilax
P_vigilax_psite_columns <- c("acanth_bigb","acanth_bkr", "nem_cona", "nem_cap",
                             "nem_dich", "nem_larv", "nem_trbk", "nem_unk",
                             "cest_godz", "cest_ntg", "cope_cl", "cope_grab",
                             "cest_pea", "cest_l", "cest_vit", "cest_unk",
                             "mono_gyro", "mono_dact", "mono_lg", "meta_hs",
                             "meta_unk", "myx_go", "cest_comp", 
                             "myxo_bc", "myx_thel", "myxo_sbad", "myxo_unk",
                             "myxo_up", "trem_met", "trem_buc", "trem_metaf",
                             "trem_metags", "trem_ms","trem_rnda", "trem_unk",
                             "trem_sss", "trem_nea")

P.vigilax.infected.smi <- infected.smi(P.vigilax.data.smi, P_vigilax_psite_columns) %>%
  mutate(state = ("infected")) %>%
  mutate(species = ("P.vigilax")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_pimvig_2) %>%
  janitor::clean_names()

P.vigilax.uninfected.smi <- uninfected.smi(P.vigilax.data.smi, P_vigilax_psite_columns) %>%
  mutate(state = ("uninfected")) %>%
  mutate(species = ("P.vigilax")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_pimvig_2) %>%
  janitor::clean_names()

P.vigilax.filtered.smi <- rbind(P.vigilax.infected.smi, P.vigilax.uninfected.smi) # fish w/ 1 and w/o parasites


# I. punctatus 

I_punctatus_psite_columns <- c("nem_cyst", "trem_meta_gg", "trem_meta_unk", "trem_lg",
                               "cest_comp", "myx_tail", "trem_unk",
                               "acan_ostr", "cest_meg", "nem_nmts",
                               "mono_ip", "mono_unk", "nem_phar", "nem_pty",
                               "nem_sp", "nem_taf", "nem_unk",
                               "trem_2", "trem_allo", "trem_ble", "trem_crep",
                               "trem_f", "trem_lip", "trem_megict", "trem_smile", "trem_phyls")


I.punctatus.infected.smi <- infected.smi(I.punctatus.data.smi, I_punctatus_psite_columns) %>%
  mutate(state = ("infected")) %>% # added for graphing purposes, see below
  filter(!(individual_fish_id %in% c("124399_02"))) %>% # removed one fish infected w/ leech
  mutate(species = ("I.punctatus")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_ictpun_2) %>%
  janitor::clean_names()


I.punctatus.uninfected.smi <- uninfected.smi(I.punctatus.data.smi, I_punctatus_psite_columns) %>%
  filter(!(individual_fish_id %in% c("182143_01"))) %>% # removed one fish infected w/ nematomorph
  mutate(state = ("uninfected")) %>% # added for graphing purposes, see below
  mutate(species = ("I.punctatus")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_ictpun_2) %>%
  janitor::clean_names()


I.punctatus.filtered.smi <- rbind(I.punctatus.infected.smi, I.punctatus.uninfected.smi) %>%  # fish w/ 1 and w/o parasites
  filter(!(individual_fish_id %in% c("182143_01"))) %>% # removed one fish infected w/ nematomorph
  filter(!(individual_fish_id %in% c("124399_02"))) %>%
  select(-nmorph, -leech_kut)

# N. atherinoides

N_atherinoides_psite_columns <- c("nem_cyst", "trem_larv", "trem_meta", "myx_sp", 
                                  "cest_comp", "cest_lvs", "myx_round", "trem_cyst",
                                  "cest_hyd", "cest_tri", "cope_a", "mono_all",
                                  "mono_unk", "nem_bran", "nem_cap", "nem_spky",
                                  "nem_tfs", "nem_twas", "nem_unkn","trem_gb",
                                  "trem_l")


N.atherinoides.infected.smi <- infected.smi(N.atherinoides.data.smi, N_atherinoides_psite_columns) %>%
  mutate(state = ("infected")) %>%
  mutate(species = ("N.atherinoides")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_notath_2) %>%
  janitor::clean_names()


N.atherinoides.uninfected.smi <- uninfected.smi(N.atherinoides.data.smi, N_atherinoides_psite_columns) %>%
  mutate(state = ("uninfected")) %>%
  mutate(species = ("N.atherinoides")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_notath_2) %>%
  janitor::clean_names()


N.atherinoides.filtered.smi <- rbind(N.atherinoides.infected.smi, N.atherinoides.uninfected.smi) # fish w/ 1 and w/o parasites


# H. nuchalis

H_nuchalis_psite_columns <- c("nem_daft", "trem_cm", "trem_meta_het", "trem_meta_sp",
                              "trem_bc", "trem_meta_go", "trem_meta_es", "meta_unk",
                              "trem_meta", "trem_meta_unk", "myx_ak", "myx_fi",
                              "myx_tin", "myx_gl", "myx_ot", "myx_eye", "myx_gc",
                              "cest_bb", "nem_unk", "trem_ass", "trem_ns",
                              "trem_diplo", "trem_unk", "cest_unk", "mono_dact",
                              "mono_gdac", "acan_ad", "cope_pood")

H.nuchalis.infected.smi <- infected.smi(H.nuchalis.data.smi, H_nuchalis_psite_columns) %>%
  mutate(state = ("infected")) %>%
  mutate(species = ("H.nuchalis")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_hybnuc_2) %>%
  janitor::clean_names()

H.nuchalis.uninfected.smi <- uninfected.smi(H.nuchalis.data.smi, H_nuchalis_psite_columns) %>%
  mutate(state = ("uninfected")) %>%
  mutate(species = ("H.nuchalis")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_hybnuc_2) %>%
  janitor::clean_names()

H.nuchalis.filtered.smi <- rbind(H.nuchalis.infected.smi, H.nuchalis.uninfected.smi) # fish w/ 1 and w/o parasites


# P. vigil

P_vigil_psite_columns <- c("nem_heath", "nem_squig", "acanth_bb", "trem_phlb",
                           "cest_bot", "mono_gyro", "nem_larv", "acanth_cya",
                           "acanth_thin", "trem_sk", "trem_chonk", "trem_endous",
                           "trem_cl", "trem_ispot")

P.vigil.infected.smi <- infected.smi(P.vigil.data.smi, P_vigil_psite_columns) %>%
  filter(!(individual_fish_id %in% c("33226_01", "72382_03")))  %>%  # removed two fish infected with nem_unk
  mutate(state = ("infected")) %>%
  mutate(species = ("P.vigil"))  %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_pervig_2) %>%
  janitor::clean_names()

P.vigil.uninfected.smi <- uninfected.smi(P.vigil.data.smi, P_vigil_psite_columns)  %>%
  filter(!(individual_fish_id %in% c("33226_01", "72382_03")))  %>%  # removed two fish infected with nem_unk
  mutate(state = ("uninfected")) %>%
  mutate(species = ("P.vigil")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_pervig_2) %>%
  janitor::clean_names()

P.vigil.filtered.smi <- rbind(P.vigil.infected.smi, P.vigil.uninfected.smi) # fish w/ 1 and w/o parasites


# C. velifer

C_velifer_psite_columns <- c("cest_ws", "cest_megan", "cest_tri", "cest_arch",
                             "cope_hny", "nem_cerl", "trem_sep", "trem_ors",
                             "trem_cv", "trem_larv", "myx_s", "myx_f",
                             "myx_e", "myx_g", "myx_cm",
                             "myx_kl", "cest_unk", "cili_fu",
                             "trem_meta_unk", "myx_unk")

C.velifer.infected.smi <- infected.smi(C.velifer.data.smi, C_velifer_psite_columns)  %>%
  mutate(state = ("infected")) %>%
  mutate(species = ("C.velifer")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_carvel_2) %>%
  janitor::clean_names()

C.velifer.uninfected.smi <- uninfected.smi(C.velifer.data.smi, C_velifer_psite_columns)%>%
  mutate(state = ("uninfected")) %>%
  mutate(species = ("C.velifer")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_carvel_2) %>%
  janitor::clean_names()

C.velifer.filtered.smi <- rbind(C.velifer.infected.smi, C.velifer.uninfected.smi) # fish w/ 1 and w/o parasites


# G. affinis

G_affinis_psite_columns <- c("acanth_nemor", "acanth_blve", "cest_mipi", "cest_unk", 
                             "nem_larv", "nem_unk", "meta_bolism","meta_z", 
                             "mono_sase", "mono_g", "myx_stwy", "myx_thin", 
                             "myx_upc", "cest_woad", "cest_y", "cest_triga", 
                             "nem_oodle", "nem_w", "nem_pretz", "nem_esis", 
                             "trem_shy", "trem_scorp", "trem_ga", "trem_pos",
                             "trem_unk")

G.affinis.infected.smi <- infected.smi(G.affinis.data.smi, G_affinis_psite_columns) %>%
  mutate(state = ("infected")) %>%
  mutate(species = ("G.affinis")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_gamaff_2) %>%
  janitor::clean_names()

G.affinis.uninfected.smi <- uninfected.smi(G.affinis.data.smi, G_affinis_psite_columns) %>%
  mutate(state = ("uninfected")) %>%
  mutate(species = ("G.affinis")) %>%
  mutate(stan_count = (standard_length_mm / psite_count)) %>% # to standardize psite count
  rename(smi = smi_gamaff_2) %>%
  janitor::clean_names()

G.affinis.filtered.smi <- rbind(G.affinis.infected.smi, G.affinis.uninfected.smi) # fish w/ 1 and w/o parasites


# binding the datasets together
full.data.filtered.smi <- bind_rows(P.vigilax.filtered.smi,
                                    I.punctatus.filtered.smi,
                                    C.velifer.filtered.smi,
                                    G.affinis.filtered.smi,
                                    H.nuchalis.filtered.smi,
                                    N.atherinoides.filtered.smi,
                                    P.vigil.filtered.smi)

# adjusting some column names for clarity before pivoting to long format
full.data.filtered.smi <- full.data.filtered.smi %>%
  rename(fish_species = species) %>%
  rename(total_psite_count = psite_count) %>%
  janitor::clean_names()


#### LONG FORMAT TRANSITION
# switching from wide format to long
full.data.filtered.smi.long <- pivot_longer(data = full.data.filtered.smi,
                                            cols = matches("^(acanth_|cest_|cope_|meta_|
                                            |mono_|myx_|myxo_|nem_|trem_|cili_|acan_)"),
                                            names_to = "psite_species",
                                            values_to = "psite_species_count") %>%
  janitor::clean_names()

# removing rows where infected fish psite_species_count = 0 or NA
remove_infected_0s <- function(data) {
    cleaned.data <- data %>%
      filter(!(state == "infected" & 
                 psite_species_count == 0)) %>% # leaves only the psite the infected fish is infected with
      filter(!(is.na(psite_species_count))) # remove NAs b/c NAs are for psites a fish species doesn't have
  return(cleaned.data)
}

full.data.filtered.smi.long.clean <- remove_infected_0s(full.data.filtered.smi.long)

## check to make sure it's removing the right rows and leaving uninfected rows alone (tested w/ 3 fish spp.)
# P. vigilax -
# should be 8 infected (so 8 infected entries)
# should be 4 uninfected (w/ 37 psite spp. so should be 4*37=148 uninfected entries)
print(full.data.filtered.smi.long.clean %>% 
        filter(state == "uninfected") %>% 
        group_by(fish_species == "P.vigilax") %>%
        count())
# 148 uninfected entries -- works!
print(full.data.filtered.smi.long.clean %>%
        filter(state == "infected") %>%
        group_by(fish_species == "P.vigilax") %>%
        count())
# 8 infected entries -- works!

# I. punctatus
# 32 infected (should be 32 entries)
# 21 uninfected (w/ 26 psite spp. should be 21*26=546 uninfected entries)
print(full.data.filtered.smi.long.clean %>% 
        filter(state == "uninfected") %>% 
        group_by(fish_species == "I.punctatus") %>%
        count())
# 546 uninfected entries -- works!
print(full.data.filtered.smi.long.clean %>%
        filter(state == "infected") %>%
        group_by(fish_species == "I.punctatus") %>%
        count())
# 32 infected entries -- works!

# G. affinis
# 35 infected
# 27 uninfected (w/ 25 psite spp. should be  27*25=675 uninfected entries)
print(full.data.filtered.smi.long.clean %>% 
        filter(state == "uninfected") %>% 
        group_by(fish_species == "G.affinis") %>%
        count())
# 675 uninfected entries -- works!
print(full.data.filtered.smi.long.clean %>%
        filter(state == "infected") %>%
        group_by(fish_species == "G.affinis") %>%
        count())
# 35 infected entries -- works!

# subsetting full long for uninfected long
full.uninfected.data.filtered.smi.long <- full.data.filtered.smi.long.clean %>%
  filter(state == "uninfected") # 3938 obs.
#check that this is right
sum(full.data.filtered.smi.long.clean$state == "uninfected") # 3938 obs. -- works!

# subsetting full long for infected long
full.infected.data.filtered.smi.long.clean <- full.data.filtered.smi.long.clean %>%
  filter(state == "infected") # 286 obs.
#check that this is right
sum(full.data.filtered.smi.long.clean$state == "infected") # 286 obs. -- works!





## Model formula for Question 1 
# K_i\sim parasite\ count\ + ci
#impact category is anything below the pulp mill
#control category is anything above the pulp mill

## Test model for question 1
#run models

## Models for new SMI (using length numerator) comparison with uninfected and infected data
summary(smi_pimvig_lm_2 <- lm(smi_pimvig_2 ~ psite_count + ci+offset(standard_length_mm), data = P.vigilax.filtered.smi)) 
summary(smi_ictpun_lm_2 <- lm(smi_ictpun_2 ~ psite_count + ci+offset(standard_length_mm), data = I.punctatus.filtered.smi)) ### @ a = 0.05 (0.0248)
## total parasite count is associated with an increase in smi when accounting for the impact category in ictpun
summary(smi_notath_lm_2 <- lm(smi_notath_2 ~ psite_count + ci+offset(standard_length_mm), data = N.atherinoides.filtered.smi)) 
summary(smi_hybnuc_lm_2 <- lm(smi_hybnuc_2 ~ psite_count + ci+offset(standard_length_mm), data = H.nuchalis.filtered.smi)) ### @ a=0.05 (ci 0.0282)
## impact category is associated with an increase in smi when accounting for total parasite count in hybnuc
summary(smi_pervig_lm_2 <- lm(smi_pervig_2 ~ psite_count + ci+offset(standard_length_mm), data = P.vigil.filtered.smi)) ### @ a = 0.05 (1.07e-06)
## impact category is associated with an increase in smi when accounting for total parasite count in pervig
summary(smi_carvel_lm_2 <- lm(smi_carvel_2 ~ psite_count + ci+offset(standard_length_mm), data = C.velifer.filtered.smi))
summary(smi_gamaff_lm_2 <- lm(smi_gamaff_2 ~ psite_count + ci+offset(standard_length_mm), data = G.affinis.filtered.smi))


conf_intervals <- confint(smi_ictpun_lm_2)
estimates <- summary(smi_ictpun_lm_2)$coefficients[,"Estimate"]
lower_ci <- conf_intervals[,1]
lower_ci <- conf_intervals[,2]


# model from Daki 10.22.24
summary(smi_pimvig_lm_2 <- lmer(smi_pimvig_2 ~ psite_count + ci+(1|IndividualfishID)+(psite_count|psite_spp.x)+(ci|psite_spp.x), data = P.vigilax.filtered.smi)) 

sjPlot::plot_model(smi_ictpun_lm_2, type="est")

# see how the models are fitting
#### significant results
s1b <- simulateResiduals(fittedModel=smi_ictpun_lm_2,n=250)
plot(s1b)

s1c <- simulateResiduals(fittedModel=smi_hybnuc_lm_2,n=250)
plot(s1c)

s1d <- simulateResiduals(fittedModel=smi_pervig_lm_2,n=250)
plot(s1d)

## Plots question 1

# Generate predictions
mydf1 <- ggpredict(smi_ictpun_lm_2, terms = c("ci","psite_count"))

# Define the theme
apatheme <- theme_bw(base_size = 11, base_family = "sans") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

# Extracted values
p_value_ciimpact1 <- 0.0726  # pulled directly from the model 

# Create the plot with annotations
ggplot(mydf1, aes(x = x, y = predicted)) + 
  geom_point(alpha = 0.6, size = 2) +
  labs(x = 'Category', y = 'Predicted K', title = NULL) +
  annotate("text", x = Inf, y = Inf, label = paste("p-value: ", format(p_value_ciimpact1, scientific = TRUE)),
           hjust = 2.5, vjust = 1.5, size = 4, color = "black") +
  apatheme

#### significant result for CI P.vigil and psite I. punctatus
# Generate predictions
mydf2 <- ggpredict(smi_pervig_lm_2, terms = c("ci"))
mydf3 <- ggpredict(smi_ictpun_lm_2, terms = c("psite_count"))

# Define the theme
apatheme <- theme_bw(base_size = 11, base_family = "sans") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

# Extracted values
p_value_ciimpact2 <- 0.00611  # pulled directly from the model 
p_value_psite_ictpun <- 0.00533 # pulled directly from the model

# Create the plot with annotations
ggplot(mydf2, aes(x = x, y = predicted)) + 
  geom_point(alpha = 0.6, size = 2) +
  labs(x = 'Category', y = 'Predicted SMI', title = NULL) +
  annotate("text", x = Inf, y = Inf, label = paste("p-value: ", format(p_value_ciimpact2, scientific = TRUE)),
           hjust = 2.5, vjust = 1.5, size = 4, color = "black") +
  apatheme

ggplot(mydf3, aes(x = x, y = predicted)) + 
  geom_point(alpha = 0.6, size = 2) +
  labs(x = 'Parasite Count', y = 'Predicted SMI', title = NULL) +
  annotate("text", x = Inf, y = Inf, label = paste("p-value: ", format(p_value_psite_ictpun, scientific = TRUE)),
           hjust = 2.5, vjust = 1.5, size = 4, color = "black") +
  apatheme

### facet wrapping graphs of SMI

# plot all fish species together
ggplot(data = full.data.filtered.smi,
       aes(x = psite_count, y = smi)) +
  geom_point() +
  facet_wrap(~species, ncol = 2) + # make an object facet_labels = species name if you need
  theme(strip.background = element_rect(fill = NA, color = NA), #label background
        panel.background = element_rect(fill = NA,
                                        color = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  scale_x_continuous(expand = c(-0.02, 50)) + # (0.01, 50) when having three columns
  scale_y_continuous(expand = c(0, 0.2)) + # (0, 0.2) when having three columns
  ylab("SMI") +
  xlab("Parasite Count")


# plot just one fish species
plot1 <- ggplot(data = P.vigilax.filtered.smi,
       aes(x = psite_count, y = smi)) +
  geom_point() +
  theme(strip.background = element_rect(fill = NA, color = NA), #label background
        panel.background = element_rect(fill = NA,
                                        color = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  scale_x_continuous(expand = c(0, 2)) +
  scale_y_continuous(expand = c(0.1, 0)) +
  labs(title = "P. vigilax",
       x = "Parasite Count",
       y = "SMI")

plot2 <- ggplot(data = I.punctatus.filtered.smi,
       aes(x = psite_count, y = smi)) +
  geom_point() +
  theme(strip.background = element_rect(fill = NA, color = NA), #label background
        panel.background = element_rect(fill = NA,
                                        color = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  scale_x_continuous(expand = c(0, 2)) +
  scale_y_continuous(expand = c(0.1, 0)) +
  annotate("text", x = Inf, y = Inf, label = paste("p-value: ", format(p_value_psite_ictpun, scientific = TRUE)),
           hjust = 1.5, vjust = 1.5, size = 4, color = "black") +
  labs(title = "I. punctatus",
       x = "Parasite Count",
       y = "SMI")

plot3 <- ggplot(data = C.velifer.filtered.smi,
       aes(x = psite_count, y = smi)) +
  geom_point() +
  theme(strip.background = element_rect(fill = NA, color = NA), #label background
        panel.background = element_rect(fill = NA,
                                        color = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  scale_x_continuous(expand = c(0, 20)) +
  scale_y_continuous(expand = c(0.1, 0)) +
  labs(title = "C. velifer",
       x = "Parasite Count",
       y = "SMI")

plot4 <- ggplot(data = G.affinis.filtered.smi,
       aes(x = psite_count, y = smi)) +
  geom_point() +
  theme(strip.background = element_rect(fill = NA, color = NA), #label background
        panel.background = element_rect(fill = NA,
                                        color = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  scale_x_continuous(expand = c(0, 1)) +
  scale_y_continuous(expand = c(0.1, 0)) +
  labs(title = "G. affinis",
       x = "Parasite Count",
       y = "SMI")

plot5 <- ggplot(data = H.nuchalis.filtered.smi,
       aes(x = psite_count, y = smi)) +
  geom_point() +
  theme(strip.background = element_rect(fill = NA, color = NA), #label background
        panel.background = element_rect(fill = NA,
                                        color = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  scale_x_continuous(expand = c(0, 2)) +
  scale_y_continuous(expand = c(0.1, 0)) +
  labs(title = "H. nuchalis",
       x = "Parasite Count",
       y = "SMI")

plot6 <- ggplot(data = N.atherinoides.filtered.smi,
       aes(x = psite_count, y = smi)) +
  geom_point() +
  theme(strip.background = element_rect(fill = NA, color = NA), #label background
        panel.background = element_rect(fill = NA,
                                        color = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  scale_x_continuous(expand = c(0, 3)) +
  scale_y_continuous(expand = c(0.1, 0)) +
  labs(title = "N. atherinoides",
       x = "Parasite Count",
       y = "SMI")

plot7 <- ggplot(data = P.vigil.filtered.smi,
       aes(x = psite_count, y = smi)) +
  geom_point() +
  theme(strip.background = element_rect(fill = NA, color = NA), #label background
        panel.background = element_rect(fill = NA,
                                        color = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  scale_x_continuous(expand = c(0, 0.3)) +
  scale_y_continuous(expand = c(0.1, 0)) +
  labs(title = "P. vigil",
       x = "Parasite Count",
       y = "SMI")

install.packages("cowplot")
library(cowplot)

plot_grid(plot1, plot2, plot3, plot4,
          plot5, plot6, plot7, ncol = 2)

# plot for I. punctatus significance
ggplot(data = I.punctatus.filtered.smi,
       aes(x = state, y = smi)) +
  geom_boxplot() +
  theme(strip.background = element_rect(fill = NA, color = NA), #label background
        panel.background = element_rect(fill = NA,
                                        color = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  scale_y_continuous(expand = c(0.2, 0.2)) +
  labs(y = "SMI",
       x = "Infection Status",
       title = "I. punctatus") +
  geom_jitter(width=0.2)

ggplot(data = I.punctatus.filtered.smi,
       aes(x = psite_count, y = smi, colour = state)) +
  geom_point() +
  theme(strip.background = element_rect(fill = NA, color = NA), #label background
        panel.background = element_rect(fill = NA,
                                        color = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  scale_x_continuous(expand = c(0, 2)) +
  scale_y_continuous(expand = c(0.1, 0)) +
  labs( y = "SMI",
        x = "Parasite Count",
        colour = "Category") +
  geom_smooth()

# model of I. punctatus
ggplot(mydf3, aes(x = x, y = predicted)) + 
  geom_point(alpha = 1, size = 2.2) +
  labs(x = 'Parasite Count', y = 'Predicted SMI', title = "I. punctatus Model",
           hjust = 2.5, vjust = 1.5, size = 4, color = "black") +
  apatheme +
  geom_smooth(linetype = "dashed", color = "navy", linewidth = 0.1) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.1))
  

# graph for P. vigil CI significance
ggplot(data = P.vigil.filtered.smi,
       aes(x = ci, y = smi)) +
  geom_boxplot() +
  theme(strip.background = element_rect(fill = NA, color = NA), #label background
        panel.background = element_rect(fill = NA,
                                        color = "darkgrey"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing.x = unit(0, "lines"),
        panel.spacing.y = unit(0, "lines")) +
  scale_y_continuous(expand = c(0.2, 0.2)) +
  labs(y = "SMI",
       x = "Category",
       title = "P. vigil - Effluent") # +
  # geom_jitter()

# model for P. vigil
ggplot(mydf2, aes(x = x, y = predicted)) + 
  geom_point(alpha = 1, size = 2) +
  labs(x = 'Category', y = 'Predicted SMI', title = "P. vigil Model") +
  apatheme +
  scale_x_discrete(limits = c("control", "impact")) +
  theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.1))
  

## Discussion
#why are fish in the impact category associated with a higher FCF compared to those in the control category when accounting for parasite count?
# higher is relative to control category, which in turn is higher than the baseline target value of 1. 
#fish in the impact category have a higher FCF than fish in the control category when accounting for total parasite burden, suggesting there is an association between control-impact category and FCF. Fish in the impact category have FCF values that are futher away from normal than control category fish, meaning that they are more 'chubby' then their control counterparts, meaning they are more negatively impacted by being in the impact category relative to control fish, when viewing their body condition relative to our measurement metric, FCF.

# Question 2

#### Import dataset 
# make sure you name all these new datasets and frames something different from what you named them above (that way we can still run all the analyses in total without losing any of the previous ones)

#JT Load Dataset

full.data <- read.csv("Full_dataset_with_psite_life_history_info_2024.09.25.csv") %>% 
  janitor::clean_names()


# CW Load Dataset

full.data <- read.csv("~/Desktop/UW/TUBRI/data/processed/Full_dataset_with_psite_life_history_info_2024.08.01.csv") %>% 
  janitor::clean_names()

### tidy up dataset variables

full.data <- full.data %>%
  mutate(weight_mg = as.numeric(weight_mg)) %>% 
  mutate(weight_g = (weight_mg/1000)) %>% 
  janitor::clean_names()

## Categorize parasites as either adult or larval. All found parasites are included in these analyses (only for the fish species listed). Add back a column for Fulton's Condition Factor. 
#### for the first question

full.data.stages <- full.data %>%
  mutate(
    psite_count = replace_na(psite_count, 0)  
  ) %>%
  mutate(
    life_stage = case_when(
      psite_count == 0 ~ "uninfected",
      psite_spp_x %in% c("ACANTH.BIGB", "ACANTH.BKR", "NEM.CONA", "NEM.CAP",
                         "NEM.DICH", "NEM.LARV", "NEM.TRBK", "NEM.UNK",
                         "CEST.GODZ", "CEST.NTG", "COPE.CL",
                         "COPE.GRAB", "CEST.L", "CEST.VIT",
                         "CEST.UNK", "MONO.GYRO", "MONO.DACT", "MONO.LG",
                         "TREM.RNDA", "TREM.UNK", "ACAN.OSTR",
                         "CEST.MEG", "LEECH.KUT", "NEM.NMTS", "MONO.IP",
                         "MONO.UNK", "NEM.PHAR", "NEM.PTY", "NEM.SP",
                         "NEM.TAF", "NEM.UNK", "NMORPH", "TREM.2",
                         "TREM.ALLO", "TREM.BLE", "TREM.CREP", "TREM.F",
                         "TREM.LIP", "TREM.MEGICT", "TREM.SMILE", "TREM.UNK",
                         "CEST.HYD", "CEST.TRI",
                         "COPE.A", "MONO.ALL", "MONO.UNK", "NEM.BRAN",
                         "NEM.CAP", "NEM.SPKY", "NEM.TFS", "NEM.TWAS",
                         "NEM.UNKN", "TREM.GB", "TREM.L", "CEST.BB", 
                         "NEM.UNK", "TREM.ASS", "TREM.NS",
                         "TREM.DIPLO", "TREM.UNK", "CEST.UNK", "MONO.DACT",
                         "MONO.GDAC", "ACAN.AD", "COPE.POOD", "CEST.WS", "CEST.MEAGAN", "CEST.TRI",
                         "CEST.ARCH", "CEST.UNKA", "COPE.HNY", "NEM.CERL", "NEM.UNK", "TREM.SEP", "TREM.ORS", "TREM.CV", 
                         "MONO.NID", "NEM.HEATH", "MEN.SQUIG", "ACANTH.BB", "TREM.PHLB",
                         "CEST.BOT", "MONO.GYRO","CEST.WOAD", "CEST.Y", "CEST.TRIGA", "NEM.OODLE", "NEM.W",
                         "NEM.PRETZ", "NEM.ESIS", "TREM.SHY", "TREM.SCORP",
                         "TREM.GA", "TREM.POS") ~ "adult",
      
      psite_spp_x %in% c("META.HS", "META.UNK", "MYX.GO", "MYX.GEOM",
                         "MYXO.BC", "MYX.THEL", "CEST.COMP", "CEST.COMP", "CEST.COMP", "CEST.PEA", "CEST.LVS",
                         "MYXO.SBAD", "MYXO.UNK",
                         "MYXO.UP", "TREM.MET", "TREM.SSS", "TREM.BUC",
                         "TREM.METAF", "TREM.METAGS", "TREM.MS", "TREM.NEA",
                         "NEM.CYST", "TREM.META.GG", "TREM.META.UNK", "TREM.LG",
                         "MYX.TAIL", "MYX.GEOM", "NEM.CYST", "NEM.DAFT",
                         "TREM.CM", "TREM.META.HET", "TREM.META.SP", "TREM.BC",
                         "TREM.META.GO", "TREM.META.ES", "MYX.AK", "MYX.FI",
                         "MYX.TIN", "MYX.GL", "MYX.OT", "MYX.EYE",
                         "MYX.GC", "NEM.CYST", "TREM.LARV", "TREM.META",
                         "MYX.SP", "MYX.ROUND", "TREM.CYST", "TREM.LARV", "META.UNK", 
                         "MYX.S", "MYX.F", "MYX.E", "MYX.G", "MYX.CM", "MYX.BNA", "MYX.KL", 
                         "MYX.GM", "NEM.LARV", "ACANTH.CYA", "ACANTH.THIN", 
                         "TREM.SK", "TREM.CHONK", "TREM.ENDOUS", "TREM.CL", "TREM.ISPOT","ACANTH.NEMOR",
                         "ACANTH.BLVE", "CEST.MIPI", "CEST.UNK", "NEM.LARV", 
                         "NEM.UNK", "META.BOLISM", "META.Z", "MONO.SASE", "MONO.G", "MYX.STWY", 
                         "MYX.THIN", "MYX.UPC") ~ "larval",
      
      is.na(psite_spp_x) ~ "0",
    )
  ) %>%
  replace_na(list(life_stage = "uninfected")) %>% 
  mutate(weight_g = (weight_mg/1000)) %>%
  janitor::clean_names()

### fix variable type
full.data.stages <- full.data.stages %>%
  mutate(life_stage = as.factor(life_stage))
full.data.stages <- full.data.stages %>%
  mutate(life_stage = relevel(life_stage,ref='uninfected'))


## New formula for calculating the condition factor for the total dataset
full.data.smi <- full.data.stages %>%
  mutate(log_weight = log(weight_g)) %>%
  mutate(log_standard_length = log(standard_length_mm)) %>%
  janitor::clean_names()

full_data_m <- lm(log_weight ~ log_standard_length, data = full.data.smi)
full_data_c <- tidy(full_data_m)
full_data_s <- full_data_c %>%
  filter(term == "log_standard_length") %>%
  pull(estimate) # 2.85147436

# should be average length not weight
# full.data.smi <- full.data.smi %>%
#  group_by(fish_sp_x) %>%
#  mutate(avg_weight = mean(weight_g, na.rm = TRUE)) %>%
#  ungroup()

full.data.smi <- full.data.smi %>%
  group_by(fish_sp_x) %>%
  mutate(avg_length = mean(standard_length_mm, na.rm = TRUE)) %>%
  ungroup()

full.data.smi <- full.data.smi %>%
  mutate(smi = weight_g*(avg_length / standard_length_mm) ^ full_data_s) %>%
  janitor::clean_names()


## Model formula for Question 2 
#  K_i\sim life\ stage + parasite\ count\ + ci\ + (1|fish\ species)\ + (1|catalog\ number)
# random effects account for things that may be happening to influence the associations we are testing. By including the random effects, we are testing for and proving that autocorrelation (things that will produce similarities/associations in the data in the background) is not playing a role in skewing the model

## Test the model for Question 2
# check how the different variables are interacting with each other 
crosstab <- full.data.smi %>%
  count(life_stage, psite_count)
crosstab_table <- crosstab %>%
  pivot_wider(names_from = psite_count, values_from = n, values_fill = list(n = 0))
view(crosstab_table)

# run the model !!!!!!!!! HAVING ISSUES HERE, ERROR MESSAGE !!!!!!!!!!!!!
library(lmerTest)

summary(q2m <- lmer(smi ~ life_stage + psite_count + ci + (1|fish_sp_x)+ (1|catalog_number), data = full.data.smi)) ## @ a = 0.05 for ci

lmer(smi ~ life_stage + psite_count + ci + (1|fish_sp_x)+ (1|catalog_number), data = full.data.smi)

# check co-linearity using vif function --> GVIF^(1/(2*Df)) under 5 is acceptable generally
vif(q2m)

# check how the model fits
s2=simulateResiduals(fittedModel=q2m,n=250)
plot(s2)

### with new smi
# run the model
summary(q2m_smi <- lmer(smi ~ life_stage + psite_count + ci + (1|fish_sp_x)+ (1|catalog_number), data = full.data.smi))

# check co-linearity using vif function --> GVIF^(1/(2*Df)) under 5 is acceptable generally
vif(q2m_smi)

# check how the model fits
s2_smi=simulateResiduals(fittedModel=q2m,n=250)
plot(s2_smi)


## Generate model predictions and plot
# Generate predictions
mydf <- ggpredict(q2m, terms = c("ci", "life_stage"))

# Define the theme
apatheme <- theme_bw(base_size = 11, base_family = "sans") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

# Extracted values
p_value_ciimpact <- 0.008166  # Example value; replace with actual extracted value

# Create the plot with annotations
ggplot(mydf, aes(x = x, y = predicted, color = group)) + 
  geom_point(alpha = 0.6, size = 2) +
  labs(x = 'Category', y = 'Predicted K', title = NULL) +
  annotate("text", x = Inf, y = Inf, label = paste("p-value control vs. impact catergory: ", format(p_value_ciimpact, scientific = TRUE)),
           hjust = 1.5, vjust = 1.5, size = 4, color = "black") +
  apatheme

#### for new smi
# Generate predictions
mydf_smi <- ggpredict(q2m_smi, terms = c("ci", "life_stage"))

# Define the theme
apatheme_smi <- theme_bw(base_size = 11, base_family = "sans") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line())

# Extracted values
p_value_ciimpact_smi <- 0.008166  # Example value; replace with actual extracted value

# Create the plot with annotations
ggplot(mydf_smi, aes(x = x, y = predicted, color = group)) + 
  geom_point(alpha = 0.6, size = 2) +
  labs(x = 'Category', y = 'Predicted K', title = NULL) +
  annotate("text", x = Inf, y = Inf, label = paste("p-value control vs. impact catergory: ", format(p_value_ciimpact_smi, scientific = TRUE)),
           hjust = 1.5, vjust = 1.5, size = 4, color = "black") +
  apatheme_smi


## Discussion
#why are fish in the impact category assocaited with a higher fcf compared to those in the control category when accounting for parasite count, and life stage of parasites?
#same inference as discussed in question 2 above, but with the addition of the life_stage and random effects variables

# Question 3
#Remove from the dataframe the parasites that can only occur as either an adult or a larvae in the host
#contains option for what parasite taxon is present

## Modify datasets to refect only parasites that are able to occur in larval and adult stages in the host. 
# we are doing this because in order to have a truly comparable dataset of parasites and their association with host body condition, we need to analyze parasites in a uniform way that does not skew the data one direction or another. For example, if you look at trematode parasites, both metacercaria and adult trematodes can live within a fish host. Now, compare that to monogenes who only occur as adults in the fish host. Because their larval stage will not occur in these fish hosts, by including them in our analyses we will be skewing the adult subset of parasites high. The monogenes have no larval counterpart in the fish, so they will always push the adult count higher. By removing monogenes from our analyses, we will have an even ground to compare larval and adult parasites. Therefore, parasites to __NOT__ include in our analyses are: 
#monogenes
#copepods
#acanthocephalans
#myxozoans

full.data.taxon <- full.data %>%
  mutate(
    psite_count = replace_na(psite_count, 0)  
  ) %>%
  mutate(
    life_stage = case_when(
      psite_count == 0 ~ "uninfected", 
      psite_spp_x %in% c("ACANTH.BIGB", "ACANTH.BKR", "NEM.CONA", "NEM.CAP",
                         "NEM.DICH", "NEM.LARV", "NEM.TRBK", "NEM.UNK",
                         "CEST.GODZ", "CEST.NTG", "COPE.CL",
                         "COPE.GRAB", "CEST.L", "CEST.VIT",
                         "CEST.UNK", "MONO.GYRO", "MONO.DACT", "MONO.LG",
                         "TREM.RNDA", "TREM.UNK", "ACAN.OSTR",
                         "CEST.MEG", "LEECH.KUT", "NEM.NMTS", "MONO.IP",
                         "MONO.UNK", "NEM.PHAR", "NEM.PTY", "NEM.SP",
                         "NEM.TAF", "NEM.UNK", "NMORPH", "TREM.2",
                         "TREM.ALLO", "TREM.BLE", "TREM.CREP", "TREM.F",
                         "TREM.LIP", "TREM.MEGICT", "TREM.SMILE", "TREM.UNK",
                         "CEST.HYD", "CEST.TRI",
                         "COPE.A", "MONO.ALL", "MONO.UNK", "NEM.BRAN",
                         "NEM.CAP", "NEM.SPKY", "NEM.TFS", "NEM.TWAS",
                         "NEM.UNKN", "TREM.GB", "TREM.L", "CEST.BB",
                         "NEM.BUD", "NEM.UNK", "TREM.ASS", "TREM.NS",
                         "TREM.DIPLO", "TREM.UNK", "CEST.UNK", "MONO.DACT",
                         "MONO.GDAC", "ACAN.AD", "COPE.POOD") ~ "adult",
      
      psite_spp_x %in% c("META.HS", "META.UNK", "MYX.GO", "MYX.GEOM",
                         "MYXO.BC", "MYX.THEL", "CEST.COMP", "CEST.COMP", "CEST.COMP", "CEST.PEA", "CEST.LVS",
                         "MYXO.SBAD", "MYXO.UNK",
                         "MYXO.UP", "TREM.MET", "TREM.SSS", "TREM.BUC",
                         "TREM.METAF", "TREM.METAGS", "TREM.MS", "TREM.NEA",
                         "NEM.CYST", "TREM.META.GG", "TREM.META.UNK", "TREM.LG",
                         "MYX.TAIL", "MYX.GEOM", "NEM.CYST", "NEM.DAFT",
                         "TREM.CM", "TREM.META.HET", "TREM.META.SP", "TREM.BC",
                         "TREM.META.GO", "TREM.META.ES", "MYX.AK", "MYX.FI",
                         "MYX.TIN", "MYX.GL", "MYX.OT", "MYX.EYE",
                         "MYX.GC", "NEM.CYST", "TREM.LARV", "TREM.META",
                         "MYX.SP", "MYX.ROUND", "TREM.CYST") ~ "larval",
      
      is.na(psite_spp_x) ~ "0",
    ),
    
    psite_taxon = case_when(
      psite_count == 0 ~ "uninfected",  
      startsWith(psite_spp_x, "CEST") ~ "cest",
      startsWith(psite_spp_x, "TREM") ~ "trem",
      startsWith(psite_spp_x, "META") ~ "trem",
      startsWith(psite_spp_x, "NEM") ~ "nem",
      TRUE ~ "other"
    )
  ) %>%
  filter(!str_detect(psite_spp_x, "^(ACAN|MONO|COPE|LEECH|NMORPH|MYX)")) %>%
  replace_na(list(life_stage = "uninfected")) %>% 
  mutate(weight_g = (weight_mg/1000)) %>%
  janitor::clean_names()

### FCF added back
full.data.taxon <- full.data.taxon %>%
  mutate(FCF = ((((full.data.taxon$weight_g)/1000)/((full.data.taxon$standard_length_mm)^3)))*(100)) %>% 
  janitor::clean_names()  # this makes a second row of fcf_2, but i am keeping this in here because I can't tell exactly where the first fcf column is coming from using this above code

### fix variable type
full.data.taxon <- full.data.taxon %>%
  mutate(life_stage = as.factor(life_stage))
full.data.taxon <- full.data.taxon %>%
  mutate(life_stage = relevel(life_stage,ref='uninfected'))
full.data.taxon <- full.data.taxon %>%
  mutate(psite_taxon = as.factor(psite_taxon))

## Plots for Question 3
# FCF ~ life_stage * taxon + psite_count (A: T vs C, L: T vs C)
# Jolee's plots - messing around
ggplot(full.data.revised %>%
         filter(psite_taxon != "CEST"),
       aes(x = life_stage,
           y = fcf,
           colour = psite_taxon,
           scale_colour_gradient(psite_count, low = "blue", high = "red", name = "Parasite Count"))) + # fix to add in psite_count possibly
  geom_boxplot() +
  xlab("Life Stage") +
  ylab("FCF") +
  labs(colour = "Parasite Taxon") +
  theme_bw()



## Model formula for Question 3 
# K_i\sim\ life\ stage\ * taxon\ + parasite\ count\ + ci + (1|fish)\ + (1|catelog) 

#### troubleshooting
unique_fish_fcf <- full.data.stages %>%
  select(individual_fish_id, fish_sp_x, standard_length_mm, weight_g, weight_mg, fcf) %>%
  distinct()

view(unique_fish_fcf)

ggplot(full.data.stages, 
       aes(x = standard_length_mm, y = weight_g)) + geom_point() + geom_smooth(method = lm) + facet_wrap(~fish_sp_x) + geom_jitter()

fdsr <- full.data.stages %>%
  mutate(ratio = (weight_mg/standard_length_mm))

ggplot(fdsr, 
       aes(x = fish_sp_x, y = ratio)) + geom_boxplot() + geom_smooth(method = lm)
plot(fdsr$ratio)

