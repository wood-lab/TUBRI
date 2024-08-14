
## for picking out total parasite counts per dissector

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

#### Load packages
library(ggeffects)
library(tidyverse)
library(ggpubr)
library(stringr)
library(tibble)
library(car)
library(tidyr)
library(readr)
library(dbplyr)
library(ggplot2)

## load datasets

pimvig <- read.csv("data/processed/Pimephales_vigilax_processed_human_readable_UPDATED_2024.08.08.csv") %>%
  mutate(
    combined_dissector = coalesce(!!!select(., starts_with("dissector")))
  ) %>%
  janitor::clean_names()

ictpun <- read.csv("data/processed/Ictalurus_punctatus_processed_human_readable.csv_UPDATED_2024.08.01.csv") %>%
  mutate(
    combined_dissector = coalesce(!!!select(., starts_with("dissector")))
  ) %>%
  janitor::clean_names()

notath <- read.csv("data/processed/Notropis_atherinoides_processed_human_readable_UPDATED_2024.08.01.csv") %>%
  mutate(
    combined_dissector = coalesce(!!!select(., starts_with("dissector")))
  ) %>%
  janitor::clean_names()

hybnuc <- read_csv("data/processed/Hybognathus_nuchalis_processed_human_readable_UPDATED_2024.08.01.csv") %>%
  mutate(
    combined_dissector = coalesce(!!!select(., starts_with("dissector")))
  ) %>%
  janitor::clean_names()

pervig <- read_csv("data/processed/Percina_vigil_processed_human_readable_UPDATED_2024.08.01.csv") %>%
  mutate(
    combined_dissector = coalesce(!!!select(., starts_with("dissector")))
  ) %>%
  janitor::clean_names()

carvel <- read_csv("data/processed/Carpiodes_velifer_processed_human_readable_UPDATED_2024.08.08.csv") %>%
  mutate(
    combined_dissector = coalesce(!!!select(., starts_with("dissector")))
  ) %>%
  janitor::clean_names()


# Define the initials and their corresponding names
initials <- c("CJW", "GMC", "DMDM", "KL", "CLW", "JT", "IJ", "SC", "DJB")
names_map <- c("CJW" = "connor", "GMC" = "gabby", "DMDM" = "daki", "KL" = "katie", "CLW" = "chelsea", "JT" = "jolee", 
               "IJ" = "imani", "SC" = "shyanne", "DJB" = "desmond")

# Function to extract initials and choose the second if two are found
extract_initials <- function(text) {
  matches <- str_extract_all(text, regex(paste(initials, collapse = "|"), ignore_case = TRUE))[[1]]
  matches <- toupper(matches) # Ensure all matches are uppercase
  if (length(matches) == 0) {
    return(NA)
  } else if (length(matches) >= 2) {
    return(matches[2])
  } else {
    return(matches[1])
  }
}

## extract initials for each dataset

pimvig <- pimvig %>%
  mutate(
    extracted_initials = sapply(combined_dissector, extract_initials),
    person = names_map[extracted_initials]
  ) %>%
  janitor::clean_names()

ictpun <- ictpun %>%
  mutate(
    extracted_initials = sapply(combined_dissector, extract_initials),
    person = names_map[extracted_initials]
  ) %>%
  janitor::clean_names()

notath <- notath %>%
  mutate(
    extracted_initials = sapply(combined_dissector, extract_initials),
    person = names_map[extracted_initials]
  ) %>%
  janitor::clean_names()

hybnuc <- hybnuc %>%
  mutate(
    extracted_initials = sapply(combined_dissector, extract_initials),
    person = names_map[extracted_initials]
  ) %>%
  janitor::clean_names()

pervig <- pervig %>%
  mutate(
    extracted_initials = sapply(combined_dissector, extract_initials),
    person = names_map[extracted_initials]
  ) %>%
  janitor::clean_names()

carvel <- carvel %>%
  mutate(
    extracted_initials = sapply(combined_dissector, extract_initials),
    person = names_map[extracted_initials]
  ) %>%
  janitor::clean_names()

## define and modulate expected columns

expected_columns <- c(
  "acanth_bigb", "acanth_bkr", "nem_cona", "nem_cap",
  "nem_dich", "nem_larv", "nem_trbk", "nem_unk",
  "cest_godz", "cest_ntg", "cope_cl", "cope_grab",
  "cest_l", "cest_vit", "cest_unk", "mono_gyro",
  "mono_dact", "mono_lg", "trem_rnda", "trem_unk",
  "acan_ostr", "cest_meg", "leech_kut", "nem_nmts",
  "mono_ip", "mono_unk", "nem_phar", "nem_pty",
  "nem_sp", "nem_taf", "nmorph", "trem_2",
  "trem_allo", "trem_ble", "trem_crep", "trem_f",
  "trem_lip", "trem_megict", "trem_smile", "cest_hyd",
  "cest_tri", "cope_a", "mono_all", "nem_bran",
  "nem_spky", "nem_tfs", "nem_twas", "nem_unkn",
  "trem_gb", "trem_l", "cest_bb", "trem_ass",
  "trem_ns", "trem_diplo", "mono_gdac", "acan_ad",
  "cope_pood", "cest_ws", "cest_meagan", "cest_arch",
  "cest_unka", "cope_hny", "nem_cerl", "trem_sep",
  "trem_ors", "trem_cv", "mono_nid", "nem_heath",
  "men_squig", "acanth_bb", "trem_phlb", "cest_bot",
  "meta_hs", "meta_unk", "myx_go", "myx_geom",
  "myxo_bc", "myx_thel", "cest_comp", "cest_pea",
  "cest_lvs", "myxo_sbad", "myxo_unk", "myxo_up",
  "trem_met", "trem_sss", "trem_buc", "trem_metaf",
  "trem_metags", "trem_ms", "trem_nea", "nem_cyst",
  "trem_meta_gg", "trem_meta_unk", "myx_tail",
  "nem_daft", "trem_cm", "trem_meta_het", "trem_meta_sp",
  "trem_bc", "trem_meta_go", "trem_meta_es", "myx_ak",
  "myx_fi", "myx_tin", "myx_gl", "myx_ot",
  "myx_eye", "myx_gc", "trem_larv", "trem_meta",
  "myx_sp", "myx_round", "trem_cyst", "meta_unk",
  "myx_s", "myx_f", "myx_e", "myx_g", "myx_cm",
  "myx_bna", "myx_kl", "myx_gm", "acanth_cya",
  "acanth_thin", "trem_sk", "trem_chonk", "trem_endous",
  "trem_cl", "trem_ispot"
)

## filter to include only existing columns per dataset

existing_columns_pimvig <- intersect(expected_columns, colnames(pimvig))

existing_columns_ictpun <- intersect(expected_columns, colnames(ictpun))

existing_columns_notath <- intersect(expected_columns, colnames(notath))

existing_columns_hybnuc <- intersect(expected_columns, colnames(hybnuc))

existing_columns_pervig <- intersect(expected_columns, colnames(pervig))

existing_columns_carvel <- intersect(expected_columns, colnames(carvel))


## select each dataset to group by the dissector

tc_pimvig <- pimvig %>%
  select(person, all_of(existing_columns_pimvig)) %>%
  pivot_longer(cols = all_of(existing_columns_pimvig), names_to = "parasite_type", values_to = "count") %>%
  group_by(person, parasite_type) %>%
  summarize(total_count = sum(count, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = "parasite_type", values_from = "total_count", values_fill = list(total_count = 0))

tc_ictpun <- ictpun %>%
  select(person, all_of(existing_columns_ictpun)) %>%
  pivot_longer(cols = all_of(existing_columns_ictpun), names_to = "parasite_type", values_to = "count") %>%
  group_by(person, parasite_type) %>%
  summarize(total_count = sum(count, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = "parasite_type", values_from = "total_count", values_fill = list(total_count = 0))

tc_notath <- notath %>%
  select(person, all_of(existing_columns_notath)) %>%
  pivot_longer(cols = all_of(existing_columns_notath), names_to = "parasite_type", values_to = "count") %>%
  group_by(person, parasite_type) %>%
  summarize(total_count = sum(count, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = "parasite_type", values_from = "total_count", values_fill = list(total_count = 0))

tc_hybnuc <- hybnuc %>%
  select(person, all_of(existing_columns_hybnuc)) %>%
  pivot_longer(cols = all_of(existing_columns_hybnuc), names_to = "parasite_type", values_to = "count") %>%
  group_by(person, parasite_type) %>%
  summarize(total_count = sum(count, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = "parasite_type", values_from = "total_count", values_fill = list(total_count = 0))

tc_pervig <- pervig %>%
  select(person, all_of(existing_columns_pervig)) %>%
  pivot_longer(cols = all_of(existing_columns_pervig), names_to = "parasite_type", values_to = "count") %>%
  group_by(person, parasite_type) %>%
  summarize(total_count = sum(count, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = "parasite_type", values_from = "total_count", values_fill = list(total_count = 0))

tc_carvel <- carvel %>%
  select(person, all_of(existing_columns_carvel)) %>%
  pivot_longer(cols = all_of(existing_columns_carvel), names_to = "parasite_type", values_to = "count") %>%
  group_by(person, parasite_type) %>%
  summarize(total_count = sum(count, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = "parasite_type", values_from = "total_count", values_fill = list(total_count = 0))

## view the tables individually

#view(tc_pimvig) 
#view(tc_ictpun) 
#view(tc_notath)
#view(tc_hybnuc)
#view(tc_pervig)
#view(tc_carvel)

## combine the tables

combined_data <- tc_pimvig %>%
  full_join(tc_ictpun, by = "person") %>%
  full_join(tc_notath, by = "person") %>%
  full_join(tc_hybnuc, by = "person") %>%
  full_join(tc_pervig, by = "person") %>%
  full_join(tc_carvel, by = "person") %>%
  janitor::clean_names()


## pivot long
long_data <- combined_data %>%
  pivot_longer(cols = -person, names_to = "variable", values_to = "value")

## pivot wide
wide_data <- long_data %>%
  pivot_wider(names_from = person, values_from = value)

view(wide_data)


## totals:
#katie myx.g, trem.cm, trem.metags
#daki myx.g, trem.cm, mono.dact
#gabby myx.g, trem.cm, mono.dact 
#connor myx.g, trem.cm, mono.dact 
#chelsea trem.cm, mono.all, mono.dact
