# Packages
library(tidyverse)
library(lubridate)

path <- fs::path("","Volumes","Gillis_Research", "Christelle Colin-Leitzinger", "merging slid ID")

# Import previous data to keep only the patient and tumor already running in the bioinformatician
Previous_germline <-
  read_csv(paste0(path, "/List tumor SLID earliest or closest to germline.csv")) %>% 
  # select(-X1) %>% 
  select(everything(), interval_germline_vs_first_tumor = "interv", -X1)
dna_selected <- paste0(Previous_germline$SLID_tumor, collapse = "|")

# Load updated tumor dna and rna
dna_rna_slid <-
  readxl::read_xlsx(paste0(path, "/Avatar_MM_03162021_OUT.xlsx"))

# Clean
dna_rna_slid2 <- dna_rna_slid %>% filter(!str_detect(DiseaseType, "germline")) %>% 
  select(avatar_id, collectiondt, SLID_tumor = "DNASequencingLibraryID", "RNASequencingLibraryID",
         disease_status_at_tumor = disease_status) %>% 
  mutate(SLID_tumor = str_replace(SLID_tumor, "unavailable", NA_character_)) %>% 
  mutate(RNASequencingLibraryID = str_replace(RNASequencingLibraryID, "unavailable", NA_character_)) %>% 
  group_by(avatar_id, collectiondt) %>% 
  fill(SLID_tumor, RNASequencingLibraryID, .direction = "updown") %>% 
  ungroup() %>% 
  distinct() %>% 
  # For RNA and DNA collected at the same time
  # Store collectiondt_dna for the tumor in the bioinf list
  mutate(selected_dna_tumor = ifelse(str_detect(SLID_tumor, dna_selected), "Yes", NA_character_)) %>% 
  mutate(collectiondt_dna = case_when(
    selected_dna_tumor == "Yes"            ~ collectiondt, 
    TRUE                                   ~ NA_POSIXct_
  )) %>% 
  # keep RNA only when selected tumor is associated
  mutate(rna_slid = ifelse(selected_dna_tumor == "Yes", RNASequencingLibraryID, NA_character_)) %>% 
  mutate(collectiondt_rna = case_when(
    !is.na(rna_slid)                       ~ collectiondt, 
    TRUE                                   ~ NA_POSIXct_
    )) %>% 
  
  # Calculate interval for different date of collection
  group_by(avatar_id) %>% 
  fill(collectiondt_dna, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(interval_rna_vs_dna = ifelse(!is.na(RNASequencingLibraryID), abs((interval(start = collectiondt_dna, end = collectiondt))/
           duration(n=1, units = "days")), NA)
         ) %>% 
  # Remove the date filling
  mutate(collectiondt_dna = case_when(
    collectiondt_dna != collectiondt          ~ NA_POSIXct_,
    TRUE                                      ~ collectiondt_dna
  )) %>% 
  # Pick the closest RNA date to DNA
  arrange(avatar_id, interval_rna_vs_dna) %>% 
  group_by(avatar_id) %>% 
  mutate(closest = dense_rank(interval_rna_vs_dna)) %>% 
  mutate(rna_slid = case_when(
    is.na(rna_slid) &
      closest == 1               ~ RNASequencingLibraryID,
    TRUE                         ~ rna_slid
  )) %>% 
  mutate(collectiondt_rna = case_when(
    is.na(collectiondt_rna) &
      closest == 1               ~ collectiondt,
    TRUE                         ~ collectiondt_rna
  )) %>% 
  # Just quick code to look at the new germline
  mutate(selected_patient = ifelse(str_detect(SLID_tumor, dna_selected), "Yes", NA_character_)) %>% 
  fill(selected_patient, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(selected_patient = ifelse(is.na(selected_patient), "No", selected_patient)) %>% 
  #
  filter((selected_dna_tumor == "Yes") |
           !is.na(rna_slid) |
           # (is.na(selected_dna_tumor) & (rna_slid == RNASequencingLibraryID)) |
           # (!is.na(selected_dna_tumor) & is.na(rna_slid) & is.na(interval_rna_vs_dna)) |
           selected_patient == "No"
         ) %>% 
  # group_by(avatar_id) %>% 
  # fill(rna_slid, collectiondt_rna, interval_rna_vs_dna, .direction = "updown") %>% 
  # ungroup() %>% 
  
  # Keep the DNA rows, rows with RNA alone
  select(-closest, -RNASequencingLibraryID) %>% 
  arrange(avatar_id, collectiondt)

dna_slid <- dna_rna_slid2 %>% 
  filter(selected_dna_tumor == "Yes") %>% 
  mutate(interval_rna_vs_dna = ifelse(is.na(interval_rna_vs_dna), "no rna", "same date")) %>% 
  select(-c(collectiondt, selected_patient))


rna_slid <- dna_rna_slid2 %>% 
  filter(is.na(selected_dna_tumor) & selected_patient == "Yes") %>% 
  mutate(selected_dna_tumor = ifelse(is.na(collectiondt_dna), "other tumor", "same date other tumor")) %>% 
  select(avatar_id, 
         SLID_tumor_2 = "SLID_tumor", disease_status_at_tumor_2 = "disease_status_at_tumor", 
         type_of_second_tumor = "selected_dna_tumor", collectiondt_dna_tumor_2 = "collectiondt_dna",
         rna_slid_tumor_2 = "rna_slid", collectiondt_rna_tumor_2 = "collectiondt_rna",
         second_tumor_interval_rna_vs_dna = "interval_rna_vs_dna")

current_germline <- dna_rna_slid %>% filter(str_detect(DiseaseType, "germline")) %>%
  select(avatar_id, collectiondt_germline = "collectiondt", SLID_germline = "DNASequencingLibraryID",
         disease_status_at_germline = "disease_status") %>%
  mutate(SLID_germline = str_replace(SLID_germline, "unavailable", NA_character_))

n_data <- full_join(dna_slid, rna_slid, by = "avatar_id")




new_data <- full_join(Previous_germline, n_data, by = c("avatar_id", "SLID_tumor")) %>% 
  left_join(., current_germline %>% select(avatar_id, SLID_germline, disease_status_at_germline),
            by = c("avatar_id", "SLID_germline")) %>% 
  select("avatar_id", "SLID_germline", "collectiondt_germline", "disease_status_at_germline", everything())
# data$SLID_tumor[(data$avatar_id == "A000180" & data$germline_slid == "SL216867")] <- "not requested"


# new_data <- left_join(Previous_germline, data, by = c("avatar_id", "SLID_germline", "SLID_tumor"))
new_data$disease_status_at_germline[(data$avatar_id == "A000428")] <- "not MM"
new_data$disease_status_at_germline[(data$avatar_id == "A000456")] <- "not MM"


write_csv(new_data, paste0(path, "/list RNA SLID corresponding to DNA SLID 03222021.csv"))

new_germline <- dna_rna_slid2 %>% filter(selected_patient == "No") %>% select(avatar_id) %>% distinct()
write_csv(new_germline, paste0(path, "/list new germline ID.csv"))



# germline$avatar_id[duplicated(germline$avatar_id)]
# [1] "A000180" "A000180" "A010797" "A010798" "A010799" "A020520" "A022605"

