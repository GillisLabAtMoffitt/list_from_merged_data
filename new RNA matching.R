# Packages
library(tidyverse)
library(lubridate)

path <- fs::path("","Volumes","Gillis_Research", "Christelle Colin-Leitzinger", "merging slid ID")

data <-
  read_csv(paste0(path, "/List tumor SLID earliest or closest to germline.csv")) %>% 
  select(-X1)
dna_selected <- paste0(data$SLID_tumor, collapse = "|")
  
dna_rna_slid <-
  readxl::read_xlsx(paste0(path, "/Avatar_MM_03162021_OUT.xlsx"))

dna_rna_slid2 <- dna_rna_slid %>% filter(!str_detect(DiseaseType, "germline")) %>% 
  select(avatar_id, collectiondt, SLID_tumor = "DNASequencingLibraryID", "RNASequencingLibraryID",
         disease_status_at_tumor = disease_status) %>% 
  mutate(SLID_tumor = str_replace(SLID_tumor, "unavailable", NA_character_)) %>% 
  mutate(RNASequencingLibraryID = str_replace(RNASequencingLibraryID, "unavailable", NA_character_)) %>% 
  group_by(avatar_id, collectiondt) %>% 
  fill(SLID_tumor, RNASequencingLibraryID, .direction = "updown") %>% 
  ungroup() %>% 
  distinct() %>% 
  mutate(selected_dna_tumor = ifelse(str_detect(SLID_tumor, dna_selected), "Yes", NA_character_)) %>% 
  mutate(collectiondt_dna = case_when(
    selected_dna_tumor == "Yes"            ~ collectiondt, 
    TRUE                                   ~ NA_POSIXct_
  )) %>% 
  # keep only slid rna when selected tumor is associated
  mutate(rna_slid = ifelse(selected_dna_tumor == "Yes", RNASequencingLibraryID, NA_character_)) %>% 
  mutate(collectiondt_rna = case_when(
    !is.na(rna_slid)                       ~ collectiondt, 
    TRUE                                   ~ NA_POSIXct_
    )) %>% 
  group_by(avatar_id) %>% 
  fill(rna_slid, collectiondt_dna, collectiondt_rna, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(interval_rna_vs_dna = ifelse(!is.na(RNASequencingLibraryID), (interval(start = collectiondt_dna, end = collectiondt))/
           duration(n=1, units = "days"), NA),
         interval_rna_vs_dna = ifelse(!is.na(rna_slid), 0, interval_rna_vs_dna)) %>% 
  
  mutate(collectiondt_dna = case_when(
    is.na(selected_dna_tumor)           ~ NA_POSIXct_,
    TRUE                                ~ collectiondt_dna
  )) %>% 
  
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
  mutate(selected_patient = ifelse(str_detect(SLID_tumor, dna_selected), "Yes", NA_character_)) %>% 
  fill(selected_patient, .direction = "updown") %>% 
  ungroup() %>% 
  mutate(selected_patient = ifelse(is.na(selected_patient), "No", selected_patient)) %>% 
  
  filter((selected_dna_tumor == "Yes") |
           (is.na(selected_dna_tumor) & (rna_slid == RNASequencingLibraryID)) |
           (!is.na(selected_dna_tumor) & is.na(rna_slid) & is.na(interval_rna_vs_dna)) |
           selected_patient == "No"
         ) %>% 
  select(-closest, -RNASequencingLibraryID) %>% 
  arrange(avatar_id, collectiondt)

germline <- dna_rna_slid %>% filter(str_detect(DiseaseType, "germline")) %>% 
  select(avatar_id, collectiondt_germline = "collectiondt", germline_slid = "DNASequencingLibraryID", 
         disease_status_at_germline = "disease_status") %>% 
  mutate(germline_slid = str_replace(germline_slid, "unavailable", NA_character_))

data <- full_join(germline, dna_rna_slid2, by = "avatar_id") %>% 
  arrange(avatar_id, collectiondt) # %>% 
data$SLID_tumor[(data$avatar_id == "A000180" & data$germline_slid == "SL216867")] <- "not requested"
  # group_by(avatar_id) %>% 
  # mutate(disease_status = factor(disease_status, levels = c(
  #   "Mgus", "Smoldering Multiple Myeloma", 
  #   "Pre Treatment Newly Diagnosed Multiple Myeloma", "Post Treatment Newly Diagnosed Multiple Myeloma",
  #   "Early Relapse Multiple Myeloma", "Late Relapse Multiple Myeloma",
  #   
  #   "Normal marrow", "POEMS", "Polyclonal Gammopathy", "Amyloidosis", "Amyloidosis- Diagnostic marrow",
  #   "MYELOFIBROSIS", "Solitary Plasmacytoma", "WALDENSTROM MACROGLOBULINEMIA",
  #   "Refractory anemia with ring sideroblasts"
  # ))) %>% 
  # mutate(progressed_disease_status = dense_rank(disease_status)) %>% 
  # ungroup() %>% 
  # mutate(progressed_disease_status = ifelse(progressed_disease_status == 1, "No", "Yes"))

# Problem A020520

write_csv(data, paste0(path, "/list RNA SLID corresponding to DNA SLID 03222021.csv"))

germline$avatar_id[duplicated(germline$avatar_id)]
[1] "A000180" "A000180" "A010797" "A010798" "A010799" "A020520" "A022605"

