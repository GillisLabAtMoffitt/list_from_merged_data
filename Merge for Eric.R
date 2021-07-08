library(tidyverse)

path <- fs::path("","Volumes","Gillis_Research", "Christelle Colin-Leitzinger", "merging slid ID")

Eric_slid <-
  readxl::read_xlsx(paste0(path, "/metadata_myeloma_r045.xlsx")) %>% 
  # select(-X1) %>% 
  select("SLID.rnaseq")


# Load updated tumor dna and rna
dna_rna_slid <-
  readxl::read_xlsx(paste0(path, "/Avatar_MM_03162021_OUT.xlsx")) %>% 
  select("RNASequencingLibraryID", "disease_status") %>% 
  filter(!is.na(RNASequencingLibraryID))


# MERGE

Eric_DS <- left_join(Eric_slid, dna_rna_slid, by = c("SLID.rnaseq" = "RNASequencingLibraryID"))
write_csv(Eric_DS, paste0(path, "/merged RNA SLID and disease status at RNA collection.csv"))




# July 2021

path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "CHIP in Avatar")
tageted_seq <- readxl::read_xlsx(paste0(path, "/Nancy's working files/M4M_MM Avatar targeted_CH from 03.03.21 paired_for Christelle.xlsx")) %>% 
  select(patient_id) %>% 
  distinct(patient_id, .keep_all = TRUE)
tarseq_patients <- paste(tageted_seq$patient_id, collapse = "|")

germline_patient_data <- readRDS("/Users/colinccm/Documents/GitHub/Gillis/list_from_merged_data/germline_patient_data.rds")
germline_patient_data <- germline_patient_data %>% 
  filter(is_patient_MM == "Yes" & is_MMDx_close_to_blood == "Yes") %>% 
  distinct(avatar_id, .keep_all = TRUE)
MM_patients <- paste(germline_patient_data$avatar_id, collapse = "|")

path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "merging slid ID")
rna_dna <- read_csv(paste0(path, "/list RNA SLID corresponding to DNA SLID 03222021.csv")) %>% 
  select(avatar_id, SLID_germline, collectiondt_germline, SLID_tumor, collectiondt_tumor, collectiondt_dna, 
         collectiondt_rna, if_no_rna_in_tumor1_use_tumor_2 = interval_rna_vs_dna, SLID_tumor_2, collectiondt_rna_tumor_2) %>% 
  mutate(if_no_rna_in_tumor1_use_tumor_2 = str_remove(if_no_rna_in_tumor1_use_tumor_2, " date")) 

tarseq_patients1 <- rna_dna %>% 
  filter(str_detect(SLID_germline, tarseq_patients))
write_csv(tarseq_patients1, paste0(path, "/list targeted sequencing patients rna and dna 07082021.csv"))
MM_patients1 <- rna_dna %>% 
  filter(str_detect(avatar_id, MM_patients))
write_csv(MM_patients1, paste0(path, "/list MM patients rna and dna 07082021.csv"))






