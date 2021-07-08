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
# ------------> 192 patients
path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "CHIP in Avatar")
targeted_seq <- readxl::read_xlsx(paste0(path, "/Nancy's working files/M4M_MM Avatar targeted_CH from 03.03.21 paired_for Christelle.xlsx")) %>% 
  distinct(patient_id, .keep_all = TRUE)
tarseq_patients <- paste(targeted_seq$patient_id, collapse = "|")

# -----------> 648 patients (180 has 2 blood)
germline_patient_data <- readRDS("/Users/colinccm/Documents/GitHub/CHIP-Avatar/germline_patient_data.rds") %>% 
  distinct(avatar_id, SLID_germline, .keep_all = TRUE)

# -----------> 480 MM
MM_data <- germline_patient_data %>% 
  distinct(avatar_id, .keep_all = TRUE) %>% 
  filter(is_patient_MM == "Yes" & is_MMDx_close_to_blood == "Yes")

# are the 192 in the 647    --------> Yes
a <- germline_patient_data %>% 
  filter(str_detect(SLID_germline, tarseq_patients))

# are the 192 in the 647    --------> No only 167
b <- MM_data %>% 
  filter(str_detect(SLID_germline, tarseq_patients))


# bb <- MM_data %>% 
#   filter(str_detect(SLID_germline, tarseq_patients)) %>% 
#   distinct(avatar_id, .keep_all = TRUE)
bb <- paste(b$SLID_germline, collapse = "|")

bbb <- targeted_seq %>% # -------------> lost 25 targ seq that are not in MM
  filter(!str_detect(patient_id, bb))

lost <- paste(bbb$patient_id, collapse = "|")
lost1 <- germline_patient_data %>% 
  filter(str_detect(SLID_germline, lost))
lost2 <- lost1 %>% select(avatar_id, mrn, MRN, moffitt_sample_id_germline, date_of_MMonly_diagnosis, MM_date_dx, Disease_Status_germline, is_patient_MM,
                          disease_stage_1, disease_stage_2, disease_stage_3, disease_stage_4, smoldering_status, is_MMDx_close_to_blood, collectiondt_germline, CH_status_TS ) 

mrn <- readRDS("/Users/colinccm/Documents/GitHub/CHIP-Avatar/mrn.rds") %>% 
  drop_na()
data <- left_join(lost2, mrn, by = "avatar_id") %>% select(MRN.y, everything())
write.csv(data, paste0(path, "/lost because no diagnosis date.csv"))


MM_patients <- germline_patient_data %>% select(avatar_id)
MM_patients <- paste(germline_patient_data$avatar_id, collapse = "|")
MM_patients_slid <- paste(germline_patient_data$SLID_germline, collapse = "|")

# Are the 9 patients we removed because of date of diag wrong (never progressed) are in the sequenced data   ------> NO
id <- paste("A022604","A027407","A029244","A007364","A000238",
            "A000530","A007146","A010533","A016764", sep = "|")
are_9_in_tarseq <- tageted_seq %>% 
  filter(!str_detect(patient_id, id))

# Are the patients we removed because of date of diag to far from blood are in the sequenced data   ------>
blood_far_to_wes <- germline_patient_data %>% 
  filter(is_MMDx_close_to_blood == "No")

blood_far_to_wes <- germline_patient_data %>% 
  distinct(avatar_id, .keep_all = TRUE) %>% 
  filter(is_patient_MM == "Yes" & is_MMDx_close_to_blood == "No") %>% 
  select(SLID_germline)
too_far <- paste(blood_far_to_wes$SLID_germline, collapse = "|")

blood_far_to_seq <- tageted_seq %>% 
  filter(str_detect(patient_id, too_far))

are_blood_far_to_wes <- germline_patient_data %>% 
  filter(is_patient_MM == "No")





MM_patients_a <- paste(germline_patient_data$SLID_germline, collapse = "|")

MM_patients_aa <- germline_patient_data_a %>% 
  filter(str_detect(SLID_germline, tarseq_patients))

germline_patient_data_b <- germline_patient_data %>% 
  filter(!is.na(date_of_MMonly_diagnosis)) %>% 
  distinct(avatar_id, .keep_all = TRUE)

MM_patients_bb <- germline_patient_data_b %>% 
  filter(!str_detect(SLID_germline, tarseq_patients)) %>% 
  filter(!str_detect(avatar_id, MM_patients))




path <- fs::path("","Volumes","Gillis_Research","Christelle Colin-Leitzinger", "merging slid ID")
rna_dna <- read_csv(paste0(path, "/list RNA SLID corresponding to DNA SLID 03222021.csv")) %>% 
  select(avatar_id, SLID_germline, collectiondt_germline, SLID_tumor, collectiondt_tumor, collectiondt_dna, 
         rna_slid, collectiondt_rna)

tarseq_patients1 <- rna_dna %>% 
  filter(str_detect(SLID_germline, tarseq_patients))
write_csv(tarseq_patients1, paste0(path, "/list targeted sequencing patients rna and dna 07082021.csv"))
MM_patients1 <- rna_dna %>% 
  filter(str_detect(avatar_id, MM_patients))
write_csv(MM_patients1, paste0(path, "/list MM patients rna and dna 07082021.csv"))

tageted_seq2 <- tageted_seq %>% 
  filter(!str_detect(patient_id, MM_patients_slid)) %>% 
  distinct(patient_id)













