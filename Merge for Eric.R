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
