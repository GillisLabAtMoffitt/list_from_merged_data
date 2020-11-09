# Packages
library(tidyverse)

path <- fs::path("","Volumes","Gillis_Research", "Christelle Colin-Leitzinger", "merging slid ID")


DNA_id <-
  read_csv(paste0(path, "/Gillis_MM Avatar_sample selection.csv")) %>% 
  select(-X1)

RNA_id <- 
  readxl::read_xlsx(paste0(path, "/Avatar_MM_Gillis_Linkage_11062020.xlsx")) %>% 
  select(c("avatar_id", "DNASequencingLibraryID", "RNASequencingLibraryID", "collectiondt")) %>% 
  mutate(RNASequencingLibraryID = str_replace(RNASequencingLibraryID, "unavailable", NA_character_))



RNA_cleaned <- RNA_id %>% 
  group_by(avatar_id, collectiondt) %>% 
  fill(DNASequencingLibraryID, RNASequencingLibraryID, .direction = "downup") %>% 
  fill(RNASequencingLibraryID, DNASequencingLibraryID, .direction = "downup") # Just to make sure but doesn't do anything


list <- right_join(RNA_cleaned, DNA_id, by = c("avatar_id", "DNASequencingLibraryID" = "SLID_tumor")) %>% 
  distinct()

# write_csv(list, paste0(path, "/list RNA SLID corresponding to DNA SLID.csv"))

# Making sure that the NA in the list are NA in the RNA df
ID <- list %>% 
  ungroup() %>% 
  filter(is.na(RNASequencingLibraryID)) %>% 
  select(avatar_id)

ID <- left_join(ID, RNA_id, by = "avatar_id")


# Making sure that the NA in the list are unavailable in the RNA df

ID1 <- list %>% 
  ungroup() %>% 
  filter(RNASequencingLibraryID == "unavailable") %>% 
  select(avatar_id)

ID1 <- left_join(ID1, RNA_id, by = "avatar_id")


