# Packages
library(tidyverse)

path <- fs::path("","Volumes","Gillis_Research", "Christelle Colin-Leitzinger", "merging slid ID")


DNA_id <-
  read_csv(paste0(path, "/Gillis_MM Avatar_sample selection.csv")) %>% 
  select(-X1)

RNA_id <- 
  readxl::read_xlsx(paste0(path, "/Avatar_MM_Gillis_Linkage_11062020.xlsx")) %>% 
  select(c("avatar_id", "DNASequencingLibraryID", "RNASequencingLibraryID", "collectiondt"))



RNA_cleaned <- RNA_id %>% 
  group_by(avatar_id, collectiondt) %>% 
  fill(DNASequencingLibraryID, RNASequencingLibraryID) %>% 
  drop_na(RNASequencingLibraryID)


list <- right_join(RNA_cleaned, DNA_id, by = "avatar_id")

write_csv(list, paste0(path, "/list RNA SLID corresponding to DNA SLID.csv"))







