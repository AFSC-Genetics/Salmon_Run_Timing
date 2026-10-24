# prepare run timing data for SRA upload
# animal.xlsx

library(tidyverse)
library(readr)

#setwd("C:/Users/Natasha.Howe/Work/Salmon-Genetics/Salmon_runtiming")

# pat found some spots where I didn't check the metadata
# in part because I used an old ABLG exports. 
# He added updates to this file (10/23/2026), so let's work off of it
meta <- read_csv("./data/raw/SRA_metadata_preparation_Updated.csv")

# need to add unique identifier outside of sample_name, so will add AltID
# let's use a more up to date export for these samples
export <- readxl::read_xlsx("../20250623_export.xlsx")
head(export)

# format is already there, so I mostly just need to change the AlternateIDs selectively
sra <- meta %>%
  mutate(ABLG = gsub('ABLG','',sample_name) %>% as.numeric()) 

max(sra$ABLG)
sra %>%
  left_join(., export %>%
              select(ABLG, AlternateID) %>%
              rename(AlternateIDex = AlternateID),by="ABLG") %>% 
  # this is very specific but there's reason 
  # bc pat didn't edit those without sex
  mutate(AlternateID = case_when(AlternateID == AlternateIDex ~ AlternateID,
                                 sex == "not collected" ~ AlternateIDex,
                                 T ~ AlternateID)) %>%
  select(-ABLG, -AlternateIDex) %>%
  write.csv("./SRA/SRA_metadata_preparation_Updated_PB_NH_20261023.csv",
          row.names = F, quote = F)
