# prepare run timing data for SRA upload
# animal.xlsx

library(tidyverse)
library(readr)

setwd("C:/Users/Natasha.Howe/Work/Salmon-Genetics/Salmon_runtiming")

meta <- read_csv("C:/Users/Natasha.Howe/Work/Salmon-Genetics/Salmon_Run_Timing/data/raw/fourspecies_runtiming_metadata_forSRA.csv")

sra <- meta %>%
  filter(startsWith(sampleID,"ABLG")) %>% # remove euclide sockeye
  mutate(organism = case_when(Species == 'Sockeye' ~ "Oncorhynchus nerka",
                              Species == 'Pink' ~ "Oncorhynchus gorbuscha",
                              Species == 'Coho' ~ "Oncorhynchus kisutch",
                              Species == 'Chum' ~ "Oncorhynchus keta",
                              TRUE ~ "error"),
         geo_loc_name = case_when(Species == 'Sockeye' ~ "Bristol Bay - Whitefish Creek",
                                  Species == 'Pink' ~ "Juneau - Auke Creek",
                                  Species == 'Coho' & Ecotype == "Early" ~ "Slilak River - Inlet of Slikak Lake",
                                  Species == 'Coho' & Ecotype == "Late" ~ "Kenai River - Outlet of Slikak Lake",
                                  Species == 'Chum' & Ecotype == "Early" ~ "Yukon River - Gisasa River",
                                  Species == 'Chum' & Ecotype == "Late" ~ "Yukon River - Pelly River",
                                  TRUE ~ "error"),
         sex = case_when(Species == 'Sockeye' ~ "Male",
                         Species == 'Pink' ~ sex,
                         Species == 'Coho' ~ "not collected",
                         Species == 'Chum' ~ "not collected",
                         TRUE ~ "error"),
         collection_date = case_when(Species == 'Sockeye' ~ "2024-08-18",
                                     Species == 'Pink' & Ecotype == "Early" ~ paste0(CollectionDate,"-08"),
                                     Species == 'Pink' & Ecotype == "Late" ~ paste0(CollectionDate,"-09"),
                                     Species == 'Coho' & Ecotype == "Early" ~ "2003-09-24",
                                     Species == 'Coho' & Ecotype == "Late" ~ "2014-01-25",
                                     Species == 'Chum' & Ecotype == "Early" ~ "2004-06-30",
                                     Species == 'Chum' & Ecotype == "Late" ~ "1993-07-01",
                                     TRUE ~ NA),
         lat_lon = case_when(Species == 'Sockeye' ~ "59.267 N 158.686 W",
                             Species == 'Pink' ~ "58.381 N 134.641 W",
                             Species == 'Coho' & Ecotype == "Early" ~ "60.392 N 150.126 W",
                             Species == 'Coho' & Ecotype == "Late" ~ "60.464 N 150.531 W",
                             Species == 'Chum' & Ecotype == "Early" ~ "65.267 N 157.667 W",
                             Species == 'Chum' & Ecotype == "Late" ~ "62.793 N 137.337 W",
                             TRUE ~ NA)) %>%
  mutate(Lineage = ifelse(is.na(Lineage), "not applicable", paste(Lineage, 'Lineage'))) %>%
  rename(sample_name = sampleID, ecotype = Ecotype) %>%
  select(-CollectionDate)

# make sure these values match
sra %>%
  group_by(organism, Lineage, ecotype) %>%
  summarize(n=n())

# need to add unique identifier outside of sample_name, so will add AltID
export <- readxl::read_xlsx("../../MarineFisheries/ak-rockfish/data/20240411_export.xlsx")
head(export)

sra_id <- sra %>%
  mutate(ABLG = as.numeric(gsub("ABLG","",sample_name))) %>%
  left_join(., export) %>%
  mutate(AlternateID = ifelse(is.na(AlternateID),paste0(Species,'_',row_number()),AlternateID),
         dev_stage = 'adult',
         tissue = ifelse(Species == "Pink", 'axillary process', 'fin')) 

outfile <- sra_id %>%
  select(sample_name, organism, ecotype, dev_stage, collection_date, geo_loc_name, sex, tissue, AlternateID, lat_lon, Lineage)

write.csv(outfile, "C:/Users/Natasha.Howe/Work/Salmon-Genetics/Salmon_Run_Timing/SRA/SRA_metadata_preparation.csv",
          row.names = F, quote = F)

#### SEQUENCING METADATA #####################################################

fastqs <- read.delim2("../Salmon_Run_Timing/data/R/fourspecies_runtiming_fastq_filenames.txt",
                      header = F) %>%
  mutate(sample_name = str_split(V1, "_", simplify = T)[,1],
         sample_name = ifelse(startsWith(sample_name,"ABLG"),sample_name,paste0('ABLG',sample_name)),
         R = sub(".*_L\\d+_([^_]+).*", "\\1", V1),
         R = sub(".fq.gz","",R), R = sub("R","",R)) %>%
  pivot_wider(id_cols = sample_name, names_from = R, values_from = V1) %>%
  left_join(sra_id) %>%
  group_by(Species, ecotype) %>%
  mutate(title = case_when(Species != 'Pink' ~ paste0(Species,'_',ecotype,'_',row_number()),
                           TRUE ~ paste0(Species,'_',gsub(' Lineage','',Lineage),'_',ecotype,'_',row_number()))) %>%
  ungroup() %>%
  select(sample_name, AlternateID, title, `1`, `2`)

write.csv(fastqs, "C:/Users/Natasha.Howe/Work/Salmon-Genetics/Salmon_Run_Timing/SRA/SRA_sequencing_preparation.csv",
          row.names = F, quote = F)
