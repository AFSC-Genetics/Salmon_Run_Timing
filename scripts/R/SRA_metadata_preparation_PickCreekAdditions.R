# prepare run timing data for SRA upload
# animal.xlsx

library(tidyverse)
library(readr)

setwd("C:/Users/Natasha.Howe/Work/Salmon-Genetics/Salmon_runtiming")

OUTDIR<-"../Salmon_Run_Timing/SRA/"

pickmeta <- read_csv("./2024_sockeye/data/raw/pick_creek_metadata.csv")

# will need to add unique identifier outside of sample_name, so will add AltID
export <- readxl::read_xlsx("./2024_sockeye/data/raw/pick_creek_ABLGexports.xlsx")
head(export)

sra <- pickmeta %>%
  mutate(sample = str_c('ABLG',sample)) %>%
  select(sample,ABLG,Runtime,CollectionDate,StartLatitudeDD,StartLongitudeDD)%>%
  mutate(organism = "Oncorhynchus nerka",
         ecotype = str_c(Runtime," Stream"),
         dev_stage = "adult",
         collection_date = CollectionDate,
         geo_loc_name = "USA: Alaska, Bristol Bay, Pick Creek",
         sex = "not collected",
         tissue = "axillary process",
         lat_lon = str_c(StartLatitudeDD," N ",-1*StartLongitudeDD," W")) %>%
  rename(sample_name = sample) %>%
  select(-c(Runtime,CollectionDate,StartLatitudeDD,StartLongitudeDD))

# make sure these values match
sra %>%
  count(ecotype)

# need an extra identifying column outside of sample_name
sra_id <- sra %>%
  left_join(., export %>% select(ABLG,AlternateID)) %>%
  rename(alternateID = AlternateID)

sra_id %>%
  write_csv(str_c(OUTDIR,"SRA_metadata_preparation_PickCreekAdditions.csv"))

#### SEQUENCING METADATA #####################################################

fastqs <- read.delim2("../Salmon_Run_Timing/data/R/pick_creek_runtiming_fastq_filenames.txt",
                      header = F) %>%
  mutate(sample_name = str_split(V1, "_", simplify = T)[,1],
         sample_name = ifelse(startsWith(sample_name,"ABLG"),sample_name,paste0('ABLG',sample_name)),
         R = sub(".*_L\\d+_([^_]+).*", "\\1", V1),
         R = sub(".fq.gz","",R), R = sub("R","",R)) %>%
  pivot_wider(id_cols = sample_name, names_from = R, values_from = V1) %>%
  left_join(sra_id) %>%
  group_by(ecotype) %>%
  mutate(title = str_c('Sockeye_',sub(' ','_',ecotype),'_Pick_Creek_',row_number())) %>%
  ungroup() %>%
  select(sample_name, alternateID, title, `1`, `2`)

write_csv(fastqs, str_c(OUTDIR,"SRA_sequencing_preparation_PickCreekAdditions.csv"))

### BASH COMMANDS ####
# srun -c 5 -t 48:00:00 --pty /bin/bash
# module load tools/lftp/4.9.2
# lftp subftp@ftp-private.ncbi.nlm.nih.gov # enter password given on SRA
# cd /uploads/natashahowe15_gmail.com_ua4ygaIi/SUB15992991_pick_runtiming
# mirror -Rv --continue --parallel=5 --use-pget-n=2 /scratch2/nhowe/runtiming/sock-pick/sednagold/Oncorhynchus_nerka/sockeye0425/ .