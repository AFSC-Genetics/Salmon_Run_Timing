# prepare run timing data for SRA upload
# animal.xlsx

library(tidyverse)
library(readr)

#OUTDIR<-"../Salmon_Run_Timing/SRA/"
PRJ="PRJNA1259009"

GTseqMeta <- read_tsv("./SRA/GTseq_meta.txt")

colnames(GTseqMeta)

GTseqMeta %>% distinct(sex)
GTseqMeta %>% distinct(ecotype)
GTseqMeta %>% distinct(geo_loc_name)
GTseqMeta %>% count(collection_date)

sra <- GTseqMeta %>%
  mutate(sex = case_when(is.na(sex) ~ 'not collected',
                         sex == 'f' ~ 'Female',
                         sex == 'm' ~ 'Male',
                         T ~ sex),
         ecotype = ifelse(is.na(ecotype),"not applicable",ecotype),
         lineage = ifelse(is.na(lineage),"not applicable",lineage),
         collection_date = ifelse(is.na(collection_date),"not collected",collection_date)) %>%
  # resolve issue with tab in some lat_lon
  separate(col = lat_lon, into = c("lat", "lon"), sep = "\t", extra = "merge") %>%
  mutate(lat_lon = ifelse(is.na(lon),lat,paste(lat,lon))) %>%
  dplyr::select(-lat,-lon) %>%
  relocate(lat_lon,.before=`lineage`)

sra %>%
  write_delim("./SRA/SRA_metadata_preparation_GTseqAdditions_SUB16623141.tsv", delim="\t")

#### SEQUENCING METADATA #####################################################

fastqs <- read_tsv("./SRA/Run_Timing_GTseq_fastq_filelist.txt", col_names = F) %>%
  filter(!grepl("__MACOSX/._",X1), # rm these hidden files
         !grepl("._",X1, fixed=TRUE)) %>%
  # ADFG data is single-end, ABL data is paired-end
  mutate(library_layout = ifelse(grepl("R1|R2",X1),"paired","single"),
         sample = gsub("./","",X1, fixed=TRUE) %>% gsub("GTseqSRA/","",., fixed=TRUE),
         sample_name = case_when(library_layout=="paired" ~ str_split(sample, "_", simplify = T)[,1],
                                 T ~ gsub('.fastq$','',sample,fixed=TRUE)))

,
         R = sub(".*_L\\d+_([^_]+).*", "\\1", X1),
         R = sub(".fastq.gz","",R), R = sub("R","",R)) %>%
  pivot_wider(id_cols = sample_name, names_from = R, values_from = V1) %>%
  left_join(sra) %>%
  group_by(organism) %>%
  mutate(title = str_c('Sockeye_',sub(' ','_',ecotype),'_Pick_Creek_',row_number())) %>%
  ungroup() %>%
  select(sample_name, alternateID, title, `1`, `2`)

write_csv(fastqs, str_c(OUTDIR,"SRA_sequencing_preparation_GTseqAdditions.csv"))

### BASH COMMANDS ####
# srun -c 5 -t 48:00:00 --pty /bin/bash
# module load tools/lftp/4.9.2
# lftp subftp@ftp-private.ncbi.nlm.nih.gov # enter password given on SRA
# cd /uploads/natashahowe15_gmail.com_ua4ygaIi/SUB15992991_pick_runtiming
# mirror -Rv --continue --parallel=5 --use-pget-n=2 /scratch2/nhowe/runtiming/sock-pick/sednagold/Oncorhynchus_nerka/sockeye0425/ .