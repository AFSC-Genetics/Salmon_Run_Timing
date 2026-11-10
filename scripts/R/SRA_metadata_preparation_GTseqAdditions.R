# ADD GTseq Data
# prepare run timing data for SRA upload
# animal.xlsx

# I REALIZED THAT SOME OF THESE SAMPLES ARE DUPLICATES... 
# e.g. they already have a biosample accession bc they were whole genome sequenced & now GTseq'ed
# So I am going to split the data into 
# -- (i) non-WGS gtseq data and 
# -- (ii) WGS gtseq data
# They will be under two separate submissions
# Well shoot... SRA only takes 1000 samples at a time and we have over 2k
# will have to break it up into additional separate submissions unfortunately


library(tidyverse)
library(readr)

PRJ="PRJNA1259009"

GTseqMeta <- read_tsv("./SRA/GTseq_meta.txt")

# decided to split by project, but
GTseqMeta %>% 
  count(Project)
# I'll combine Auke Pinks 2019 & Chum Yukon River bc below 1000

SRAmeta <- GTseqMeta %>%
  mutate(bioproject_accession = PRJ,
         sex = case_when(is.na(sex) ~ 'not collected',
                         sex == 'f' ~ 'Female',
                         sex == 'm' ~ 'Male',
                         T ~ sex),
         ecotype = ifelse(is.na(ecotype),"not applicable",ecotype),
         lineage = ifelse(is.na(lineage),"not applicable",lineage)) %>%
  # resolve issue with tab in some lat_lon
  separate(col = lat_lon, into = c("lat", "lon"), sep = "\t", extra = "merge") %>%
  mutate(lat_lon = ifelse(is.na(lon),lat,paste(lat,lon))) %>%
  dplyr::select(-lat,-lon) %>%
  relocate(lat_lon,.before=`lineage`)

# these pink samples were dropped from analyses bc no collection date recorded
dropNoDateData <- SRAmeta %>%
  filter(is.na(collection_date)) %>%
  mutate(sample_name = str_c("ABLG",sample_name)) %>% 
  pull(sample_name)

# FROM the SRA downloadable form w/ relevant cols for our data
# useful for later
metaCols <- c("sample_name","bioproject_accession","organism",
              "ecotype","dev_stage","collection_date","geo_loc_name",
              "sex","tissue","AlternateID","lat_lon","lineage","Project")

seqCols <- c("sample_name", "library_ID","title","library_strategy",
             "library_source","library_selection","library_layout","platform",
             "instrument_model","design_description","filetype","filename","filename2" )




#### SEQUENCING METADATA #####################################################
# these samples need to be dropped
DropDuplicates <- read_tsv("./SRA/DroppedJuvenileChumDuplicates.txt")

# I am just going to get rid of those fastqs and recreate the fastq filelist without them
DropDuplicates %>%
  select(Fastq) %>%
  write_tsv("./SRA/DroppedJuvenileChumDuplicates_fastq_filelist.txt",col_names = F)


# okay, after moving those files... 
fastqfiles <- read_tsv("./SRA/Run_Timing_GTseq_fastq_filelist_updated.txt", col_names = F) %>%
  filter(!grepl("unused_fastq",X1, fixed = TRUE), # the above duplicated files were moved here
         !grepl("txt$",X1)) %>%
  # ADFG data is single-end, ABL data is paired-end
  mutate(sample = basename(X1),
         sample = sub("\\q$", "q.gz",sample), # will have to compress
         library_layout = if_else(grepl("(_|\\.)R[12](_|\\.)", X1), "paired", "single"),
         sample_name = case_when(
           grepl("^P10",sample) ~ str_c(str_split(sample, "_", simplify = T)[,2],"_",str_split(sample, "_", simplify = T)[,3]),
           T ~ str_split(sample, "_", simplify = T)[,1])) %>%
  mutate(R = case_when(library_layout=="single" ~ '1',
                       T ~ sub(".*_L\\d+_([^_]+).*", "\\1", sample)),
         R = sub(".fastq.gz","",R) %>% sub("R","",.))

# now format for the sra
SRA <- fastqfiles %>%
  tidyr::pivot_wider(id_cols = sample_name, names_from = R, values_from = sample) %>%
  dplyr::full_join(SRAmeta, by="sample_name") %>%
  mutate(sample_name = case_when(grepl('^P',sample_name) | grepl('^G',sample_name) ~ sample_name, 
                                 T ~ str_c('ABLG',sample_name)),
         library_ID = AlternateID,
         title = str_c("GT-sequencing of ",organism,": "),
         title = case_when(lineage != "not applicable" & sex == "not collected" ~ str_c(title,lineage,'-lineage ',dev_stage),
                           lineage != "not applicable" ~ str_c(title,lineage,'-lineage ',dev_stage,' ',sex),
                           sex == "not collected" ~ str_c(title,dev_stage),
                           T ~ str_c(title,dev_stage,' ',sex)),
         library_strategy="AMPLICON", library_source="GENOMIC", library_selection="PCR",
         library_layout = ifelse(is.na(`2`),"paired","single"),
         platform = "ILLUMINA",
         instrument_model = ifelse(is.na(`2`),"NextSeq 500","Illumina MiSeq"),
         design_description = "short read GT-sequencing for run timing markers",
         filetype = "fastq",
         filename = `1`, filename2 = ifelse(is.na(`2`),'',`2`)
         ) %>%
  filter(!(sample_name %in% dropNoDateData)) %>%
  dplyr::select(-c(`1`,`2`))

SRA %>%
  select(-colnames(SRAmeta)[-1]) %>%# remove metadata cols except sample_name
  filter(sample_name %in% dropNoDateData) %>%
  select(sample_name, filename,filename2) %>%
  pivot_longer(!sample_name) %>%
  select(value) %>%
  write_tsv("./SRA/DropDatelessSamples_fastq_filelist.txt",col_names = F)

# Separate data by if it was whole genome sequencer or not. 
# I will first upload the data that was not.
# this will include all the ADFG data,
# and expectantly some of the sockeye locations that were also WGS, but a subset of samples

# gonna have to submit saying that they already have biosample_accessions
# this seems simple. No metedata upload, just the sequence data and I downloaded the accessions from the other submissions and will pull them in here
pick_acc <- read_tsv('SRA/attributes_SUB15992991_PickCreekAdditionsWGS.tsv') %>%
  select(accession,sample_name)
og_acc <- read_tsv('SRA/attributes_SUB15296226_OriginalWGSUpload.tsv') %>%
  select(accession,sample_name)


# pickWGS <- read_tsv("./SRA/Salmon_Run_Timing_SRA_Sequences_PickCreekAdditions_SUB15992991.tsv") %>%
#   pull(sample_name)
# originalWGS <- read_csv("./SRA/SRA_metadata_preparation.csv") %>%
#   pull(sample_name)

justGTseqSRA <- SRA %>%
  filter(!(sample_name %in% pick_acc$sample_name),   # 192 pick creek samples
         !(sample_name %in% og_acc$sample_name)) %>% # actually none of the same whitefish samples
  dplyr::group_split(., Project)

# im gonna combine the two lowest sample size groups
justGTseqSRA[[1]] <- bind_rows(justGTseqSRA[[1]], justGTseqSRA[[4]])
justGTseqSRA <- setNames(justGTseqSRA[-4], c("PinkSockeyeABL","YukonChum","RangewidePink"))

# Metadata Output
iwalk(justGTseqSRA, ~ .x %>%
        dplyr::select(all_of(metaCols)) %>%
        write_tsv(str_c("./SRA/Salmon_Run_Timing_metadata_GTseqAdditions_",.y, ".tsv")))
# Sequence Output
iwalk(justGTseqSRA, ~.x %>%
        dplyr::select(all_of(seqCols)) %>%
        write_tsv(str_c("./SRA/Salmon_Run_Timing_sequence_GTseqAdditions_",.y, ".tsv")))

# Now I will need to create the proper subfolders for each of these datasets
# and zip unzipped fastqs
fastqfiles %>% 
  mutate(fastqfile = basename(X1),
         sample_name=case_when(grepl('^P',sample_name) | grepl('^G',sample_name) ~ sample_name,
                               T ~ str_c('ABLG',sample_name))) %>%
  dplyr::select(sample_name,fastqfile) %>%
  left_join(SRA %>% distinct(sample_name,Project), by='sample_name') %>%
  filter(!(sample_name %in% pick_acc$sample_name),   # 192 pick creek samples
         !(sample_name %in% og_acc$sample_name),     # this shouldn't drop anything, but to be sure
         !(sample_name %in% dropNoDateData)) %>%     # drop pink samples without collectiond date that weren't used in analyses
  mutate(Project=gsub(" ",'',Project)) %>%
  group_by(Project) %>%
  group_walk(~ { 
    write_tsv(select(.x,fastqfile), str_c("SRA/FastqFilesByProject_",.y$Project, ".txt"),col_names = F)
  })

# in sedna, i created subfolders for each fastq and moved them there
# while IFS= read -r f || [ -n "$f" ]; do [ -n "$f" ] && [ -e "$f" ] && { case "$f" in *.gz) cp -v -- "$f" ../PinkRangewide/ ;; *) \gzip -v -- "$f" && cp -v -- "${f}.gz" ../PinkRangewide/ ;; esac; }; done < ../FastqFilesByProject_PinksRangewide.txt
WGSandGTseqSamples <- bind_rows(og_acc,pick_acc) %>%
  inner_join(SRA, by="sample_name")

# why are there duplicates?????????????????????????????????

# Sequence columns but replace sample_name with accession
WGSandGTseqSamples %>%
  dplyr::select(all_of(c('accession',seqCols))) %>%
  dplyr::select(-sample_name) %>%
  relocate(accession, .before=1) %>%
  rename(biosample_accession=accession) %>%
  write_tsv("./SRA/SRA_sequencing_preparation_PickCreekGTseqAdditions.tsv")

# don't think i need this since the metadata was already imported in the WGS submission
WGSandGTseqSamples %>%
  dplyr::select(all_of(metaCols)) %>%
  write_tsv("./SRA/SRA_metadata_preparation_WGSandGTseqAdditions.tsv")


### BASH COMMANDS ####
# srun -c 5 -t 48:00:00 --pty /bin/bash
# module load tools/lftp/4.9.2
# lftp subftp@ftp-private.ncbi.nlm.nih.gov # enter password given on SRA
# cd /uploads/natashahowe15_gmail.com_ua4ygaIi/SUB15992991_pick_runtiming
# mirror -Rv --continue --parallel=5 --use-pget-n=2 /scratch2/nhowe/runtiming/sock-pick/sednagold/Oncorhynchus_nerka/sockeye0425/ .

# Pink Rangewide
# lftp subftp@ftp-private.ncbi.nlm.nih.gov:/uploads/natashahowe15_gmail.com_LZS2pxfx/PinkRangewideGTseq
# PW: B1ckyadixdaycsuD

# just to note, but i think ive had this before
# chmod: Access failed: 550 CHMOD 644 ./P104_PSITU14_8_65727_H7.fastq.gz: Permission denied
# chmod: Access failed: 550 CHMOD 644 ./P104_PSITU14_9_65727_A8.fastq.gz: Permission denied
# New: 951 files, 0 symlinks
# 5067173056 bytes transferred in 176 seconds (27.51 MiB/s)
# 951 errors detected