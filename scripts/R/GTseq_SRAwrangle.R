# GTseq SRA submission Wrangling
# Must have run the GTseq script for environmental variables to populate

dir.create("./data/raw/GTseq/SRA") #make a directory to copy all the files

Files2Cp <- list.files("./data/raw/GTseq",recursive = T,full.names = T) %>%
  grep(pattern=".fastq",.,value=T)%>%
  {.[!grepl(pattern="Undetermined|NTC|QC|empty",.)]}
  
#2019 Auke Pinks
Pink2019_meta <- Pinks2019 %>%
  select(ABLG) %>%
  left_join(.,Pinkmeta,by="ABLG") %>%
  mutate(sample_name= as.character(ABLG),
         organism = 'Oncorhynchus gorbuscha',
         dev_stage = "adult",
         ecotype=NA,
         geo_loc_name = "USA:Alaska, Juneau, Auke Creek",
         sex= 'Female',
         tissue = 'axillary process', 
         lineage = 'Odd',
         collection_date = format(Date,"%Y-%m-%d"),
         lat_lon = "58.381 N 134.641 W") %>%
  select(sample_name,organism,ecotype,dev_stage,collection_date,geo_loc_name,sex,tissue,AlternateID,lat_lon,lineage)
  
#Rangewide pinks
PinksRW_meta <- Pinks_RW %>%
  filter(!is.na(`samp date`)) %>%
  mutate(sample_name = paste(SILLY, Samp,sep="_"),
         organism = 'Oncorhynchus gorbuscha',
         dev_stage = "adult",
         ecotype=NA,
         sex=NA,
         geo_loc_name = ifelse(Region=="Russia",paste(Region, ": ", `Pop name`,sep=""),paste("USA: Alaska, ",`Pop name`,sep="")),
         tissue = 'fin',
         collection_date = format(`samp date`,"%Y-%m-%d"),
         lat_lon = ifelse(Region=="Russia",paste(lat," N ",long," E",sep=""),paste(lat," N ",-1*long," W",sep="")),
         AlternateID = indiv) %>%
  select(sample_name,organism,ecotype,dev_stage,collection_date,geo_loc_name,sex,tissue,AlternateID,lat_lon,lineage)


#Chum Yukon River  
ChumYukon_meta<-MixGenosDeDup %>%
  select(ABLG,SampleID) %>%
  left_join(.,JYC24meta %>% mutate(SampleID=gsub(pattern="-",replacement="_",SampleID)), by="SampleID") %>%
  filter(!grepl(pattern="Undetermined|QC|NTC|empty|^X",SampleID))%>%
  mutate(sample_name = as.character(ABLG),
         organism = 'Oncorhynchus keta',
         dev_stage = "juvenile",
         ecotype=NA,
         sex=NA,
         geo_loc_name = "USA: Alaska, Yukon River",
         tissue = 'fin',
         collection_date = format(Date,"%Y-%m-%d"),
         lat_lon = "62.595 N	163.967 W",
         AlternateID = ABLG,
         lineage=NA) %>%
  select(sample_name,organism,ecotype,dev_stage,collection_date,geo_loc_name,sex,tissue,AlternateID,lat_lon,lineage)

#Sockeye Bristol Bay
SockeyeBB_meta <- SockeyeGenos %>%
  filter(!grepl(pattern="Undetermined|QC|NTC|empty|^X",indiv))%>%
  mutate(sample_name = as.character(ABLG),
         organism = 'Oncorhynchus nerka',
         dev_stage = "adult",
         ecotype='stream',
         sex=SEX,
         geo_loc_name = paste("USA: Alaska, ",Location,sep=""),
         tissue = 'fin',
         collection_date = format(Date,"%Y-%m-%d"),
         lat_lon = ifelse(grepl(pattern="Pick",geo_loc_name),"59.549 N	159.07 W",ifelse(grepl(pattern="Whitefish",geo_loc_name),"59.267 N	158.686 W",NA)),
         AlternateID = paste(Location,sample_name,sep=""),
         lineage=NA) %>%
  select(sample_name,organism,ecotype,dev_stage,collection_date,geo_loc_name,sex,tissue,AlternateID,lat_lon,lineage)%>%
  mutate(sex = gsub(pattern="<NA>",NA,sex))


GtseqMeta <- bind_rows(Pink2019_meta %>%mutate(AlternateID = as.character(AlternateID)) %>% 
                         mutate(Project = "Auke Pinks 2019"),
                   PinksRW_meta%>%mutate(AlternateID = as.character(AlternateID))%>% 
                     mutate(Project = "Pinks Rangewide"),
                   ChumYukon_meta%>%mutate(AlternateID = as.character(AlternateID))%>% 
                     mutate(Project = "Chum Yukon River"),
                   SockeyeBB_meta%>%mutate(AlternateID = as.character(AlternateID))%>% 
                     mutate(Project = "Sockeye Bristol Bay"))

#GTseqMeta has only the samples we want

FileNameWrangle <- tibble('File_FullPath' = Files2Cp)%>%
  mutate(FileName_shrt = basename(File_FullPath))%>%
  mutate(Destn = paste("./data/raw/GTseq/SRA/",FileName_shrt,sep=""),
         Cmd = paste('cp', File_FullPath, Destn,sep=" "))

# FileNameWrangle has some samples that didn't genotype 
# We want to subset down the files to upload

FileNameWrangle <- FileNameWrangle %>%
  mutate(SeqID = ifelse(grepl(pattern = "Chum",File_FullPath), lapply(str_split(FileName_shrt,"_"),"[[",1) %>% unlist(), ifelse(
    grepl(pattern = "Pink/MiSeq",File_FullPath), lapply(str_split(FileName_shrt,"_"),"[[",1) %>% unlist(), ifelse(
    grepl(pattern = "Pink/NovaSeq",File_FullPath), word(FileName_shrt, start = 2, end = 3, sep = "_"), ifelse(
  grepl(pattern = "Sockeye/MiSeq",File_FullPath), lapply(str_split(FileName_shrt,"_"),"[[",1) %>% unlist(), NA)))))
  
    
FN2Rm <- FileNameWrangle$SeqID[!FileNameWrangle$SeqID%in%GtseqMeta$sample_name] 

FileNameWrangle <- FileNameWrangle %>%
  filter(!(SeqID %in% FN2Rm))

for (i in 1:nrow(FileNameWrangle)){
  SysCmd <- FileNameWrangle$Cmd[i] 
  system(SysCmd,wait=T)
}

#Check for missing data - esp. dates

MissingDates <- GtseqMeta %>%
  filter(is.na(collection_date))

# There are some pink salmon that do not have a sampling date, but 
# The chum salmon are just not binding cause they have a letter attached to their identifier

GtseqMeta<-GtseqMeta %>%
  mutate(ChumID = ifelse(grepl(pattern="Chum Yukon River",Project),gsub(pattern="a$",replacement="",sample_name),NA))%>%
  left_join(.,JYC24meta %>% mutate(ChumID = SampleID), by="ChumID") %>%
  mutate(collection_date = ifelse((is.na(collection_date) & grepl(pattern="Chum Yukon River",Project)), format(Date,"%Y-%m-%d"),collection_date)) %>%
  select(1:12)

GtseqMeta %>%
  filter(grepl(pattern="a$",sample_name))


write_delim(GtseqMeta,
            file="./data/raw/GTseq/GTseq_meta.txt",
            delim="\t")
