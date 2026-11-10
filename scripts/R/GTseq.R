#### GTseq Analysis #### 
library(tidyverse)

#### Plotting params
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                "#0072B2", "#D55E00", "#CC79A7", "#999999")

# for science, consistent font sizing
width2col <- 7.25  # 2 col, in inches
width1col <- 3.6  # 1 col, in inches

largeFS <- 10 # Panel Labels
mediumFS <- 8  # Axis titles, legend titles
smallFS <- 6  # Axis text, annotations, captions

theme_set(
  theme(legend.text=element_text(size=mediumFS),
        legend.title = element_text(size = mediumFS),
        axis.text = element_text(size = smallFS, angle = 0),
        axis.title.x = element_text(size = mediumFS),
        axis.title.y = element_text(size = mediumFS, angle=90, margin=margin(0,0,0,0)),
        plot.title = element_text(size = mediumFS, hjust = 0.5, margin=margin(0,0,2,0)),
        panel.background = element_rect(fill = "white"), 
        panel.spacing = unit(0,"lines"),
        strip.text.y = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )
)

pointSize = 1.5

####

#Define the directories
Path2Flash <- normalizePath("./bin/FLASH2-master/") 
scripts_dir <- "./scripts"
pl_script <- normalizePath(file.path(scripts_dir,"R/functions/AmpliconReadCounter.pl"))

source(file.path(scripts_dir,"R/functions/GTscore.R"))
source(file.path(scripts_dir,"R/functions/ScoreGTSeq_Fncts.R"))

dir.create("./results/GTseq")

#### Pink Salmon ####
#### Auke Creek 2019 Samples ####
dir.create("./results/GTseq/Pink")

# Primer probe file has 1 locus: OgoE1_LG10_72730177
Pink_pp_filepath <-normalizePath(file.path("./data/raw/GTseq/Pink/PinkSalmon_Adaptive1L_v1.0.0.txt"))

Pinkmeta <- read_delim("./data/raw/AukeCreekPink_Recovery2019.txt",
                               delim="\t") %>%
  mutate(Date = as.Date(MarkDate),
         JD = format(Date, "%j"),
         Year = format(Date, "%Y"))

OutDir = "./results/GTseq/Pink"

MiSeqDir = "./data/raw/GTseq/Pink/MiSeqRuns"

ScoreMiSeqRuns(MiSeqDir = MiSeqDir, #MiSeq Sample Sheets should be in a SampleSheets subdirectory
               OutDir = OutDir, 
               primer_probe_filepath = Pink_pp_filepath)

MiSeqRuns<-list.dirs(MiSeqDir,recursive = F) %>%
  {fs::path_file(fs::path_file(.))}

#Read in the genotypes
MixFiles<-list.files(OutDir,recursive = T) %>%
  grep(pattern="polygenResults_singleSNP_rubias.txt",x=.,value=T)%>%
  {grep(pattern=paste(MiSeqRuns,collapse="|"),.,value=T)}

DepthFiles <- list.files(OutDir,recursive = T) %>%
  grep(pattern="AlleleReads_singleSNPs",x=.,value=T)%>%
  {grep(pattern=paste(MiSeqRuns,collapse="|"),.,value=T)}

MixGenos<-lapply(1:length(MixFiles),function(x) read_delim(file.path(OutDir,MixFiles[x]),
                                                           delim="\t",
                                                           col_types = cols(.default = "c"))%>%
                   mutate(collection = lapply(str_split(MixFiles[x],"/"),"[[",1) %>%
                            gsub(pattern="Analyzed_",replacement="",.)) %>%
                   as.data.frame()) %>%
  do.call(rbind,.)

Depths <- lapply(1:length(DepthFiles),function(x) read_delim(file.path(OutDir,DepthFiles[x]),
                                                             delim="\t",
                                                             col_types = cols(.default = "c"))%>%
                   t()%>%
                   .[-1,] %>%
                   {tibble(RdCnts = .,
                           ABLG = names(.) %>% as.numeric())}%>%
                   mutate(TotalReadsLst = lapply(RdCnts, function(x) {str_split(x, ",") %>% 
                                                                         unlist() %>% 
                                                                         as.numeric() %>% 
                                                                         sum()
                                                                     }) %>% unlist())
) %>% do.call(rbind,.)

Pinks2019 <- MixGenos %>%
  filter(!grepl(pattern="Undetermined|QC|NTC|empty",indiv))%>% #remove the 
  mutate(ABLG = gsub(pattern="^X",replacement="",indiv) %>% as.numeric())%>%
  left_join(.,Pinkmeta, by="ABLG") %>%
  select(collection,ABLG,starts_with("Ogo"),Date,JD,Year)%>%
  mutate(Week = lubridate::isoweek(Date)) %>%
  left_join(.,Depths,by="ABLG")

#
Pinks2019 %>%
  filter(!is.na(OgoE1_LG10_72730177_1)) %>%
  filter(!is.na(JD)) %>%
  select(TotalReadsLst) %>% unlist() %>%
  as.numeric() %>%
  mean()

# Try a weekly estimate
Pink_WeeklyFreq <- Pinks2019 %>%
  select(Week, Year,ABLG, starts_with("Ogo")) %>%
  tidyr::pivot_longer(
    cols = -c( Week, Year,ABLG), 
    names_to = "locus", 
    values_to = "allele"
  ) %>%
  mutate(locus = gsub("_1.1", replacement="_1", locus)) %>% 
  ungroup() %>%
   mutate(Week = case_when(
     Week %in% c(32:33) ~ "32-33",
     TRUE ~ as.character(Week)
   )) %>%
  filter(!is.na(allele)) %>%
  count(locus, allele,Year,Week) %>%                       
  group_by(locus,Year,Week) %>%                            
  mutate(frequency = n / sum(n)) %>%            
  ungroup()

# Decrease in pink salmon lrrc9 frequency through week in 2019 at Auke Creek Weir
# not included in manuscipt
p_Pinks <- Pink_WeeklyFreq %>%
  filter(!is.na(Week))%>%
  mutate(locus = gsub(pattern="_1$",replacement="",locus))%>%
  filter(locus == 'OgoE1_LG10_72730177' & allele == 'A') %>%
  mutate(locus_rn = ifelse(locus == 'OgoE1_LG10_72730177','lrrc9',NA))%>%
  ggplot(.,aes(x=as.factor(Week),y=frequency,group=Year,color=as.factor(Year)))+
  geom_point(aes(fill=as.factor(Year),shape=as.factor(Year)),color="black",size=3)+
  geom_smooth(method = "lm", se = T)+
  #scale_x_discrete(breaks = WeekLabs$Week, labels = WeekLabs$lab)+
  ylab('Frequency of Early Allele')+
  xlab("Week")+
  scale_color_manual(values=cb_palette[c(3,4,5,6)])+
  scale_fill_manual(values=cb_palette[c(3,4,5,6)])+
  scale_shape_manual(values=c(21,22,23,24))+
  scale_y_continuous(labels=scales::percent)+
  theme_bw()+
  labs(color = "Year", 
       fill = "Year", 
       shape = "Year")+
  coord_cartesian(ylim=c(0,1))

#Pink_WeeklyFreq %>%
#  filter(!is.na(Week))%>%
#  mutate(locus = gsub(pattern="_1$",replacement="",locus))%>%
#  filter(locus == 'OgoE1_LG10_72730177' & allele == 'A')

# Lets take compute the first and last quarter of the return to auke creek in 2019

AukeCreekRun2019 <- read_delim("./data/raw/Auke_pink_returns_80to22.csv",
                             delim=",") %>%
  mutate(Date_fmt = as.Date(Date,format = "%d-%b-%y")) %>%
  mutate(TotalReturn = `Pink Male Release`+`Pink Female Release`) %>%
  filter(!is.na(TotalReturn)) %>%
  mutate(Year = year(Date_fmt),
         JD = yday(Date_fmt)) %>%
  filter(Year==2019) %>%
  mutate(CumRun=cumsum(TotalReturn),
         CumProp = CumRun/sum(TotalReturn)) %>%
  filter(TotalReturn!=0)

# Early is prior to 23 Aug, Late is post 28 August
AukeCreekRun2019 %>% 
  filter(Date_fmt %in% c("2019-08-23","2019-08-28"))

#Corresponds to JD 235 and 240

PinkSalmon_2019Q1_Q4 <- Pinks2019 %>%
  filter(!is.na(OgoE1_LG10_72730177_1)) %>% #filter out unscored samples
  mutate(Pheno = ifelse(JD<235,'Early',ifelse(JD>240,'Late',NA)))%>%
  filter(!is.na(Pheno)) %>%
  mutate(Geno = paste(OgoE1_LG10_72730177_1, OgoE1_LG10_72730177_1.1,sep=""))%>%
  group_by(Pheno,Geno) %>%
  summarize(nInd = n())

write.csv(file="./results/GTseq/Pink/PinkSalmon_2019Q1_Q4.csv",
          PinkSalmon_2019Q1_Q4,
          row.names = F,
          quote=F)

#SeqDepth
Pinks2019 %>%
  filter(!is.na(OgoE1_LG10_72730177_1)) %>% #filter out unscored samples
  mutate(Pheno = ifelse(JD<235,'Early',ifelse(JD>240,'Late',NA)))%>%
  filter(!is.na(Pheno)) %>%
  group_by(Pheno) %>%
  summarize(RdDepth = mean(TotalReadsLst))


#### Rangewide Collections ####

PinkPopMeta <- readxl::read_xlsx('./data/raw/pink_rangewide_adaptive_popmetadata.xlsx')

NextSeqDir = "./data/raw/GTseq/Pink/NextSeqRuns"

ScoreNextSeqRuns(NextSeqDir = NextSeqDir, #No sample sheets here
               OutDir = OutDir, 
               primer_probe_filepath = Pink_pp_filepath)

NextSeqRuns<-list.dirs(NextSeqDir,recursive = F) %>%
  {fs::path_file(fs::path_file(.))}

MixFiles<-list.files(OutDir,recursive = T) %>%
  grep(pattern="polygenResults_singleSNP_rubias.txt",x=.,value=T)%>%
  {grep(pattern=paste(NextSeqRuns,collapse="|"),.,value=T)}

DepthFiles <- list.files(OutDir,recursive = T) %>%
  grep(pattern="AlleleReads_singleSNPs",x=.,value=T)%>%
  {grep(pattern=paste(NextSeqRuns,collapse="|"),.,value=T)}

MixGenos<-lapply(1:length(MixFiles),function(x) read_delim(file.path(OutDir,MixFiles[x]),
                                                           delim="\t",
                                                           col_types = cols(.default = "c"))%>%
                   mutate(collection = lapply(str_split(MixFiles[x],"/"),"[[",1) %>%
                            gsub(pattern="Analyzed_",replacement="",.)) %>%
                   as.data.frame()) %>%
  do.call(rbind,.)

Depths <- lapply(1:length(DepthFiles),function(x) read_delim(file.path(OutDir,DepthFiles[x]),
                                                             delim="\t",
                                                             col_types = cols(.default = "c"))%>%
                   t()%>%
                   .[-1,] %>%
                   {tibble(RdCnts = .,
                           SILLY_CODE = names(.) %>% as.character())}%>%
                   mutate(SILLY = lapply(str_split(pattern="_",SILLY_CODE),"[[",2) %>% unlist(),
                          Samp = lapply(str_split(pattern="_",SILLY_CODE),"[[",3) %>% unlist())%>%
                   mutate(TotalReadsLst = lapply(RdCnts, function(x) {str_split(x, ",") %>% 
                       unlist() %>% 
                       as.numeric() %>% 
                       sum()
                   }) %>% unlist())
) %>% do.call(rbind,.)

Pinks_RW <- MixGenos %>%
  filter(!grepl(pattern="Undetermined|QC|NTC|empty",indiv))%>% #remove the 
  mutate(SILLY = lapply(str_split(pattern="_",indiv),"[[",2) %>% unlist(),
         Samp = lapply(str_split(pattern="_",indiv),"[[",3) %>% unlist())%>%
  left_join(.,PinkPopMeta,by="SILLY")%>%
  left_join(.,Depths,by=c("SILLY",'Samp'))

#Depths
Pinks_RW %>%
  mutate(JD = yday(`samp date`))%>%
  filter(!is.na(OgoE1_LG10_72730177_1)) %>%
  filter(!is.na(JD)) %>%
  select(TotalReadsLst) %>% unlist() %>%
  as.numeric() %>%
  mean()

# Try a weekly estimate
Pink_RWpopFreq <- Pinks_RW %>%
  select(`Pop name`, `samp date`, lineage,Region, starts_with("Ogo")) %>%
  tidyr::pivot_longer(
    cols = -c( `Pop name`, `samp date`, lineage,Region), 
    names_to = "locus", 
    values_to = "allele"
  ) %>%
  mutate(locus = gsub("_1.1", replacement="_1", locus)) %>% 
  ungroup() %>%
  filter(!is.na(allele)) %>%
  count(locus, allele,`Pop name`, `samp date`, Region,lineage) %>%                       
  group_by(locus,`Pop name`, `samp date`, Region,lineage) %>%                            
  mutate(frequency = n / sum(n)) %>%            
  ungroup() %>%
  mutate(JD = yday(`samp date`)) 

pointSize<-2

p_PinksRW <- Pink_RWpopFreq %>%
  filter(locus == 'OgoE1_LG10_72730177_1' & allele == 'A') %>% #subset on early allele
  filter(!is.na(JD))%>%
  mutate(Region = factor(Region,levels=c('Russia','Western AK','Kodiak','Prince William Sound','Southeast AK')))%>%
  ggplot(.,aes(x=JD,y=frequency))+
  geom_point(aes(shape=Region,color=lineage),size=pointSize)+
  geom_smooth(method='lm', formula= y~x)+
  ylab(expression(paste("Frequency of Early ",italic("lrrc9"), " Allele",sep="")))+
  xlab("Julian Day")+
  scale_y_continuous(labels=scales::percent)+
  scale_shape_manual(values = c(25,15,1,17,4))+
  geom_point(data=Pink_RWpopFreq %>%
               filter(locus == 'OgoE1_LG10_72730177_1' & allele == 'A') %>%
               filter(grepl(pattern='Duck',`Pop name`)),aes(x=JD,y=frequency),pch=2,size=pointSize,color='black') +
  coord_cartesian(ylim=c(0,1))+
  scale_color_manual(values=cb_palette[c(2,6)])+
  labs(color = "Lineage", 
       shape = "Region")
  
# for import into MS
png(str_c("./figures/FigS3_Pinklrrc9Rangewide",format(Sys.Date(),"%Y%m%d"),".png"),
    width= width2col, height= 0.75*width2col, res = 500, units = "in")
p_PinksRW
dev.off()

# svg for science 
svglite::svglite(str_c("./figures/FigS3_Pinklrrc9Rangewide",format(Sys.Date(),"%Y%m%d"),".svg"),
        width= width2col, height= 0.75*width2col)
print(p_PinksRW)
dev.off()

#Rangewide lm 

fit_RWv1 <- lm(frequency ~ JD+lineage, data = Pink_RWpopFreq %>%
                                          filter(locus == 'OgoE1_LG10_72730177_1' & allele == 'A') %>%
                                          filter(!is.na(JD)))
summary(fit_RWv1)
PinkRW_model <- summary(fit_RWv1)$coefficients

fit_RWv1.1 <- lm(frequency ~ JD+Region+lineage, data = Pink_RWpopFreq %>%
                 filter(locus == 'OgoE1_LG10_72730177_1' & allele == 'A') %>%
                 filter(!is.na(JD)))
summary(fit_RW)
PinkRW_model <- summary(fit_RWv1.1)$coefficients

# Test for JD:Region interaction to see if the allele freq shift changes depending on region
fit_RWv1.2 <- lm(frequency ~ JD*Region+lineage, data = Pink_RWpopFreq %>%
               filter(locus == 'OgoE1_LG10_72730177_1' & allele == 'A') %>%
               filter(!is.na(JD)))

summary(fit_RWv1.2)

fit_RWv1.3 <- lm(frequency ~ JD*Region, data = Pink_RWpopFreq %>%
                 filter(locus == 'OgoE1_LG10_72730177_1' & allele == 'A') %>%
                 filter(!is.na(JD)))

summary(fit_RWv1.3)

#test for lineage by region effect
fit_RWv2 <- lm(frequency ~ JD+Region*lineage, data = Pink_RWpopFreq %>%
                 mutate(total = round(n/frequency,0))%>%
                 filter(locus == 'OgoE1_LG10_72730177_1' & allele == 'A') %>%
                 filter(!is.na(JD)))

summary(fit_RWv2)


  
#### Chum Salmon ####
# We will analyze samples from the 2024 outmigration of chum salmon 
# from the Yukon River.

dir.create("./results/GTseq/Chum")

# Primer probe file has 2 loci: OkeV2_LG29_25455833, and OkeV2_LG35_28165076

Chum_pp_filepath <-normalizePath(file.path("./data/raw/GTseq/Chum/ChumSalmon_Adaptive2L_v1.0.0.txt"))

JYC24meta <- readxl::read_xlsx("./data/raw/2024 Yukon River Juvenile Chum Metadata 111324.xlsx",
                               sheet = "Sheet1") %>%
  select(sin,station_name,visit_date,length,mass_g) %>%
  rename('SampleID' = 1,
         'Station' = 2,
         'Date' = 3,
         'Length' = 4,
         'Weight' = 5)

OutDir = "./results/GTseq/Chum"

ScoreMiSeqRuns(MiSeqDir = "./data/raw/GTseq/Chum/MiSeqRuns", #MiSeq Sample Sheets should be in a SampleSheets subdirectory
               OutDir = OutDir, 
               primer_probe_filepath = Chum_pp_filepath)

MiSeqDir = "./data/raw/GTseq/Chum/MiSeqRuns"

MiSeqRuns<-list.dirs(MiSeqDir,recursive = F) %>%
  {fs::path_file(fs::path_file(.))}

#Read in the genotypes
MixFiles<-list.files(OutDir,recursive = T) %>%
  grep(pattern="polygenResults_singleSNP_rubias.txt",x=.,value=T)

MixGenos<-lapply(1:length(MixFiles),function(x) read_delim(file.path(OutDir,MixFiles[x]),
                                                           delim="\t",
                                                           col_types = cols(.default = "c"))%>%
                   mutate(collection = lapply(str_split(MixFiles[x],"/"),"[[",1) %>%
                            gsub(pattern="Analyzed_",replacement="",.)) %>%
                   as.data.frame()) %>%
  do.call(rbind,.)

DepthFiles <- list.files(OutDir,recursive = T) %>%
  grep(pattern="AlleleReads_singleSNPs",x=.,value=T)%>%
  {grep(pattern=paste(MiSeqRuns,collapse="|"),.,value=T)}

Depths <- lapply(1:length(DepthFiles),function(x) read_delim(file.path(OutDir, DepthFiles[x]), 
                                                             delim = "\t", 
                                                             col_types = cols(.default = "c")) %>%
                   t() %>%
                   .[-1, ] %>%
                   { tibble(
                     RdCnts_esrb = .[, 1],
                     RdCnts_lrrc9 = .[, 2],
                     ABLG = rownames(.))} %>%
                   mutate(TotalReads_esrb = lapply(RdCnts_esrb, function(x) {str_split(x, ",") %>% 
                       unlist() %>% 
                       as.numeric() %>% 
                       sum()
                   }) %>% unlist(),
                   TotalReads_lrrc9 = lapply(RdCnts_lrrc9, function(x) {str_split(x, ",") %>% 
                       unlist() %>% 
                       as.numeric() %>% 
                       sum()}) %>% unlist(),
                   MiSeqRun = MiSeqRuns[x])) %>% do.call(rbind,.)

#Filter Mix genos to remove any sample run 2x in different miseq runs,
#Keep the run with the highest read depth

MixGenosDeDup <- MixGenos %>%
  filter(!grepl(pattern="Undetermined|QC|NTC|empty",indiv))%>% #remove the 
  #mutate(indiv = gsub(pattern="[a-b]$",replacement="",indiv))%>%
  mutate(SampleID = gsub(pattern="\\.",replacement="_",indiv),
         MiSeqRun = collection)%>%
  left_join(Depths %>% mutate(SampleID=gsub(pattern="-",replacement="_",ABLG)),by=c("SampleID","MiSeqRun")) %>%
  mutate(DDindiv = gsub(pattern="[a-b]$",replacement="",indiv))%>%
  group_by(DDindiv) %>%
  mutate(TotalReads = TotalReads_esrb+TotalReads_lrrc9)%>%
  slice_max(.,TotalReads) %>%
  ungroup()

DrpdInds <- MixGenos %>%
  filter(!grepl(pattern="Undetermined|QC|NTC|empty",indiv))%>% #remove the 
  #mutate(indiv = gsub(pattern="[a-b]$",replacement="",indiv))%>%
  mutate(SampleID = gsub(pattern="\\.",replacement="_",indiv),
         MiSeqRun = collection)%>%
  left_join(Depths %>% mutate(SampleID=gsub(pattern="-",replacement="_",ABLG)),by=c("SampleID","MiSeqRun")) %>%
  mutate(DDindiv = gsub(pattern="[a-b]$",replacement="",indiv))%>%
  group_by(DDindiv) %>%
  filter(n() > 1) %>%
  mutate(TotalReads = TotalReads_esrb+TotalReads_lrrc9)%>%
  slice_min(.,TotalReads,with_ties = F) %>%
  ungroup()

# Need to make a table with full path to each file and the MiSeq run 
# so we know which to drop
ChumFiles <- list.files("./data/raw/GTseq/Chum/",recursive = T,full.names = T) %>%
  grep(pattern=".fastq",.,value=T)%>%
  {.[!grepl(pattern="Undetermined|NTC|QC|empty",.)]}

DeDupTib <- tibble(FileName = ChumFiles) %>%
  mutate(MiSeqRun = lapply(str_split(FileName,pattern="/"),"[[",8) %>% unlist(),
         Sample = lapply(str_split(FileName,pattern="/"),"[[",9) %>% unlist(),
         Fastq = lapply(str_split(FileName,pattern="/"),"[[",10) %>% unlist(),
         indiv = lapply(str_split(Sample,pattern="_S"),"[[",1) %>% unlist(),
         collection = MiSeqRun)


DroppedTib<-DeDupTib %>%
  left_join(.,DrpdInds%>%
              mutate(indiv=gsub(pattern="\\.",replacement="-",indiv)),
            by=c('indiv','collection','MiSeqRun')) %>%
  filter(!is.na(SampleID))

write_delim(DroppedTib,
            "./results/GTseq/DroppedJuvenileChum.txt",
            delim="\t")

#

MixGenosDeDup %>%
  filter(!is.na(OkeV2_LG29_25455833_1)) %>%
  filter(!is.na(OkeV2_LG35_28165076_1))%>%
  filter(!grepl(pattern="^X",DDindiv))%>%
  select(DDindiv,TotalReads_esrb,TotalReads_lrrc9) %>%
  pivot_longer(cols = starts_with('TotalReads'), 
               names_to = "variable", 
               values_to = "value") %>%
  ungroup() %>%
  group_by(variable) %>%
  summarize(MeanDepth = mean(value))

#How many samples genotyped for at least one locus 
MixGenosDeDup %>%
       filter(!grepl(pattern="Undetermined|QC|NTC|empty",indiv))%>% #remove the 
       mutate(indiv = gsub(pattern="[a-b]$",replacement="",indiv))%>%
       mutate(SampleID = gsub(pattern="\\.",replacement="_",indiv))%>%
       left_join(.,JYC24meta %>% mutate(SampleID=gsub(pattern="-",replacement="_",SampleID)), by="SampleID") %>%
       filter(!is.na(Date))%>%filter(!grepl(pattern="^X",DDindiv))%>%filter(!(is.na(OkeV2_LG29_25455833_1)&is.na(OkeV2_LG35_28165076_1))) %>%
       mutate(Year = substr(Date, start = 1, stop = 4),
                          JD = format(Date, "%j")) %>%
       mutate(Week = lubridate::isoweek(Date))%>%filter(!is.na(JD)) %>%
  nrow()

#What was the sequencing depth of those samples?
MixGenosDeDup %>%
  filter(!grepl(pattern="Undetermined|QC|NTC|empty",indiv))%>% #remove the 
  mutate(indiv = gsub(pattern="[a-b]$",replacement="",indiv))%>%
  mutate(SampleID = gsub(pattern="\\.",replacement="_",indiv))%>%
  left_join(.,JYC24meta %>% mutate(SampleID=gsub(pattern="-",replacement="_",SampleID)), by="SampleID") %>%
  filter(!is.na(Date))%>%filter(!grepl(pattern="^X",DDindiv))%>%filter(!(is.na(OkeV2_LG29_25455833_1)&is.na(OkeV2_LG35_28165076_1))) %>%
  mutate(Year = substr(Date, start = 1, stop = 4),
         JD = format(Date, "%j")) %>%
  mutate(Week = lubridate::isoweek(Date))%>%filter(!is.na(JD))%>%
  select(TotalReads) %>%
  unlist()%>%
  mean()

JuvChum <- MixGenosDeDup %>%
  filter(!grepl(pattern="Undetermined|QC|NTC|empty",indiv))%>% #remove the 
  mutate(indiv = gsub(pattern="[a-b]$",replacement="",indiv))%>%
  mutate(SampleID = gsub(pattern="\\.",replacement="_",indiv))%>%
  left_join(.,JYC24meta %>% mutate(SampleID=gsub(pattern="-",replacement="_",SampleID)), by="SampleID") %>%
  filter(!is.na(Date))%>%
  mutate(Year = substr(Date, start = 1, stop = 4),
         JD = format(Date, "%j")) %>%
  mutate(Week = lubridate::isoweek(Date)) %>%
  ungroup() %>%
  select(Year, Week, indiv, starts_with("Oke")) %>%
  tidyr::pivot_longer(
    cols = -c(Year, Week, indiv), 
    names_to = "locus", 
    values_to = "allele"
  ) %>%
  mutate(locus = gsub("_1.1", replacement="_1", locus)) %>% 
  ungroup() %>%
  mutate(Week = case_when(
    Week %in% c(30, 31) ~ "30+",
    TRUE ~ as.character(Week)
  )) 

WeeklyFreqs <- JuvChum %>%
  filter(!is.na(allele)) %>%
  count(locus, allele,Week) %>%                       
  group_by(locus,Week) %>%                            
  mutate(frequency = n / sum(n)) %>%            
  ungroup()

WeekLabs <- WeeklyFreqs %>% 
  select(locus,Week,n) %>%
  group_by(locus,Week) %>%
  summarize(nInds = sum(n)/2)%>%
  tidyr::pivot_wider(names_from = locus, values_from = nInds)%>%
  mutate(lab = paste(Week,"\n","(n=",OkeV2_LG29_25455833_1,"/",OkeV2_LG35_28165076_1,")",sep=""))

p_Chum <- WeeklyFreqs %>%
  mutate(locus = gsub(pattern="_1$",replacement="",locus))%>%
  filter(locus == 'OkeV2_LG29_25455833' & allele == 'T' | locus == "OkeV2_LG35_28165076" & allele == 'C') %>%
  mutate(locus_rn = ifelse(locus == 'OkeV2_LG29_25455833','esr2b',ifelse(locus == "OkeV2_LG35_28165076",'lrrc9',NA)))%>%
  ggplot(.,aes(x=as.factor(Week),y=frequency,group=locus_rn,color=locus_rn))+
  geom_point(aes(shape=locus_rn,fill=locus_rn),color="black",size=pointSize)+
  geom_smooth(method = "lm", se = TRUE)+
  scale_x_discrete(breaks = WeekLabs$Week, labels = WeekLabs$lab)+
  ylab('Frequency of Early Allele')+
  xlab("Week")+
  scale_color_manual(values=cb_palette[c(5,6)])+
  scale_fill_manual(values=cb_palette[c(5,6)])+
  scale_shape_manual(values=c(21,22))+
  scale_y_continuous(labels=scales::percent)+
  labs(color = "Locus", 
       fill = "Locus", 
       shape = "Locus")
  
#Model Fits  
Locus_Fits <-  WeeklyFreqs %>%
  mutate(locus = gsub(pattern="_1$",replacement="",locus))%>%
  filter(locus == 'OkeV2_LG29_25455833' & allele == 'T' | locus == "OkeV2_LG35_28165076" & allele == 'C') %>%
  mutate(NumWeek = as.numeric(as.factor(Week))) %>%
  group_by(locus) %>%
  do(fit = lm(frequency ~ NumWeek, data = .))

summary(Locus_Fits$fit[[1]]) #OkeV2_LG29_25455833
summary(Locus_Fits$fit[[2]]) #OkeV2_LG35_28165076

 
#### Sockeye Salmon ####
dir.create("./results/GTseq/Sockeye")

# Primer probe file has 2 loci: One2_LG13_27674674, and One_RAD27165

Sockeye_pp_filepath <-normalizePath(file.path("./data/raw/GTseq/Sockeye/SockeyeSalmon_Adaptive_2L_v1.0.0Clean.txt"))

# There are three metadata files for these samples, and the database columns are different between each,
# So we have to deal with them individually, merge them, then put the nicely formatted version in 
# the data/raw directory

Meta1 <- readxl::read_xlsx("./data/raw/GTseq/Sockeye/ABLG_exports_Pick_WhitefishSockeye.xlsx",
                               sheet = "Sheet1",col_types = "text") 

Meta2 <- readxl::read_xlsx("./data/raw/GTseq/Sockeye/ABLG_exports_PickCreekSockeye.xlsx",
                           sheet = "Sheet1",col_types = "text") 

Meta3 <- readxl::read_xlsx("./data/raw/GTseq/Sockeye/20260728_SWHIT11_P985.xlsx",
                           sheet = "ABLG-DB_working-copy",col_types = "text") 

PickWhitefish_meta <- bind_rows(Meta1,Meta2,Meta3) %>%
  select(ABLG,AlternateID,
         CollectionYear,CollectionMonth,CollectionDay,
         FreshwaterWatershed,Location1,Location2,   
         Sex, SexMorF) %>%
  mutate(SEX = ifelse(!is.na(Sex),Sex,SexMorF))%>%
  select(-c(Sex,SexMorF))%>%
  mutate(Location = ifelse((is.na(Location1)&grepl(pattern="WHIT11",AlternateID)),"Whitefish Creek",Location1))%>%
  select(-c(Location1,Location2))%>%
  filter(!is.na(ABLG)) %>%
  mutate(CollectionDate = str_c(CollectionYear,'-',CollectionMonth,'-',CollectionDay),
         LocaleDate = str_c(sub(' Creek','',Location),'.',CollectionDate))

rm(Meta1,Meta2,Meta3)

OutDir = "./results/GTseq/Sockeye"

MiSeqDir = "./data/raw/GTseq/Sockeye/MiSeqRuns"

MiSeqRuns<-list.dirs(MiSeqDir,recursive = F) %>%
  {fs::path_file(fs::path_file(.))}

ScoreMiSeqRuns(MiSeqDir = "./data/raw/GTseq/Sockeye/MiSeqRuns", #MiSeq Sample Sheets should be in a SampleSheets subdirectory
               OutDir = OutDir, 
               primer_probe_filepath = Sockeye_pp_filepath)

#Read in the genotypes
MixFiles<-list.files(OutDir,recursive = T) %>%
  grep(pattern="polygenResults_singleSNP_rubias.txt",x=.,value=T)

MixGenos<-lapply(1:length(MixFiles),function(x) read_delim(file.path(OutDir,MixFiles[x]),
                                                           delim="\t",
                                                           col_types = cols(.default = "c"))%>%
                   mutate(collection = lapply(str_split(MixFiles[x],"/"),"[[",1) %>%
                            gsub(pattern="Analyzed_",replacement="",.)) %>%
                   as.data.frame()) %>%
  do.call(rbind,.)

DepthFiles <- list.files(OutDir,recursive = T) %>%
  grep(pattern="AlleleReads_singleSNPs",x=.,value=T)%>%
  {grep(pattern=paste(MiSeqRuns,collapse="|"),.,value=T)}

Depths <- lapply(1:length(DepthFiles),function(x) read_delim(file.path(OutDir, DepthFiles[x]), 
                                                             delim = "\t", 
                                                             col_types = cols(.default = "c")) %>%
                   t() %>%
                   .[-1, ] %>%
                   { tibble(
                     RdCnts_esr1 = .[, 1],
                     RdCnts_lrrc9 = .[, 2],
                     ABLG = rownames(.))} %>%
                   mutate(TotalReads_esr = lapply(RdCnts_esr1, function(x) {str_split(x, ",") %>% 
                       unlist() %>% 
                       as.numeric() %>% 
                       sum()
                   }) %>% unlist(),
                   TotalReads_lrrc9 = lapply(RdCnts_lrrc9, function(x) {str_split(x, ",") %>% 
                       unlist() %>% 
                       as.numeric() %>% 
                       sum()}) %>% unlist(),
                   MiSeqRun = MiSeqRuns[x])) %>% do.call(rbind,.)

Depths %>%
  as.data.frame() %>%  # Strips hidden, corrupt object attributes
  filter(!grepl("Undetermined|QC|NTC|empty", ABLG)) %>%
  select(ABLG, starts_with('TotalReads')) %>%
  pivot_longer(cols = starts_with('TotalReads'), 
               names_to = "variable", 
               values_to = "value") %>%
  group_by(variable) %>%
  summarize(MeanDepth = mean(value, na.rm = TRUE))

SockeyeGenos <- MixGenos %>%
  filter(!grepl(pattern="Undetermined|QC|NTC|empty",indiv))%>% #remove the 
  mutate(indiv = gsub(pattern="^X",replacement="",indiv))%>%
  left_join(.,PickWhitefish_meta %>% 
              mutate(indiv=as.character(ABLG)) , by="indiv") %>%
  filter(!is.na(CollectionDate))%>%
  mutate(Date = as.Date(CollectionDate), #convert to date obj to pull the Julian
         JD = format(Date, "%j")) %>%
  ungroup() 

SockeyeGenos %>%
  filter(!(is.na(One2_LG13_27674674_16)&is.na(One_RAD27165_52)))%>% #remove any ind. that aren't genotyped for both
  left_join(.,Depths %>% mutate(TtlRds = TotalReads_esr+TotalReads_lrrc9),by="ABLG") %>%
  group_by(LocaleDate) %>% 
  summarize(nSamps = n(),
    AvgDepth = mean(TtlRds))

SockeyeGenos %>%
  group_by(Location) %>%
  dplyr::summarize(nSamps = n())

DailyFreq <- SockeyeGenos %>%
  select(indiv,Location,JD,starts_with("One")) %>%
  tidyr::pivot_longer(
    cols = -c(Location,JD,indiv), 
    names_to = "locus", 
    values_to = "allele"
  ) %>%
  mutate(locus = gsub(".1$", replacement="", locus)) %>% 
  ungroup() %>%
  filter(!is.na(allele)) %>%
  count(locus, allele,JD,Location) %>%                       
  group_by(locus,JD,Location) %>%                            
  mutate(frequency = n / sum(n)) %>%            
  ungroup()

#filter by early allele 
p_Sockeye <- DailyFreq %>%
  mutate(locus_rn = ifelse(locus == "One2_LG13_27674674_16",'esr1',ifelse(locus == "One_RAD27165_52",'lrrc9',NA)))%>%
  mutate(PopLoc = paste(Location, locus_rn, sep=", ")) %>%
  filter(locus_rn=="lrrc9" & allele=="T" | locus_rn== "esr1" & allele=="G") %>% #filter(locus=="One_RAD27165_52" & allele=="T" |locus== "One2_LG13_27674674_16" & allele=="G") %>%
  ggplot(aes(x = JD, y = frequency, group=PopLoc, color=PopLoc)) +
  geom_point(aes(shape=PopLoc, fill=PopLoc), color='black',size=pointSize) +
  geom_smooth(method="lm", se=F) +
  scale_shape_manual(values=c(22, 22, 21, 21)) + 
  ylab('Frequency of Early Allele')+
  xlab('Julian Date')+
  scale_color_manual(values = c(cb_palette[c(3,8,6,7)]))+
  scale_fill_manual(values = c(cb_palette[c(3,8,6,7)]))+
  scale_y_continuous(labels=scales::percent)+
  labs(color = "Population, Locus", 
       fill = "Population, Locus", 
       shape = "Population, Locus")

#Model fits weren't explored because of very low sample sizes by location and locus (n=3)



### Plot for S5 ### 

theme_set(
    theme(
      legend.text = element_text(size = mediumFS),
      legend.title = element_text(size = mediumFS),
      axis.text = element_text(size = smallFS),
      axis.text.x = element_text(angle = 45, hjust = 0.3), 
      axis.title.x = element_text(size = mediumFS),
      axis.title.y = element_text(size = mediumFS, angle = 90, margin = margin(0,0,0,0)),
      plot.title = element_text(size = mediumFS, hjust = 0.5, margin = margin(0,0,2,0)),
      panel.background = element_rect(fill = "white"), 
      panel.spacing = unit(0,"lines"),
      strip.text.y = element_text(angle = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
)
)

pointSize = 2



pChum_Sockeye <- cowplot::plot_grid(p_Chum,
                                    p_Sockeye,
                                    labels = c("A","B"),
                                    nrow=2,
                                    align = 'vh')
# for import into MS
png(str_c("./figures/FigS5_ChumSockeyeGTseq",format(Sys.Date(),"%Y%m%d"),".png"),
    width= width2col, height= width2col, res = 500, units = "in")
pChum_Sockeye
dev.off()

# svg for science 
svglite::svglite(str_c("./figures/FigS5_ChumSockeyeGTseq",format(Sys.Date(),"%Y%m%d"),".svg"),
                 width= width2col, height= width2col)
print(pChum_Sockeye)
dev.off()


