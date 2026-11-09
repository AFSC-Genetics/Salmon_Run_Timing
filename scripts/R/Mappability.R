# Mappability 
library(tidyverse)

Species_v <- c("Chum","Pink","Sockeye","Coho")

Map_Summary <- tibble(Species = rep(Species_v,times=c(6,2,4,4)),
                      Chr = rep(c(8,29,35,10,12,13,7,21),each=2),
                      Bounds = c('23100059 - 23997201',
                                 '23525841 - 23552881',
                                 '24700042 - 263998921',
                                 '25451985 - 25483494',
                                 '27861977 - 28299479',
                                     '28128954 - 28169980',
                                     '72650536 - 72875058',
                                     '72710811 - 72735126',
                                     '41162738 - 41448689',
                                     '41184975 - 41207762',
                                     '7206200  - 7645274',
                                     '7621613 - 7645274',
                                     '29455765 - 30509378',
                                     '30264668 - 30291493',
                                     '24716700 - 24902281',
                                     '24719119 - 24743990'),
                      BoundDes = c('esr1 region',
                                   'esr1 gene',
                                   'esr2b region',
                                   'esr2b gene',
                                   'lrrc9 region',
                                   'lrrc9 gene',
                                   'lrrc9 region',
                                   'lrrc9 gene',
                                   'lrrc9 region',
                                   'lrrc9 gene',
                                   'esr1 region',
                                   'esr1 gene',
                                   'esr2b region',
                                   'esr2b gene',
                                   'esr1 region',
                                   'esr1 gene'))

SumList<-list()
for (i in 1:length(Species_v)){

Sp_ChrNames <- read_delim(paste("./results/Mappability/",Species_v[i],"/",Species_v[i],"_ChrNames.txt",sep=""),
                          delim="\t",
                          col_names = F)%>%
  unlist() %>%
  as.vector()

SumTemp <- Map_Summary %>%
  filter(Species == Species_v[i])%>%
  mutate(ChrName = Sp_ChrNames[Chr])%>%
  mutate(BndStrt = lapply(str_split(Bounds," - "),"[[",1) %>% unlist(),
         BndEnd = lapply(str_split(Bounds," - "),"[[",2)%>% unlist())
         
  
Sp_Maply_RsOI <- read_delim(paste("./results/Mappability/",Species_v[i],"/GenMap_Res/",Species_v[i],"_ROI_Mappability.bedgraph",sep=""),
                             delim="\t",
                             col_names = F) %>%
  rename('Chr' = 1,
         'Strt' = 2,
         'End' = 3,
         'Maply' = 4) 

Sp_Maply_RsOI_Expanded <- Sp_Maply_RsOI %>%
  mutate(Pos = map2(Strt, End, ~ .x:.y)) %>% #make list of positions
  unnest(Pos) %>%
  select(Chr, Pos, Maply)

SumList[[i]] <- lapply(1:nrow(SumTemp),function(x)
  Sp_Maply_RsOI_Expanded %>%
    filter(Chr == SumTemp$ChrName[x],
           Pos>= SumTemp$BndStrt[x] & Pos<= SumTemp$BndEnd[x])%>%
    summarize(Mean = mean(Maply),
              Median = median(Maply),
              SD = sd(Maply),
              Min = min(Maply),
              Max = max(Maply),
              FracUniq = mean(Maply==1))
)%>%
  do.call(rbind,.)%>%
  cbind(SumTemp,.)
}

write.csv(x=do.call(rbind,SumList),
          file="./results/mappability/Mappability_AllSp.csv",
          row.names=F,
          quote=F)
