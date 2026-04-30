
#### Identify Genes #### 

#Define which regions we want to look at
#The boundaries can be chosen by looking at the peaks and manually defining peaks,
#or from local score, but those might not encompass the entire peak we may be interested in.

Genes2Check = data.frame(Species = c('Pink','Chum','Chum','Sockeye','Sockeye','Coho','Coho'),
                         LG = c('NC_060182.1','NC_068449.1','NC_068455.1','NC_042546.1', 'NC_042547.1','NC_034180.2', 'NC_034194.2'),
                         Gene= c('lrrc9','esr2b','lrrc9','lrrc9','esr1','esr2b','esr1'),
                         GeneStart = c(72650536,25440358,28033205,41162738,7206200,29455765,24716700),
                         GeneEnd = c(72875058,25545548,28206418,41448689,7643690,30509378,24902281))

#Get the paths of the .gff files to search through
paths2gff <- list.files(path="./",recursive = T,full.names = T) %>% grep(.,pattern='gff$',value = T) %>%
  as_tibble(.)%>%
  rename('gffPath' = 'value')%>%
  mutate(Species = str_extract(gffPath, "(?<=Genomes/)[^/]+") %>% gsub(pattern='Genome',replacement="",.))

#Join those to the table of genes
Genes2Check <- Genes2Check %>%
  left_join(.,paths2gff,by='Species')

#Write out the gene information for the regions of interest.
cmds2Pass <- list() 
for (l in 1:nrow(Genes2Check)){
  cmds2Pass[[l]] <- paste("awk '$1 == ", '"',Genes2Check[l,2],'"',' && $3 == "gene" && $5 > ', Genes2Check[l,4],' && $4 < ',Genes2Check[l,5], "' ", Genes2Check[l,6]," > ./results/mutations/",Genes2Check[l,1],"_", Genes2Check[l,2],"_",Genes2Check[l,3],".txt",sep="")
  system(cmds2Pass[[l]],wait=T)
  
  }

#### Gene Sequences for each gene with the regions #### 
#BiocManager::install("Biostrings")
dir.create("./results/mutations/GeneSequences")
#library(Biostrings) #don't load biostrings because it masks dplyr::rename and that rename function is set up opposite to dplyrs

#Need all the Genes by species to read in.
GR <- list()
for (l in 1:nrow(Genes2Check)){
  GR[[l]] <- read_delim(paste("./results/mutations/",Genes2Check[l,1],"_", Genes2Check[l,2],"_",Genes2Check[l,3],".txt",sep=""),
                   delim="\t",
                   col_names =F) %>%
    mutate(Name = str_extract(X9, "(?<=Name=)[^;]+"))%>%
    mutate(Species = Genes2Check[l,1])
}

GR <- do.call(rbind,GR)

paths2genomes <- list.files(path="./",recursive = T,full.names = T) %>% grep(.,pattern='genomic.fna$',value = T) %>%
  as_tibble(.)%>%
  dplyr::rename('genomePath' = 'value')%>% #which rename are we using now?
  mutate(Species = str_extract(genomePath, "(?<=Genomes/)[^/]+") %>% gsub(pattern='Genome',replacement="",.))%>%
  filter(!grepl(pattern="data/GCA_",genomePath)) 

for (sp in 1:length(unique(GR$Species))){
Species2subset <- unique(GR$Species)[sp]

SpGenomePath<-paths2genomes %>%
  filter(Species==Species2subset)%>%
  select(genomePath) %>%
  unlist()%>%
  file.path(.)

GRspcs <- GR %>% filter(Species==Species2subset)

#Loop over species
#Load the FASTA file (indexes the file for fast access)
genome <- Biostrings::readDNAStringSet(SpGenomePath)

#names(genome)

#Loop over genes within species
for (gns in 1:nrow(GRspcs)){
Chr <- GRspcs$X1[gns]
GenName <- names(genome)[grep(pattern=Chr,names(genome))]
#Extract the sequence
gene_seq <- Biostrings::subseq(x=genome[[GenName]], start = GRspcs$X4[gns], end = GRspcs$X5[gns])

#Write seq as a text file
Biostrings::writeXStringSet(x = Biostrings::DNAStringSet(gene_seq),
                filepath=paste("./results/mutations/GeneSequences/",Species2subset,"_",GRspcs$X1[gns],"_",GRspcs$Name[gns],".fasta",sep=""),
                format="fasta")
} #over genes within species
} # over species we are interested in

#### Mutations ####

# Now that we a list of all the genes within the larger regions of interst (object GRspcs) we
# can check to see if the polymorphic sites in our species differ by run timing
# group and if those differences lead to changes in the protein coding sequence (CDS)

GRSpGeneChr <- GR %>%
  select(Species,Name,X1)

#Paths to Fsts
SpFsts <- list.files("./results/fst/SpeciesAlign/",recursive = T) %>%
  grep(pattern="\\.fst\\.txt",value=T)%>%
  grep(pattern="Sockeye/euclide|Pink/PhenotypeBased/|early-late_minDepthHalf|Chum/chumrun_",.,value=T)%>%
  grep(pattern="EE-LL",.,value=T,invert = T)%>%
  as_tibble() %>%
  rowwise()%>%
  mutate(Species = lapply(str_split(value,pattern="/"),"[[",1)%>%unlist())%>%
  dplyr::rename('FstPATH' = 'value')%>%
  mutate(Chr =str_extract(FstPATH, "NC_[0-9.]+"))%>%
  mutate(FstPATH = paste("./results/fst/SpeciesAlign/",FstPATH,sep=""))
  

#Path2MAFearly #chromosome specific maf files
PATH2MAFearly <- list.files("./results/fst/SpeciesAlign/",recursive = T) %>%
  grep(pattern="\\.mafs",value=T)%>%
  grep(pattern="early|Early|summer",.,value=T)%>%
  grep(pattern="AEarly_LEven",.,value=T,invert=T) %>%
  as_tibble() %>%
  rowwise()%>%
  mutate(Species = lapply(str_split(value,pattern="/"),"[[",1)%>%unlist())%>%
  dplyr::rename('mafPATHearly' = 'value')%>%
  mutate(Chr = str_extract(mafPATHearly, "NC_[0-9.]+")) %>%
  mutate(mafPATHearly = paste("./results/fst/SpeciesAlign/",mafPATHearly,sep=""))

#Path2MAFlate #chromosome specific maf files
PATH2MAFlate <- list.files("./results/fst/SpeciesAlign/",recursive = T) %>%
  grep(pattern="\\.mafs",value=T)%>%
  grep(pattern="late|Late|Fall|fall",.,value=T)%>%
  grep(pattern="ALate_LEven",.,value=T,invert=T) %>%
  as_tibble() %>%
  rowwise()%>%
  mutate(Species = lapply(str_split(value,pattern="/"),"[[",1)%>%unlist())%>%
  dplyr::rename('mafPATHlate' = 'value')%>%
  mutate(Chr = str_extract(mafPATHlate, "NC_[0-9.]+"))%>%
  mutate(mafPATHlate = paste("./results/fst/SpeciesAlign/",mafPATHlate,sep=""))


#Path2GeneSeq paste together OutputDir, Species, Chr, gene .fasta 
SlmnMutPrtns <- GRSpGeneChr %>%
  dplyr::rename('Chr'='X1')%>%
  left_join(.,SpFsts,by=c("Species","Chr")) %>%
  left_join(.,PATH2MAFearly,by=c("Species","Chr")) %>%
  left_join(.,PATH2MAFlate,by=c("Species","Chr")) %>%
  mutate(Path2GeneSeq = paste("./results/mutations/GeneSequences/",Species,"_",Chr,"_",Name,".fasta",sep=""))%>%
  left_join(.,paths2gff,by='Species')

#Sum Gene Protein Function


SlmnMutPrtns_Out<- list()

source("./scripts/R/SumGeneProtein.R")

for (x in 1:nrow(SlmnMutPrtns)){
  SlmnMutPrtns_Out[[x]] <- SumGeneProtein(Species = SlmnMutPrtns[x,1]%>%unlist(),
                 gene = SlmnMutPrtns[x,2]%>%unlist(),
                 Chrom = SlmnMutPrtns[x,3]%>%unlist(),
                 Path2Fsts = SlmnMutPrtns[x,4]%>%unlist(),
                 Path2MAFearly = SlmnMutPrtns[x,5]%>%unlist(),
                 Path2MAFlate = SlmnMutPrtns[x,6]%>%unlist(), 
                 Path2GeneSeq = SlmnMutPrtns[x,7]%>%unlist(),
                 Path2gff = SlmnMutPrtns[x,8]%>%unlist(),
                 SaveTVrds = T,
                 UseFstThreshold = T)
}

write_delim(x = SlmnMutPrtns_Out %>% bind_rows(),
            file="./results/mutations/MutationsResults_v2.txt",
            delim="\t")

SlmnMutPrtns_Out<-read_delim("./results/mutations/MutationsResults.txt")


SlmnMutPrtns_Out %>%
  bind_rows(.) %>%
  group_by(Species,Gene) %>%
  filter(Type=='CDS') %>%
  slice_head(n=1) %>%
  ungroup()%>%
  summarize(nGenes = length(unique(Gene)),
            nSpeices = length(unique(Species)),
            nSitesGenes = sum(nSitesGene,na.rm=T),
            nSitesGenePoly = sum(nSitesGenePoly,na.rm=T),
            nSitesCDS = sum(nSitesCDS,na.rm=T),
            nSitesCDSPoly = sum(nSitesCDSPoly,na.rm=T),
            nNSynMuts = sum(nNSynmuts,na.rm=T))
  
SlmnMutPrtns_Out %>%
  bind_rows(.) %>%
  group_by(Species,Gene) %>%
  filter(Type=='CDS') %>%
  slice_head(n=1) %>%
  ungroup()%>%
  select(Fst) %>% unlist()%>%
  str_split(.,"-")%>%
  unlist() %>% as.numeric()%>%
  .[!is.na(.)] %>%
  hist(.,main="Fst of Non-synonymous mutations")
  

SlmnMutPrtns_Out %>%
  bind_rows(.) %>%
  group_by(Species,Gene) %>%
  filter(Type=='CDS') %>%
  summarize(AvgnSites = mean(nSitesGene,na.rm=T),
            meannSitesGenePoly = mean(nSitesGenePoly,na.rm=T),
            meannStiesCDS=mean(nSitesCDS,na.rm=T),
            meannSitesCDSPoly = mean(nSitesCDSPoly,na.rm=T),
            meannNSynMuts = mean(nNSynmuts,na.rm=T))%>%
  ungroup()%>%
  group_by(Species) %>%
  summarize(nSitesInGenes = sum(AvgnSites),
            nSitesInGenesPoly = sum(meannSitesGenePoly),
            nSitesInCDS = sum(meannStiesCDS),
            nSitesInCDSPoly = sum(meannSitesCDSPoly),
            nSitesNSynMuts = sum(meannNSynMuts,na.rm=T)) %>%
  ungroup() %>%
  summarize(nSites = sum(nSitesInGenes),
            nSitesInGenesPoly = sum(nSitesInGenesPoly),
            nSitesInCDS = sum(nSitesInCDS),
            nSitesInCDSPoly = sum(nSitesInCDSPoly),
            nSitesNSynMuts = sum(nSitesNSynMuts))



#### Compare TV across species for which there are NonSyn mutations

GenesShrdNSmuts <- c("lrrc9","rtn1a","syne2b")

Files2Read <- list.files("./results/mutations/TranscriptVariants/",full.names = T)%>%
  grep(pattern=paste(GenesShrdNSmuts,collapse="|"),.,value = T)

Files2Read<-list(Files2Read[c(1,5,7)],#lrrc9
                 Files2Read[c(2,6)], #rtn1a
                 Files2Read[c(3,4)]) #syne2b

#Over genes 3
for (Gns in 1:3) {
GenSeq_ShrdNSmuts <- lapply(1:length(Files2Read[[Gns]]), function(x) {
  
  Temp <- readRDS(Files2Read[[Gns]][x]) %>%
    unlist()
  
  Sp <- Files2Read[[Gns]][x] %>% str_split("/|_") %>% lapply(.,"[[",6) %>% unlist()
  Gene <- Files2Read[[Gns]][x] %>% str_split("/|_") %>% lapply(.,"[[",7) %>% unlist()
  
  TempNames <- paste(Sp,Gene,names(Temp),sep=" ")
  
  for(nms in 1:length(TempNames)){
    write.table(x=paste(">",gsub(pattern=" ",replacement="",TempNames[nms]),sep=""),file=paste("./results/mutations/TranscriptVariants/",GenesShrdNSmuts[Gns],"_TranscriptVariants_AllSp.fasta",sep=""),append=T,quote=F,row.names=F,col.names = F)
    write.table(x=Temp[nms],file=paste("./results/mutations/TranscriptVariants/",GenesShrdNSmuts[Gns],"_TranscriptVariants_AllSp.fasta",sep=""),append=T,quote=F,row.names=F,col.names = F)
  } #over the names in the transcript variants
}) #finish the lapply
} #over the Gns genes




