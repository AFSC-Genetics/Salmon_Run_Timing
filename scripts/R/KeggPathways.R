#Kegg Pathways

Files2Kegg <- list.files("./results/mutations/",full.names = T)%>%
  grep(pattern="CDS",.,value = T)

Genes2Kegg <- lapply(1:length(Files2Kegg),function(x) 
  read_delim(Files2Kegg[x],delim="\t",col_names = F) %>%
    filter(X3=='gene')%>%
    summarize(GeneID = str_extract(X9,"(?<=GeneID:)[^;]+"),
              Name = str_extract(X9,"(?<=Name=)[^;]+"),
              Gene = str_extract(X9,"(?<=gene=)[^;]+"),
              GeneDesc = str_extract(X9,"(?<=description=)[^;]+"),
              Species = sub("_.*", "",basename(Files2Kegg[x])))
  )%>%
  do.call(rbind,.) %>%
  mutate(Prefix = case_when(
    Species == 'Chum' ~ 'oke',
    Species == 'Pink' ~ 'ogo',
    Species == 'Sockeye' ~ 'one',
    Species == 'Coho' ~ 'oki',
    TRUE ~ NA_character_))%>%
  mutate(KeggID = paste(Prefix,GeneID,sep=":"))

#remove uncharacterized genes
  # 109893920 = ncRNA https://www.ncbi.nlm.nih.gov/gene/109893920
  # 109893921 = ncRNA http://ncbi.nlm.nih.gov/gene/?term=109893921
  # 116356176 = ncRNA https://www.ncbi.nlm.nih.gov/gene/?term=116356176
  # 115139260 = lncRNA https://www.ncbi.nlm.nih.gov/gene/115139260

Genes2Kegg_cln <- Genes2Kegg %>%
  filter(!GeneID %in% c("109893920","109893921","115139260","116356176"))

#BiocManager::install("KEGGREST")
Kegg_KOs <- sapply(1:nrow(Genes2Kegg_cln), function(x) 
                KEGGREST::keggGet(Genes2Kegg_cln$KeggID[x]) %>%
               {.[[1]]$ORTHOLOGY} %>%
                 {names(.)[1]}
)

names(Kegg_KOs)<-Genes2Kegg_cln$GeneID
  
Kegg_KOsTib <- enframe(Kegg_KOs, 
        name = "entry_name", 
        value = "entry_value") %>%
  unnest_longer(entry_value) %>%
  rename('GeneID' = 'entry_name') %>%
  full_join(., Genes2Kegg_cln %>% select(GeneID,Name,Gene,GeneDesc,Species))

# Manual lookup
Harbi1 <- 'K23222'
Syne2 <- 'K19346' 
#lrrc9 not found
#pcnx4 not found, K28396 is pecanex-like protein
#protein_LBH not found,
exocyst_complex_component_3_like_protein <-'K19987' 
#uncharacterized si:dkey-229e3.2 looks to be ZDB-GENE-131121-296 in zebrafish
#MYCT1 not found
PAK6 <- 'K05735'  #serine/threonine-protein kinase PAK 6-like == LOC115138567
                             
Kegg_KOsTib <- Kegg_KOsTib %>%
  mutate(entry_value = case_when(
    Name == 'syne2b' ~ 'K19346',
    Name == 'LOC127915932' ~ 'K23222',
    Name == 'LOC109893917' ~ 'K19987',
    Name == 'LOC115138567' ~ 'K05735', #serine/threonine-protein kinase PAK 6-like / PAK6 == K05735
    .default = entry_value  # Use the column name directly here
  )) 

# We could pull the pathways by hand by 
# create a 2 col. file that has gene name, K# and paste into the reconstruct tool here:
# https://www.genome.jp/kegg/mapper/reconstruct.html

# pull from the API

mykos<-Kegg_KOsTib$entry_value[!is.na(Kegg_KOsTib$entry_value)] %>%
  gsub(pattern="K",replacement="ko:K") 

KEGG_PathMap <-  KEGGREST::keggLink("pathway", "ko") %>%
  {data.frame(
    ko = names(.),
    pathway = as.character(.)
  )} %>%
  filter(ko %in% mykos) %>%
  group_by(pathway) %>%
  summarise(
    count = n(),
    genes = paste(ko, collapse = ", ")
  )

target_pathways <- KEGG_PathMap %>%
  filter(grepl('ko', pathway)) %>%
  pull(pathway)

get_kegg_info <- function(ids) {
  chunks <- split(ids, ceiling(seq_along(ids) / 10)) # Split IDs into chunks of 10
  
  map_df(chunks, function(chunk) {
    res <-  KEGGREST::keggGet(chunk)
    map_df(res, ~tibble(
      pathway = .x$ENTRY,
      Class = ifelse(!is.null(.x$CLASS), .x$CLASS, NA),
      PathName = ifelse(!is.null(.x$NAME), .x$NAME, NA)
    ))
  })
}

KEGG_Info_Table <- get_kegg_info(target_pathways)

KEGG_Final <- KEGG_PathMap %>%
  left_join(KEGG_Info_Table %>% mutate(pathway = paste("path:",pathway,sep="")) , by = "pathway") %>%
  filter(!grepl(pattern="map", pathway)) #these are redundant

#Now Join the better names of the genes than their ko numbers
Kegg_KOsTib <- Kegg_KOsTib %>%
  mutate(GeneName = ifelse(!grepl(pattern="LOC",Gene),Gene,GeneDesc))%>%
  mutate(GeneName = ifelse(is.na(GeneName),Gene,GeneName))

KEGG_Final_Mapped <- KEGG_Final %>%
  separate_rows(genes, sep = ",\\s*") %>% 
  mutate(ko_clean = gsub("ko:", "", genes)) %>%
  left_join(Kegg_KOsTib, by = c("ko_clean" = "entry_value")) %>%
  group_by(pathway, Class, PathName,count) %>% # Add any other columns you want to keep
  summarise(
    GeneNames_List = paste(unique(GeneName), collapse = ", "),
    .groups = "drop"
  )
  

write_delim(KEGG_Final_Mapped,file="./KeggRes.txt",delim="\t")
