# Create table with chromosome positions for each FST value across all species 
# for the two regions of interest

library(tidyverse)
library(magrittr)
install.packages("biomaRt")
library(BiocManager)
BiocManager::install("biomaRt")
library(biomaRt)

# reference genome files
GFFDIR <- "./data/R/genomic.gff"
GAFDIR <- "./data/R/GCF_023373465.1-RS_2023_03_gene_ontology.gaf"

# table for chromosome name and number
chrom_df <- read.table("./data/R/chrom_meta.txt", header = TRUE)

################ GFF FILE #############################################
# find exons from gff file for genes of interest (from NCBI chum reference genome)
gff_df <- read.delim(GFFDIR, header = F, comment.char = "#")
  gff_df <- gff_df[,c(1:5,9)] # remove excess columns
    colnames(gff_df) <- c("chrName", "RefSeq","exon","start.pos","fin.pos", "ID") # rename remaining columns

# only keep chr29 & 35
gff_subset <- gff_df %>%
  lst(chrom_df) %>%
  reduce(left_join) %>%
  filter(chr == "29" & fin.pos >= 24.7*10^6 & start.pos <= 26.4*10^6
           | chr == "35" & fin.pos >= 27.86*10^6 & start.pos <= 28.3*10^6,
         exon == "gene") %>%
  dplyr::select(-c(chrName, RefSeq, exon)) %>% dplyr::select(chr, everything())
rm(gff_df)

# prep pattern for str_match below
gene_pattern <- "gene=\\s*(.*?)\\s*;"    # keep string btwn "gene=" & ":product" 
exon_pattern <- "ID=exon-\\s*(.*?)\\s*;Parent"  # keep string btwn "exon=" & ";Parent" 
descr_pattern <- ";description=\\s*(.*?)\\s*;"

# for first creating gene file
gene_write <- gff_subset %>% 
  mutate(geneID = str_match(ID, gene_pattern)[,2],
         geneName = str_match(ID, descr_pattern)[,2])

rm(gene_pattern, exon_pattern, descr_pattern)

########### IMPORT GENE ASSOCIATION FILE  ######################################

go_df <- read.delim(GAFDIR, header = F, comment.char = "!")
  go_df <- go_df[,c(3:5,9:10)]
  colnames(go_df) <- c("geneID", "qualifier", "go_id", "ontology", "geneProduct")

GO_subset <- go_df %>%
  filter(geneID %in% unique(gene_write$geneID),
         ontology != "C") 

coho_ensembl <- useMart(biomart = "ensembl", dataset = "okisutch_gene_ensembl")

go_fxn <- getBM(attributes = c("go_id", "name_1006"),
                filters = "go",
                values = unique(GO_subset$go_id),
                mart = coho_ensembl)
rm(coho_ensembl, go_df)

GO_subset1 <- GO_subset %>% 
  left_join(go_fxn, by = "go_id") %>%
  filter(!is.na(name_1006)) %>% # not all GO terms have annotations
  mutate(qualifier = gsub("_"," ",qualifier),
         go_Annotation = paste0(ontology,": ",qualifier," ",name_1006)) %>%
  group_by(geneID) %>%
  dplyr::summarize(gene_fxns = ifelse(all(is.na(go_Annotation)), NA_character_, paste(na.omit(go_Annotation), collapse = "; "))) %>%
  ungroup()

geneFxn <- gene_write %>%
  left_join(GO_subset1, by = "geneID")
rm(GO_subset, GO_subset1, go_fxn)

################    Read in each species FST file     #############################

########## PINK
pink_Fst <- read.delim2("./results/fst/allele/pink-chum_NC_068455.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")
colnames(pink_Fst) <- c("region", "chrName", "midPos", "Nsites", "Pink.Fst")
pink_Fst <- pink_Fst[,c(2:3,5)] 

########## SOCKEYE
sock_Fst <- read.delim("./results/fst/sock-all_NC_068455.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                       row.names = NULL,sep = "\t")
  colnames(sock_Fst) <- c("region", "chrName", "midPos", "Nsites", "Sockeye.Fst")
  sock_Fst <- sock_Fst[,c(2:3,5)]

########## COHO
coho_Fst <- read.delim2("./results/fst/allele/coho-chum_NC_068449.1_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")
  colnames(coho_Fst) <- c("region", "chrName", "midPos", "Nsites", "Coho.Fst")
  coho_Fst <- coho_Fst[,c(2:3,5)]

########## CHUM
chum29 <- read.delim2("./results/fst/allele/chumrun_NC_068449.1_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")
chum35 <- read.delim2("./results/fst/allele/chumrun_NC_068455.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")
chum_Fst <- rbind(chum29, chum35)
  colnames(chum_Fst) <- c("region", "chrName", "midPos", "Nsites", "Chum.Fst")
  chum_Fst <- chum_Fst[,c(2:3,5)]

  rm(chum29, chum35)
  
############## COMBINE FOUR SPECIES INTO ONE DATAFRAME #######################
# bind rows together
four_df <- lst(pink_Fst, sock_Fst, chum_Fst, coho_Fst) %>% 
  reduce(full_join) %>%
  lst(chrom_df) %>%
  reduce(left_join) %>%
  dplyr::select(chr, everything()) %>% dplyr::select(-chrName) %>%
  filter(Pink.Fst > 0 | Sockeye.Fst > 0 | Chum.Fst > 0 | Coho.Fst > 0) %>%
  mutate(HighFst = case_when(
    as.numeric(Pink.Fst) > 0.5 | Sockeye.Fst > 0.5 | Chum.Fst > 0.5 | Coho.Fst > 0.5 ~ T,
    TRUE ~ NA))

rm(pink_Fst, sock_Fst, coho_Fst, chum_Fst)

# Change negative Fst SNPs to 0s
four_df[four_df < 0] <- 0

chr29 <- four_df %>%
  filter(chr == 29,
         as.numeric(midPos) > 24.7*10^6,
         as.numeric(midPos) < 26.4*10^6) 

chr35 <- four_df %>%
  filter(chr == 35,
         midPos > 27.86*10^6,
         midPos < 28.3*10^6) 

four_subset <- rbind(chr29, chr35); rm(chr29, chr35, four_df)

table_df <- four_subset %>%
  left_join(geneFxn, by = "chr", relationship = "many-to-many") %>%
  # Filter rows based on gene range
  mutate(gene = ifelse(midPos >= start.pos & midPos <= fin.pos, geneID, NA),
         geneName = ifelse(midPos >= start.pos & midPos <= fin.pos, geneName, NA),
         geneAnnotation = ifelse(midPos >= start.pos & midPos <= fin.pos, gene_fxns, NA)) %>%
  # Group by 'chr' and 'position' to avoid duplicates
  group_by(chr, midPos) %>%
  # remove duplicates by summarizing & keeping the first non-NA gene
  dplyr::summarize(gene = ifelse(all(is.na(gene)), NA_character_, paste(na.omit(gene), collapse = "; ")),
            geneName = ifelse(all(is.na(geneName)), NA_character_, paste(na.omit(geneName), collapse = "; ")),
            geneAnnotation = ifelse(all(is.na(geneAnnotation)), NA_character_, paste(na.omit(geneAnnotation), collapse = "; "))) %>%
  ungroup() %>%
  left_join(four_subset, by = c("chr", "midPos")) %>% # re-add FST columns
  mutate(sharedSNPs = rowSums(!is.na(dplyr::select(., Pink.Fst:Coho.Fst))) > 1,
         sharedSNPs = ifelse(sharedSNPs == T, "Multiple Spp", NA),
         geneName = gsub("%","-",geneName)) %>%
  relocate(gene, geneName, geneAnnotation, .after = last_col())

write.csv(table_df, "./data/R/supplemental_table2_FST_and_genes.csv", row.names = F)

