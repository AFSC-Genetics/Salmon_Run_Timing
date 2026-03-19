# SUBSET INDIVIDUALS TO SIMPLIFY PHYLOGENETIC TREE
# SCRIPT FOR LRRC9 ALLELES AND ESRB ALLELES, SEPARATELY
# OUTPUT IS A BAMLIST TO UPLOAD TO USE IN ANGSD IBS RUN

# METHOD
# choose individuals with highest depth of coverage from each 
#  homozygous allele group for each species
#   5 individuals from each groups, therefore 20 per species.
# alternative top 5 individuals for each group and use the outgroups other top 5 for this.

library(here)
library(stringr)
library(tools)
library(tidyverse)

### INPUT FILES #################

#### Metadata ##########################
bam_df <- read.table("./data/bams/fourspp_bamslist_all.txt", header = F)

# metadata
pop_df <- read.csv("./data/raw/fourspecies_runtiming_metadata.csv", header = T)

# depths from each species folder
pink_depth <- read.csv("./results/depth/pink-chum_depths.csv", header = F)
sock_depth <- rbind(read.csv("./results/depth/sock-chum_depths.csv", header = F, sep = "\t"),
                    read.csv("./results/depth/anvb_depths.csv", header = F, sep = "\t"))
chum_depth <- read.csv("./results/depth/chumrun_depths.csv", header = F)
coho_depth <- read.csv("./results/depth/coho-chum_depths.csv", header = F)

#### Genotypes #################

# LRRC9
lrrc9_al <- rbind(read.delim2("./data/R/threespp_lrrc9_alleles/sock-all_lrrc9_alleles_all.txt", sep = "\t", header = F),
                  read.delim2("./data/R/threespp_lrrc9_alleles/chumrun_minInd0.3_minDepthHalf_lrrc9_alleles_all.txt", sep = "\t", header = F),
                  read.delim2("./data/R/threespp_lrrc9_alleles/pink-chum_minInd0.3_lrrc9_alleles_all.txt", sep = "\t", header = F))
colnames(lrrc9_al) <- c("sampleID", "RuntimeAl")

# ESRB
esrb_al <- rbind(read.delim2("./data/R/twospp_esrb_alleles/coho-chum_esrb_alleles.txt", sep = "\t", header = F),
                 read.delim2("./data/R/twospp_esrb_alleles/chumrun_esrb_alleles.txt", sep = "\t", header = F))
colnames(esrb_al) <- c("sampleID", "RuntimeAl")

# ESR1
er1_al <- rbind(read.delim2("./data/R/twospp_esr1_alleles/coho-chum_esr1_alleles.txt", sep = "\t", header = F),
                 read.delim2("./data/R/twospp_esr1_alleles/euclide_esr1_alleles.txt", sep = "\t", header = F))
colnames(er1_al) <- c("sampleID", "RuntimeAl")

#### Edit Input data #####################

# add column with just sampleID
bam_df <- bam_df %>%
  dplyr::mutate(temp = basename(file_path_sans_ext(bam_df$V1)),
                temp = gsub("^[^_]*_", "", temp),  # remove everything after 1st underscore
                sampleID = str_extract(temp, "[^_]+")) %>%    
  select(-temp)

depth_df <- rbind(pink_depth, sock_depth, chum_depth, coho_depth)
colnames(depth_df) <- c("sampleID", "depth")

sample_df <- bam_df %>%
  inner_join(depth_df, by = "sampleID") %>%
  inner_join(pop_df, by = "sampleID")

### LRRC9 ######################

lrrc9_df <- inner_join(sample_df, lrrc9_al, by = "sampleID") %>%
  mutate(sppAl = paste0(Species,"-",RuntimeAl))

lrrc9_top5 <- lrrc9_df %>%
  filter(RuntimeAl != "EL") %>%
  group_by(sppAl) %>%
  top_n(5, depth)

lrrc9_top5 <- lrrc9_top5[order(lrrc9_top5$sppAl, lrrc9_top5$depth, decreasing = T),]

### ESRB ##################################

esrb_df <- inner_join(sample_df, esrb_al, by = "sampleID") %>%
  mutate(sppAl = paste0(Species,"-",RuntimeAl))

esrb_top5 <- esrb_df %>%
  filter(RuntimeAl != "EL") %>%
  group_by(sppAl) %>%
  top_n(5, depth)

esrb_top5 <- esrb_top5[order(esrb_top5$sppAl, esrb_top5$depth, decreasing = T),]

### ESR1 ##################################

er1_df <- inner_join(sample_df, er1_al, by = "sampleID") %>%
  mutate(sppAl = paste0(Species,"-",RuntimeAl))

er1_top5 <- er1_df %>%
  filter(RuntimeAl != "EL") %>%
  group_by(sppAl) %>%
  top_n(5, depth)

er1_top5 <- er1_top5[order(er1_top5$sppAl, er1_top5$depth, decreasing = T),]
View(er1_top5)

### 40 Individuals - 10/spp #############################
# use the top 5 from represented gene for unrepresented gene group
# i.e. for coho at lrrc9, use coho's top 5 from esrb.

# LRRC9
fourspp_lrrc9 <- rbind(lrrc9_top5,
                       filter(esrb_top5, Species == "Coho"))

  write.table(fourspp_lrrc9[,"V1"], file = "./data/R/fourspp_lrrc9_top5_ibs_input.txt",
              quote = F, row.names = F, col.names = F)

# ESRB
fourspp_esrb <- rbind(esrb_top5,
                      filter(lrrc9_top5, Species == "Pink"),
                      filter(lrrc9_top5, Species == "Sockeye"))
  
  write.table(fourspp_esrb[,"V1"], file = "./data/R/fourspp_esrb_top5_ibs_input.txt",
              quote = F, row.names = F, col.names = F)

# ER1
fourspp_er1 <- rbind(er1_top5,
                     filter(lrrc9_top5, Species == "Pink"),
                     filter(lrrc9_top5, Species == "Chum"))
   
  write.table(fourspp_er1[,"V1"], file = "./data/R/fourspp_esr1_top5_ibs_input.txt",
              quote = F, row.names = F, col.names = F)
