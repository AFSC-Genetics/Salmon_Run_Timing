# Create allele-base NJ Tree for esrb and lrrc9 data


library(BiocManager)
BiocManager::install("ggtree", force = T)
library(ggtree)
library(treeio)
packages_needed <- c("tidyverse", "tools", "cowplot", "ape", "ggstance", 
                     "pals", "phytools", "magrittr")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

##### INPUT DATA #########################
pop_df <- read.csv("./data/raw/fourspecies_runtiming_metadata.csv", header = T)

# call in allele data and add to sample data
lrrc9_df <- rbind(read.delim2("./data/R/threespp_lrrc9_alleles/sock-all_lrrc9_alleles_all.txt", sep = "\t", header = F),
                  read.delim2("./data/R/threespp_lrrc9_alleles/chumrun_minInd0.3_minDepthHalf_lrrc9_alleles_all.txt", sep = "\t", header = F),
                  read.delim2("./data/R/threespp_lrrc9_alleles/pink-chum_minInd0.3_lrrc9_alleles_all.txt", sep = "\t", header = F))
  colnames(lrrc9_df) <- c("sampleID", "RuntimeAl")

esrb_df <- rbind(read.delim2("./data/R/twospp_esrb_alleles/coho-chum_esrb_allele_groups.txt", sep = "\t", header = F),
                 read.delim2("./data/R/twospp_esrb_alleles/chumrun_esrb_alleles.txt", sep = "\t", header = F))
  colnames(esrb_df) <- c("sampleID", "RuntimeAl")

options(ignore.negative.edge=TRUE)

######################  LRRC9: TOP FIVE PER GROUP   ###########################

ibs_mat_lrrc9 <- as.matrix(read.table(paste0("./results/ibs/fourspp_top5_GL1_MM4_lrrc9.ibsMat")))
bam_lrrc9 <- read.table(paste0("./data/R/fourspp_lrrc9_top5_ibs_input.txt"), header = F)

# convert bam to FID with ABLG
FIDlrrc9 <- bam_lrrc9 %>%
  dplyr::mutate(FID = row_number(),
                temp = basename(file_path_sans_ext(V1)),      # remove file path
                temp = gsub("^[^_]*_", "", temp),             # remove everything after 1st underscore
                sampleID = str_extract(temp, "[^_]+")) %>%
  select(FID, sampleID)

sample_table <- left_join(pop_df, lrrc9_df, by = "sampleID") %>%
  inner_join(FIDlrrc9, by = "sampleID")

colnames(ibs_mat_lrrc9) <- sample_table$sampleID # set IBS matrix colnames to sampleID

lrrc9AlleleTree <- nj(ibs_mat_lrrc9) %>%
  phytools::midpoint.root() %>%
  as.treedata() %>%
  dplyr::mutate(FID = row_number()) %>%
  left_join(sample_table, by = "FID")

legend <- get_legend(ggtree(lrrc9AlleleTree, aes(color = RuntimeAl), size = 1.5) + 
                       scale_color_manual(name = "Genotype", values = c("goldenrod2", "royalblue2")) +
                       theme(legend.title = element_text(size = 22),
                             legend.text = element_text(size = 20)))
plot(legend)

# look at node numbers to choose colors below
ggtree(lrrc9AlleleTree, aes(color = RuntimeAl), size = 1.5) + 
  scale_color_manual(name = "lrrc9", values = c("goldenrod2", "royalblue2")) +
  ggnewscale::new_scale_colour() +
  geom_tiplab(align = T, aes(color = Species)) +
  scale_color_manual(name = "Species", values = palette.colors(palette = "R4")[1:4]) +
  geom_treescale(x = 0.05, y = 1) +
  geom_text(aes(label=node), hjust=0)

## binary choices of colors
d_lrrc9 <- data.frame(node=1:Nnode2(lrrc9AlleleTree), colour = "black")

# choose allele colors
d_lrrc9[c(16:20,76:79), 2] <- "goldenrod2" # pink
d_lrrc9[c(11:15,72:75), 2] <- "royalblue2" # pink
d_lrrc9[c(6:10,58:61), 2] <- "goldenrod2" # sockeye
d_lrrc9[c(1:5,54:57), 2] <- "royalblue2" # sockeye
d_lrrc9[c(26:30,67:70), 2] <- "goldenrod2" # chum
d_lrrc9[c(21:25,63:66), 2] <- "royalblue2" # chum
d_lrrc9[c(31:40,42:49), 2] <- "gray40" # coho NA

## multiple choices of line types
p_lrrc9 <- ggtree(lrrc9AlleleTree, size = 1.5)

c_lrrc9 <- flip(p_lrrc9, 72, 76) %>% 
  flip(53, 50) %>%
  scaleClade(node=50, scale = 0.3)

finalplot_lrrc9 <- c_lrrc9 %<+% d_lrrc9 + aes(colour=I(colour))+
  geom_treescale(x = 0.03, y = 6, fontsize = 5)
finalplot_lrrc9

ggsave(plot = finalplot_lrrc9, filename = "./figures/tree/fourspp_lrrc9_tree_top5_flip2.pdf",
       width = 15, height = 15, limitsize = F)

# visualize branch lengths
ggtree(lrrc9AlleleTree, aes(color = RuntimeAl), size = 1.5) + 
  ggnewscale::new_scale_colour() +
  geom_text2(aes(label = branch.length), hjust = 0.5)

rm(ibs_mat_lrrc9, bam_lrrc9, lrrc9AlleleTree, c_lrrc9, d_lrrc9, p_lrrc9, finalplot_lrrc9)


######################  ESRB: TOP FIVE PER GROUP   ###########################

ibs_mat_esrb <- as.matrix(read.table(paste0("./results/ibs/fourspp_GL1_MM4_3_top5_esrb_expand.ibsMat")))
bam_esrb <- read.table(paste0("./data/R/fourspp_esrb_top5_ibs_input.txt"), header = F)

# convert bam to FID with ABLG
FID_esrb <- bam_esrb %>%
  dplyr::mutate(FID = row_number(),
                temp = basename(file_path_sans_ext(V1)),      # remove file path
                temp = gsub("^[^_]*_", "", temp),             # remove everything after 1st underscore
                sampleID = str_extract(temp, "[^_]+")) %>%
  select(FID, sampleID)

sample_table <- left_join(pop_df, esrb_df, by = "sampleID") %>%
  inner_join(FID_esrb, by = "sampleID")

colnames(ibs_mat_esrb) <- sample_table$sampleID # set IBS matrix colnames to sampleID

esrbAlleleTree <- nj(ibs_mat_esrb) %>%
  phytools::midpoint.root() %>%
  as.treedata() %>%
  dplyr::mutate(FID = row_number()) %>%
  left_join(sample_table, by = "FID")

## multiple choices of line types
p_esrb <- ggtree(esrbAlleleTree, size = 1.5)
ggtree(esrbAlleleTree, aes(color = RuntimeAl), size = 1.5) + 
  ggnewscale::new_scale_colour() +
  geom_tiplab(align = T, aes(color = Species)) +
  geom_treescale(x = 0.1, y = 1) +
  geom_text(aes(label=node), hjust=0)

## binary choices of species colors
d_esrb <- data.frame(node=1:Nnode2(esrbAlleleTree), colour = "black")

# choose allele colors
d_esrb[c(6:10,47:50), 2] <- "goldenrod2" # coho
d_esrb[c(1:5,42:45), 2] <- "royalblue2" # coho
d_esrb[c(16:20,63:66), 2] <- "goldenrod2" # chum
d_esrb[c(11:15,67:70), 2] <- "royalblue2" # chum
d_esrb[c(21:30,31:40), 2] <- "gray40" # sockeye and pink NA

c_esrb <- flip(p_esrb, 63, 67) %>% 
  flip(52, 46) %>%
  scaleClade(node=52, scale = 0.3) %>%
  scaleClade(node=71, scale = 0.3)

finalplot_esrb <- c_esrb %<+% d_esrb + aes(colour=I(colour)) +
  geom_treescale(x = 0.04, y = 4, fontsize=5)
finalplot_esrb

ggsave(plot = finalplot_esrb, filename = "./figures/tree/fourspp_esrb_tree_top5_flip2.pdf",
       width = 15, height = 15, limitsize = FALSE)

# visualize branch lengths
ggtree(esrbAlleleTree, aes(color = RuntimeAl), size = 1.5) + 
  ggnewscale::new_scale_colour() +
  geom_text2(aes(label = branch.length), hjust = -.3)

rm(ibs_mat_esrb, bam_esrb, esrbAlleleTree, c_esrb, d_esrb, p_esrb, finalplot_esrb)

