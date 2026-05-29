# COMBINE PCAs and FST for LRRC9, ESRB, and ESR1 
# FIGURE 3:
#  1A) LRRC9 COMBINE PCAS 
#  1B) ESRB COMBINE PCAS
#  1C) ESR1 COMBINE PCAS

# this part is unfinished for esr1 because they might be split figures !
#  2A) LRRC9 ALLELE FST
#  2B) ESRB ALLELE FST
#  2C) ESR1 ALLELE FST

packages_needed <- c("ggplot2", "scales", "ggpubr", "ggrepel", "stringr", 
                     "data.table", "plyr","tools","gtools","reshape2", 
                     "patchwork", "cowplot", "tidyverse")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}
rm(i, packages_needed)

# setwd("../../Salmon_Run_Timing")

# metadata
allmeta_df <- read.csv("./data/raw/fourspecies_runtiming_metadata.csv", header = T)

# Plotting
mypalette <- c("goldenrod1", "royalblue3")

theme_set(
  theme( 
    legend.text=element_text(size=16), legend.title = element_text(size = 18),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text = element_text(angle = 0, size = 14), axis.title = element_text(size = 16),
    legend.position = "top", panel.background = element_rect(fill = "white"), 
    panel.spacing = unit(0,"lines"), strip.text.y = element_text(angle = 0)
  )
)

#### 1A) LRRC9 #################################################################

##### Pink Even ######################

# read in the covariance matrix
pink_cov <- as.matrix(read.table("./results/pca/pink-chum_NC_068455.1_lrrc9_minInd0.3.cov"))
  pink_e <- eigen(pink_cov) # calculate eigenvalues
  pink_e_vectors <- as.data.frame(pink_e$vectors)
  # determine the variance explained as a percent
  pink_varPC1 <- (pink_e$values[1]/sum(pink_e$values))*100 # PC1 variance
  pink_varPC2 <- (pink_e$values[2]/sum(pink_e$values))*100 # PC2 variance
  
# call in bams & metadata
pink_bam_df <- read.table("./data/bams/pink-chum_filtered_bamslist.txt", header = F) %>%
  mutate(ABLG = as.numeric(gsub("[^0-9]","", V1))) %>% select(ABLG)

pink_meta <- allmeta_df %>%
  filter(Species == 'Pink') %>%
  mutate(ABLG = as.numeric(sub('ABLG','',sampleID))) %>%
  inner_join(pink_bam_df, .)

# combine row names (population info) with the covariance matrix
pink_pca.vectors = as_tibble(cbind(pink_meta, pink_e_vectors))

# create dataframe for color designation
pink_pca.vectors$Runtime <- factor(pink_pca.vectors$Runtime, levels = c("Early", "Late"))
  names(mypalette) <- levels(pink_pca.vectors$Runtime)

### Plot
pink_lrrc9 <- ggplot(data = pink_pca.vectors, 
                     aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = 3, color = "gray20") +
  scale_fill_manual(name = expression('Pink'~italic(lrrc9)), values = mypalette) +
  scale_shape_manual(name = expression('Pink'~italic(lrrc9)), values = c(21,22)) +
  scale_x_continuous(breaks = c(-0.05, 0, 0.05, 0.1, 0.15)) +
  geom_vline(xintercept = 0.07, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle("Pink") + 
  #ggtitle(expression('Pink'~italic(lrrc9))) +
  labs(x = paste0("PC1 (",round(pink_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(pink_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5, margin=margin(0,0,10,0)))

##### Sockeye - Euclide ####################

# read in the covariance matrix
sock_cov <- as.matrix(read.table("./results/pca/euclide_NC_068455.1_minInd0.3_lrrc9.cov"))
  sock_e <- eigen(sock_cov) # calculate eigenvector values
  sock_e_vectors <- as.data.frame(sock_e$vectors)
  # determine the variance explained as a percent
  sock_varPC1 <- (sock_e$values[1]/sum(sock_e$values))*100 # PC1 variance
  sock_varPC2 <- (sock_e$values[2]/sum(sock_e$values))*100 # PC2 variance

# call in some metadata  
sock_bam_df <- read.table("./data/bams/euclide_bamslist.txt", header = F) %>%
  dplyr::mutate(temp = basename(file_path_sans_ext(V1)),
                temp = gsub("^[^_]*_", "", temp), # remove everything after 1st underscore
                sampleID = str_extract(temp, "[^_]+")) %>%
  select(sampleID)

sock_meta <- allmeta_df %>%
  filter(Species == "Sockeye",
         Runtime != "Late Stream") %>% # remove whitefish
  mutate(Runtime = if_else(Runtime == 'Late Beach', 'Late', 'Early')) %>%
  inner_join(sock_bam_df, ., by = "sampleID")

# combine row names (population info) with the covariance matrix
sock_pca.vectors = as_tibble(cbind(sock_meta, sock_e_vectors))

# create dataframe for color designation
sock_pca.vectors$Runtime <- factor(sock_pca.vectors$Runtime, levels = c("Early", "Late"))
  names(mypalette) <- levels(sock_pca.vectors$Runtime)

### SOCKEYE
sock_lrrc9 <- ggplot(data = sock_pca.vectors, 
                     aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = 3, color = "gray20") +
  scale_fill_manual(name = expression('Sockeye'~italic(lrrc9)), values = mypalette) +
  scale_shape_manual(name = expression('Sockeye'~italic(lrrc9)), values = c(21,22)) +
  geom_vline(xintercept = 0.1, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle("Sockeye") + 
  #ggtitle(expression('Sockeye'~italic(lrrc9))) +
  scale_x_reverse(breaks = c(-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15)) +
  labs(x = paste0("PC1 (",round(sock_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(sock_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5, margin=margin(0,0,10,0)))

##### Chum ###################

# read in the covariance matrix
chum_cov <- as.matrix(read.table("./results/pca/chumrun_NC_068455.1_lrrc9_minInd0.3_minDepthHalf.cov"))
  chum_e <- eigen(chum_cov) # calculate eigenvector values
  chum_e_vectors <- as.data.frame(chum_e$vectors)
  # determine the variance explained as a percent
  chum_varPC1 <- (chum_e$values[1]/sum(chum_e$values))*100 # PC1 variance
  chum_varPC2 <- (chum_e$values[2]/sum(chum_e$values))*100 # PC2 variance
  
# call in bams & metadata
chum_bam_df <- read.table("./data/bams/chumrun_bamslist.txt", header = F) %>%
  mutate(ABLG = as.numeric(gsub('[^0-9]','', V1))) %>% select(ABLG)

chum_meta <- allmeta_df %>%
  filter(Species == 'Chum') %>%
  mutate(ABLG = as.numeric(sub('ABLG','',sampleID))) %>%
  inner_join(chum_bam_df, ., by = "ABLG")

##combine row names (population info) with the covariance matrix
chum_pca.vectors = as_tibble(cbind(chum_meta, chum_e_vectors))

# plot colors
chum_pca.vectors$Runtime <- factor(chum_pca.vectors$Runtime, levels = c("Early", "Late"))
  names(mypalette) <- levels(chum_pca.vectors$Runtime)

### CHUM PLOT
chum_lrrc9 <- ggplot(data = chum_pca.vectors, 
                     aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = 3, color = "gray20") +
  scale_fill_manual(name = expression('Chum'~italic(lrrc9)), values = mypalette) +
  scale_shape_manual(name = expression('Chum'~italic(lrrc9)), values = c(21,22)) +
  geom_vline(xintercept = 0.1, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  # ggtitle(expression('Chum'~italic(lrrc9))) +
  ggtitle('Chum') +
  scale_x_reverse(breaks = c(-0.05, 0, 0.05, 0.1, 0.15)) +
  labs(x = paste0("PC1 (",round(chum_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(chum_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5, margin=margin(0,0,10,0)))

#### 1B) ESRB ###################################################################

###### COHO ###########

# read in the covariance matrix
coho_cov <- as.matrix(read.table("./results/pca/coho-chum_NC_068449.1_s25414060_e25501622_esrb_minInd0.3_minDepthHalf.cov"))
  coho_e <- eigen(coho_cov)
  coho_e_vectors <- as.data.frame(coho_e$vectors)
  # determine the variance explained as a percent
  coho_varPC1 <- (coho_e$values[1]/sum(coho_e$values))*100 # PC1 variance
  coho_varPC2 <- (coho_e$values[2]/sum(coho_e$values))*100 # PC2 variance
  
# call in metadata
coho_bam_df <- read.table("./data/bams/coho-chum_filtered_bamslist.txt", header = F) %>%
  mutate(ABLG = as.numeric(gsub("[^0-9]", "", V1))) %>% select(-V1)

coho_meta <- allmeta_df %>%
  filter(Species == 'Coho') %>%
  mutate(ABLG = as.numeric(gsub('ABLG','',sampleID))) %>%
  inner_join(coho_bam_df, ., by = "ABLG")

# combine row names (population info) with the covariance matrix
coho_pca.vectors = as_tibble(cbind(coho_meta, coho_e_vectors))

coho_pca.vectors$Runtime <- factor(coho_pca.vectors$Runtime, levels = c("Early", "Late"))
  names(mypalette) <- levels(coho_pca.vectors$Runtime)

### COHO
coho_esrb <- ggplot(data = coho_pca.vectors, 
                    aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = 3, color = "gray20") +
  scale_fill_manual(name = expression('Coho'~italic(esrb)), values = mypalette) +
  scale_shape_manual(name = expression('Coho'~italic(esrb)), values = c(21,22)) +
  scale_x_continuous(breaks = c(-.1,-0.05, 0, 0.05, 0.1, 0.15)) +
  geom_vline(xintercept = 0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  # ggtitle(expression('Coho'~italic(esrb))) +
  ggtitle('Coho') +
  theme(legend.position = "bottom") +
  labs(x = paste0("PC1 (",round(coho_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(coho_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5, margin=margin(0,0,10,0)))

##### CHUM ####################

# read in the covariance matrix
chumE_cov <- as.matrix(read.table("./results/pca/chumrun_NC_068449.1_s25414060_e25501622_esrb_minInd0.3_minDepthHalf.cov"))
  chumE_e <- eigen(chumE_cov)
  chumE_e_vectors <- as.data.frame(chumE_e$vectors)
  # determine the variance explained as a percent
  chumE_varPC1 <- (chumE_e$values[1]/sum(chumE_e$values))*100 # PC1 variance
  chumE_varPC2 <- (chumE_e$values[2]/sum(chumE_e$values))*100 # PC2 variance

# combine row names (population info) with the covariance matrix
chumE_pca.vectors = as_tibble(cbind(chum_meta, chumE_e_vectors))

# plot colors
chumE_pca.vectors$Runtime <- factor(chumE_pca.vectors$Runtime, levels = c("Early", "Late"))
  names(mypalette) <- levels(chumE_pca.vectors$Runtime)

### CHUM PLOT
chum_esrb <- ggplot(data = chumE_pca.vectors, 
                    aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = 3, color = "gray20") +
  theme(legend.position = "bottom") +
  scale_fill_manual(name = 'Run Timing', values = mypalette) +
  scale_shape_manual(name = 'Run Timing', values = c(21,22)) +
  scale_x_continuous(breaks = c(-0.1,-0.05, 0, 0.05, 0.1, 0.15)) +
  geom_vline(xintercept = 0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  # ggtitle(expression('Chum'~italic(esrb))) +
  ggtitle('Chum') +
  labs(x = paste0("PC1 (",round(chumE_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(chumE_varPC2, digits = 1),"%)")) +
  theme(legend.position = "right",
        legend.text=element_text(size=32), legend.title = element_text(size=35),
        plot.title = element_text(size=20, hjust=0.5, margin=margin(0,0,10,0))) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5,
                             override.aes = list(size=8)))

#### 1C) ESR1 #################################################################
###### COHO ###########

# read in the covariance matrix
cohoE1_cov <- as.matrix(read.table(paste0("./results/pca/coho-chum_NC_068428.1_s23264160_e23645630_minInd0.3_minDepthHalf_esr1.cov")))
  cohoE1_e <- eigen(cohoE1_cov)
  cohoE1_e_vectors <- as.data.frame(cohoE1_e$vectors)
  # determine the variance explained as a percent
  cohoE1_varPC1 <- (cohoE1_e$values[1]/sum(cohoE1_e$values))*100 # PC1 variance
  cohoE1_varPC2 <- (cohoE1_e$values[2]/sum(cohoE1_e$values))*100 # PC2 variance

##combine row names (population info) with the covariance matrix
cohoE1_pca.vectors = as_tibble(cbind(coho_meta, cohoE1_e_vectors))

cohoE1_pca.vectors$Runtime <- factor(cohoE1_pca.vectors$Runtime, levels = c("Early", "Late"))
  names(mypalette) <- levels(cohoE1_pca.vectors$Runtime)

### COHO
coho_er1 <- ggplot(data = cohoE1_pca.vectors, 
                   aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = 3, color = "gray20") +
  scale_fill_manual(name = expression('Coho'~italic(er1)), values = mypalette) +
  scale_shape_manual(name = expression('Coho'~italic(er1)), values = c(21,22)) +
  scale_x_continuous(breaks = c(-.15,-.1,-0.05, 0, 0.05, 0.1, 0.15)) +
  geom_vline(xintercept = 0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  # ggtitle(expression('Coho'~italic(esr1))) +
  ggtitle('Coho') +
  theme(legend.position = "bottom") +
  labs(x = paste0("PC1 (",round(cohoE1_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(cohoE1_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5, margin=margin(0,0,10,0)))

##### SOCKEYE EXPAND ####################

# read in the covariance matrix
sockE1_cov <- as.matrix(read.table(paste0("./results/pca/euclide_NC_068428.1_s23264160_e23645630_minInd0.3_esr1.cov")))
  sockE1_e <- eigen(sockE1_cov)
  sockE1_e_vectors <- as.data.frame(sockE1_e$vectors)
  # determine the variance explained as a percent
  sockE1_varPC1 <- (sockE1_e$values[1]/sum(sockE1_e$values))*100 # PC1 variance
  sockE1_varPC2 <- (sockE1_e$values[2]/sum(sockE1_e$values))*100 # PC2 variance

##combine row names (population info) with the covariance matrix
sockE1_pca.vectors = as_tibble(cbind(sock_meta, sockE1_e_vectors))

# plot colors
sockE1_pca.vectors$Runtime <- factor(sockE1_pca.vectors$Runtime, levels = c("Early", "Late"))
  names(mypalette) <- levels(sockE1_pca.vectors$Runtime)

### SOCKEYE PLOT
sock_er1 <- ggplot(data = sockE1_pca.vectors, 
                   aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = 3, color = "gray20") +
  theme(legend.position = "bottom") +
  scale_fill_manual(name = 'Run Timing', values = mypalette) +
  scale_shape_manual(name = 'Run Timing', values = c(21,22)) +
  scale_x_reverse(breaks = c(-0.1, 0,  0.1, 0.2,0.3)) +
  geom_vline(xintercept = 0.15, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.1, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle('Sockeye') +
  labs(x = paste0("PC1 (",round(sockE1_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(sockE1_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none",
        legend.text=element_text(size=32), legend.title = element_text(size=35),
        plot.title = element_text(size=20, hjust=0.5, margin=margin(0,0,10,0))) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5,
                             override.aes = list(size=8)))
# sock_er1 
#### Plotting ##########################################################

legend <- get_legend(chum_esrb)
chum_esrb <- chum_esrb + theme(legend.position = "none")

# # cowplot of leftmost 6 plot grid
# pca_six <- plot_grid(NULL,chum_lrrc9, pink_lrrc9,
#                       NULL,NULL,NULL,
#                       NULL,chum_esrb, coho_esrb,
#                       NULL,NULL,NULL,
#                       NULL,sock_er1, coho_er1,
#                       rel_heights = c(2,0.1,2,0.1,2), rel_widths = c(0.09,1,1),
#                       nrow = 5, align = 'v',
#                       labels = c('A','','','','','','B','','','C','','','','',''),
#                       label_fontfamily = "helvetica",
#                       label_size = 30, label_colour = "black")
# 
# # add sockeye lrrc9 and legend to the right
# pca_seven_legend <- plot_grid(pca_six,
#                               plot_grid(sock_lrrc9, legend, NULL, ncol=1),
#                               rel_widths=c(2, 1),
#                               align = 'hv')
# pca_seven_legend
# 
# # jpeg("./figures/pca/combine/fourspp_genes_pca_allele_cutoff_sevenPanels_legend_20260310.jpg",
# #      width = 14, height = 9, res = 300, units = "in")
# jpeg("../Salmon_runtiming/2024_fourspecies/figures/pca/combine/fourspp_genes_pca_allele_cutoff_sevenPanels_legend_20260310.jpg",
#      width = 15, height = 13, res = 300, units = "in")
# pca_seven_legend
# dev.off()

# --- Final Version: ESR1 and ESRB first --------

# cowplot of leftmost 6 plot grid
pca_flip <- plot_grid(NULL,NULL,NULL,
                      sock_er1, coho_er1, coho_esrb,
                      NULL,NULL,NULL,
                      chum_lrrc9, pink_lrrc9,sock_lrrc9,
                      rel_heights = c(0.1,1,0.1,1), rel_widths = c(1,1,1),
                      nrow = 4, align = 'v',
                      labels = c('A) esr1','','B) esr2b',
                                 '','','',
                                 'C) lrrc9','',''),
                      label_fontfamily = "helvetica",
                      label_size = 25, label_colour = "black")

extras_flip <- plot_grid(NULL, chum_esrb, NULL,legend, rel_heights = c(0.1,1,0.1,1), ncol=1)

# add sockeye lrrc9 and legend to the right
pca_three_genes_flip <- plot_grid(pca_flip, extras_flip,
                                  rel_widths=c(3, 1),
                                  align = 'hv')

jpeg("../Salmon_runtiming/2024_fourspecies/figures/pca/combine/fourspp_genes_pca_threegenes_allele_cutoff_20260528.jpg",
#jpeg("./figures/pca/combine/fourspp_genes_pca_allele_cutoff_sevenPanels_legend_reconfigure3_20260528.jpg",
     width = 20, height = 10, res = 500, units = "in")
pca_three_genes_flip
dev.off()

### END OF SCRIPT ############################################################
#
#
# No longer including the FST plots in the same figure as PCA?
#
# 
#
##### GFF NCBI Data ###########################################################

# find exons from gff file for genes of interest (from NCBI chum reference genome)
gff_df <- read.delim('./data/R/genomic.gff', header = F, comment.char = "#")
gff_df <- gff_df[,c(1:5,9)] # remove excess columns
colnames(gff_df) <- c("chrName", "RefSeq","exon","start.pos","fin.pos", "ID")

# only keep chr29 & 35 & 8
gff_chr35 <- gff_df %>%
  filter(chrName == "NC_068455.1")
gff_chr29 <- gff_df %>%
  filter(chrName == "NC_068449.1")
gff_chr8 <- gff_df %>%
  filter(chrName == "NC_068428.1")
rm(gff_df)

# prep pattern for str_match below
gene_pattern <- "gene=\\s*(.*?)\\s*;"           # keep string btwn "gene=" & ":product" 
exon_pattern <- "ID=exon-\\s*(.*?)\\s*;Parent"  # keep string btwn "exon=" & ";Parent" 
descr_pattern <- ";description=\\s*(.*?)\\s*;"  # keep description


##### 2A) LRRC9 FST #############################################################

# Define Boundaries
  xstart.lrrc9 = 27.86
  xend.lrrc9 = 28.24
  pca.start.lrrc9 = 28128954
  pca.end.lrrc9 = 28169980

###### Pink lrrc9 FST ###################
pink_Fst <- read.delim2("./results/fst/allele/pink-chum_NC_068455.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")
  colnames(pink_Fst) <- c("region", "chrName", "midPos", "Nsites", "Fst")

pink_Fst <- pink_Fst %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst),
         chr = 35) %>%
  select(chr, midPos, Fst)

pink_Fst$Fst[pink_Fst$Fst < 0] <- 0 # remove negative Fst values

# Cut Fst to start and end
pink_df <- pink_Fst %>%
  filter(midPos > xstart.lrrc9, midPos < xend.lrrc9) 

###### Sockeye lrrc9 FST #####################
sock_fst <- read.delim("./results/fst/allele/sock-all_NC_068455.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                       row.names = NULL,sep = "\t")
  colnames(sock_fst) <- c("region", "chrName", "midPos", "Nsites", "Fst")

# change midPos to Mb, and remove unnecessary columns
sock_fst <- sock_fst %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst),
         chr = 35) %>%
  select(chr, midPos, Fst)

sock_fst$Fst[sock_fst$Fst < 0] <- 0 # remove negative Fst values

sock_df <- sock_fst %>%
  filter(midPos > xstart.lrrc9, midPos < xend.lrrc9)

###### Chum lrrc9 FST ######################
chum_Fst <- read.delim2("./results/fst/allele/chumrun_NC_068455.1_EE-LL_minInd0.3.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")
  colnames(chum_Fst) <- c("region", "chrName", "midPos", "Nsites", "Fst")

chum_Fst <- chum_Fst %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst),
         chr = 35) %>%
  select(chr, midPos, Fst)

chum_Fst$Fst[chum_Fst$Fst < 0] <- 0 # remove negative Fst values

# filter both pops to desired start and end point
chum_df <- chum_Fst %>%
  filter(midPos > xstart.lrrc9, midPos < xend.lrrc9) 

####### Add Genes ########################################################

# only the region of interest
gff_region <- gff_chr35 %>%
  mutate(start.pos = start.pos/1e6, fin.pos = fin.pos/1e6) %>% 
  filter(fin.pos > xstart.lrrc9, start.pos < xend.lrrc9)

# create new columns for genes and exons from ID
gff_region_exon <- gff_region %>%
  filter(exon == "exon") %>%
  mutate(gene = str_match(ID, gene_pattern)[,2],  # gene abbr.      
         exonID = str_match(ID,exon_pattern)[,2]) # mRNA name and exon number

# this exon in six6a is too small that it doesn't even plot
# make slightly larger so it is visible in plot
gff_region_exon$fin.pos[which(gff_region_exon$gene == "six6a")[1]] <- 27.994500 # changed from 27.994277

exons_to_plot <- gff_region_exon[,c(4,5,7:8)] # only retain columns of interest
exons_to_plot <- exons_to_plot[!grepl("^LOC",exons_to_plot$gene),] # remove uncharacterized loci
exons_to_plot$gene <- factor(exons_to_plot$gene, levels = unique(exons_to_plot$gene)) # factor based on gene name

# The below file is manually edited for plotting purposes
# Colors are assigned to each gene
genes_df <- read.delim2("./data/R/chr35_chum_genes_exons.txt", header = T, 
                        sep = "\t", row.names = NULL) %>%
  mutate(beg.pos = as.numeric(beg.pos), end.pos = as.numeric(end.pos),
         y.min = as.numeric(y.min), y.max = as.numeric(y.max))

# set factors for plotting columns
genes_df$gene <- factor(genes_df$gene, levels = genes_df$gene)
  mypalette <- genes_df$color
  names(mypalette) <- levels(genes_df$gene)

###### Make genes with FST < 0.5 gray w/o legend ####################
highfst_lrrc9 <- rbind(pink_df, sock_df, chum_df) %>%
  filter(Fst > 0.5) %>%
  distinct(midPos)

lowfst_lrrc9_genes <- gff_region %>%
  filter(exon == 'gene') %>%
  mutate(gene = str_match(ID, gene_pattern)[,2]) %>%  # gene abbr.      
  rowwise() %>%
  mutate(highfst = any(highfst_lrrc9$midPos >= start.pos & highfst_lrrc9$midPos <= fin.pos)) %>%
  ungroup() %>%
  filter(highfst == F)

highfst_genes <- genes_df %>%
  filter(!(gene %in% lowfst_lrrc9_genes$gene))
  
lowfst_genes <- genes_df %>%
  filter((gene %in% lowfst_lrrc9_genes$gene))

# only keep exons from genes that have color codes
highfst_exons <- filter(exons_to_plot, gene %in% highfst_genes$gene)
lowfst_exons <- filter(exons_to_plot, gene %in% lowfst_genes$gene)

######## Plotting ###########################
# Set the general themes
theme_set(
  theme( 
    panel.grid.major = element_line(color = "gray90"),  panel.grid.minor.x = element_blank(),
    axis.text.y = element_text(angle = 0, size = 18, color = "black", vjust = 0.5),
    axis.title.y = element_text(size = 22, angle = 90,
                                margin = margin(t = 0, r = 8, b = 0, l = 0)),
    strip.text.y = element_text(angle = 0), axis.line = element_line(),
    panel.background = element_rect(fill = "white"),
    panel.border = element_rect(color = "black", fill = "NA")
  )
)

pink_plot <- ggplot() +
  geom_point(data = pink_df, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") + #8C510A
  labs(y = "Pink", x = "Chromosome Position (Mb)") +
  geom_vline(xintercept = pca.start.lrrc9/1e6, linetype = "dashed", color = "navyblue") +
  geom_vline(xintercept = pca.end.lrrc9/1e6, linetype = "dashed", color = "navyblue") +
  scale_y_continuous(limits = c(-0.02, 1.02),
                     breaks = seq(0, 1, by = 0.5),
                     expand = expansion(mult = c(0.001, 0.01))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.lrrc9, xend.lrrc9),
                     breaks = seq(0, 100, by = 0.1)) +  
  theme(axis.text.x = element_text(angle = 0, size = 18, color = "black", vjust = 0.5),
        axis.title.x = element_text(size = 20), strip.text.y = element_text(angle = 0),
        plot.margin = unit(c(0.1,0.15,0.1,0.05), "cm"))

sock_plot <- ggplot() +
  geom_point(data = sock_df, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") + 
  geom_vline(xintercept = pca.start.lrrc9/1e6, linetype = "dashed", color = "navyblue") +
  geom_vline(xintercept = pca.end.lrrc9/1e6, linetype = "dashed", color = "navyblue") +
  labs(y = "Sockeye") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(-0.02, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.lrrc9, xend.lrrc9),
                     breaks = seq(0, 100, by = 0.1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        plot.margin = unit(c(0.1,0.15,0.1,0.05), "cm"))

chum_plot <- ggplot() +
  geom_point(data = chum_df, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") +
  geom_vline(xintercept = pca.start.lrrc9/1e6, linetype = "dashed", color = "navyblue") +
  geom_vline(xintercept = pca.end.lrrc9/1e6, linetype = "dashed", color = "navyblue") +
  labs(x="Chromosome Position (Mb)", y = "Chum") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(-0.01, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.lrrc9, xend.lrrc9),
                     breaks = seq(0, 100, by = 0.1)) +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        strip.text.y = element_text(angle = 0),
        plot.margin = unit(c(0.1,0.15,0.1,0.05), "cm"))

####### plot genes/exons
gene_plot <- ggplot() +
  geom_rect(data = highfst_genes, aes(xmin = beg.pos, xmax = end.pos, ymin = y.min, ymax = y.max,
                                 fill = gene), alpha = 1) +
  geom_rect(data = highfst_exons, aes(xmin = start.pos, xmax = fin.pos, ymin = 0, ymax = 0.1,
                                      fill = gene)) +
  # gray blocks without legend
  geom_rect(data = lowfst_genes, aes(xmin = beg.pos, xmax = end.pos, ymin = y.min, ymax = y.max),
                                 fill = "gray70", alpha = 0.7) +
  geom_rect(data = lowfst_exons, aes(xmin = start.pos, xmax = fin.pos, ymin = 0, ymax = 0.1),
                                      fill = "gray70", alpha = 0.7) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.lrrc9, xend.lrrc9)) +
  scale_fill_manual(values = mypalette) +
  guides(fill = guide_legend(title = "Genes")) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    strip.text.y = element_blank(), panel.spacing = unit(0.1,"lines"),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    legend.text = element_text(size = 24), legend.title = element_text(size = 26),
    axis.line = element_blank(), plot.margin = unit(c(0.05,0,0,0), "cm"))
gene_plot

######### Combine Plots
multiplot_temp <- gene_plot / sock_plot / chum_plot / pink_plot + 
  plot_layout(heights = c(0.25, 1, 1, 1), guides = "collect") &
  theme(legend.position = 'right',
        legend.justification = 'top', legend.justification.right = c(0,0.8))

# add FST as separate label (it will be it's own plot)
y_lab <- ggplot() + 
  annotate(geom = "text", x = 1, y = 1, label = expression(italic(F[ST])), angle = 90, size = 10) + 
  coord_cartesian(clip = "off") +
  theme_void() 

# combine FST label to other plots 
multiplot_lrrc9 <- (y_lab - multiplot_temp) + # patchwork uses hyphen to allow for lefthand additions
  plot_layout(widths = c(1, 12))
multiplot_lrrc9


##### 2B) ESRB FST ##############################################################

# which region of chr29 to plot (first and last position)
# panel spanning larger region (Fig 5)
xstart.esrb = 24.7
xend.esrb = 26.4

# where were boundaries for allele-based PCA
pca.start.esrb = 25414060
pca.end.esrb = 25501622

###### Chum esrb FST ##############################

chum_Fst <- read.delim2("./results/fst/allele/chumrun_NC_068449.1_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")
  colnames(chum_Fst) <- c("region", "chrName", "midPos", "Nsites", "Fst")

chum_Fst <- chum_Fst %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst),
         chr = 29) %>%
  dplyr::select(chr, midPos, Fst)

chum_Fst$Fst[chum_Fst$Fst < 0] <- 0 # remove negative Fst values

# filter both pops to desired start and end point
chum_df <- chum_Fst %>%
  filter(midPos > xstart.esrb, midPos < xend.esrb)

###### Coho esrb FST ###########################

coho_Fst <- read.delim2("./results/fst/allele/coho-chum_NC_068449.1_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")
  colnames(coho_Fst) <- c("region", "chrName", "midPos", "Nsites", "Fst")

coho_Fst <- coho_Fst %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst),
         chr = 29) %>%
  dplyr::select(chr, midPos, Fst)

coho_Fst$Fst[coho_Fst$Fst < 0] <- 0 # remove negative Fst values

# filter both pops to desired start and end point
coho_df <- coho_Fst %>%
  filter(midPos > xstart.esrb, midPos < xend.esrb) 

###### Add Genes ##############################################################

# only the region of interest
gff_region <- gff_chr29 %>%
  mutate(start.pos = start.pos/1e6, fin.pos = fin.pos/1e6) %>% 
  filter(fin.pos > xstart.esrb, start.pos < xend.esrb)

# create new columns for genes and exons from ID
gff_region_exon <- gff_region %>%
  filter(exon == "exon") %>%
  mutate(gene = str_match(ID, gene_pattern)[,2],  # gene abbr.      
         exonID = str_match(ID,exon_pattern)[,2]) # mRNA name and exon number

exons_df <- gff_region_exon[,c(4,5,7:8)] # only retain columns of interest

# this file was manually edited for plotting purposes based on above description
# This filtered file removed all uncharacterized loci except those of great interest (esrb, six genes)
genes_df <- read.csv("./data/R/chr29_chum_genes_exons_filtered.csv", 
                     header = T, row.names = NULL) %>%
  mutate(beg.pos = as.numeric(beg.pos), end.pos = as.numeric(end.pos),
         y.min = as.numeric(y.min), y.max = as.numeric(y.max))

genes_df$geneAbbr <- factor(genes_df$geneAbbr, levels = genes_df$geneAbbr) # set factors for plotting columns 
exons_to_plot <- inner_join(exons_df, genes_df, by = "gene")
exons_to_plot$geneAbbr <- factor(exons_to_plot$geneAbbr, levels = unique(exons_to_plot$geneAbbr)) # factor based on gene name

unique(exons_to_plot$gene)

######### edit exons to plot - some are so small they don't register
exons_to_plot <- exons_to_plot %>% 
  mutate(beg.pos = ifelse(beg.pos < xstart.esrb, xstart.esrb, beg.pos),
         end.pos = ifelse(end.pos > xend.esrb, xend.esrb, end.pos),
         exon_length = (fin.pos - start.pos)*1e6,
         fin.pos = ifelse(exon_length > 1200, fin.pos, fin.pos + 0.5*(1200 - exon_length)/1e6),
         beg.pos = ifelse(exon_length > 1200, beg.pos, beg.pos - 0.5*(1200 - exon_length)/1e6),
         # run again
         fin.pos = ifelse(exon_length > 1200, fin.pos, fin.pos + 0.5*(1200 - exon_length)/1e6),
         beg.pos = ifelse(exon_length > 1200, beg.pos, beg.pos - 0.5*(1200 - exon_length)/1e6))

# this exon in esrb is too small even with the above edit
# make slightly larger so it is visible in plot
exons_to_plot$fin.pos[which(exons_to_plot$geneAbbr == "esrb")[length(which(exons_to_plot$geneAbbr == "esrb"))]] <- 25.456 # changed from 27.994277

###### Make genes with FST < 0.5 gray w/o legend ####################
highfst_esrb <- rbind(coho_df, chum_df) %>%
  filter(Fst > 0.5) %>%
  distinct(midPos)

lowfst_esrb_genes <- gff_region %>%
  filter(exon == 'gene') %>%
  mutate(gene = str_match(ID, gene_pattern)[,2]) %>%  # gene abbr.      
  rowwise() %>%
  mutate(highfst = any(highfst_esrb$midPos >= start.pos & highfst_esrb$midPos <= fin.pos)) %>%
  ungroup() %>%
  filter(highfst == F)

highfst_genes <- genes_df %>%
  filter(!(gene %in% lowfst_esrb_genes$gene))

lowfst_genes <- genes_df %>%
  filter((gene %in% lowfst_esrb_genes$gene))

# only keep exons from genes that have color codes
highfst_exons <- filter(exons_to_plot, gene %in% highfst_genes$gene)
lowfst_exons <- filter(exons_to_plot, gene %in% lowfst_genes$gene)

mypalette <- genes_df$color # color based on color column
  names(mypalette) <- levels(genes_df$geneAbbr)

####### Plotting #################
# Set the general themes
theme_set(
  theme( 
    axis.text.x = element_text(angle = 0, size = 18, color = "black", vjust = 0.5),
    axis.title.x = element_text(angle = 0, size = 20, color = "black"),
    axis.text.y = element_text(angle = 0, size = 18, color = "black", vjust = 0.5),
    axis.title.y = element_text(size = 22, angle = 90,
                                margin = margin(t = 0, r = 8, b = 0, l = 0)),
    strip.text.y = element_text(angle = 0), panel.grid.major = element_line(color = "gray90"),
    axis.line = element_line(), panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = "NA"),
    legend.position = "none", panel.background = element_rect(fill = "white")
  )
)

coho_plot <- ggplot() +
  geom_point(data = coho_df, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") + 
  geom_vline(xintercept = pca.start.esrb/1e6, linetype = "dashed", color = "navyblue") +
  geom_vline(xintercept = pca.end.esrb/1e6, linetype = "dashed", color = "navyblue") +
  labs(y = "Coho") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(-0.01, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.esrb, xend.esrb),
                     breaks = seq(0, 100, by = 0.4)) +
  theme(
    plot.margin = unit(c(0.1,0.15,0.05,0.05), "cm"),
    axis.title.x = element_blank(), axis.text.x = element_blank())

chum_esrb_plot <- ggplot() +
  geom_point(data = chum_df, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") + 
  geom_vline(xintercept = pca.start.esrb/1e6, linetype = "dashed", color = "navyblue") +
  geom_vline(xintercept = pca.end.esrb/1e6, linetype = "dashed", color = "navyblue") +
  labs(x="Chromosome Position (Mb)", y = "Chum") +
  scale_y_continuous(limits = c(-0.02, 1.02),
                     breaks = seq(0, 1, by = 0.5),
                     expand = expansion(mult = c(0.001, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.esrb, xend.esrb),
                     breaks = seq(0, 100, by = 0.4)) +          
  theme(plot.margin = unit(c(0.15,0.05,0.1,0.05), "cm"))

####### plot genes/exons
gene_plot <- ggplot() +
  geom_rect(data = highfst_exons, aes(xmin = beg.pos, xmax = end.pos, ymin = y.min, ymax = y.max,
                                      fill = geneAbbr)) +
  geom_rect(data = highfst_exons, aes(xmin = start.pos, xmax = fin.pos, ymin = 0, ymax = 0.1,
                                      fill = geneAbbr)) +
  geom_rect(data = lowfst_exons, aes(xmin = beg.pos, xmax = end.pos, ymin = y.min, ymax = y.max),
                                     fill = "gray70", alpha = 0.7) +
  geom_rect(data = lowfst_exons, aes(xmin = start.pos, xmax = fin.pos, ymin = 0, ymax = 0.1),
                                     fill = "gray70", alpha = 0.7) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.esrb, xend.esrb)) +
  scale_fill_manual(values = mypalette) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    strip.text.y = element_blank(), axis.line = element_blank(),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.spacing = unit(0.1,"lines"),
    legend.text = element_text(size = 24), legend.title = element_text(size = 26),
    legend.position = "right", plot.margin = unit(c(0.05,0,0,0), "cm")) +
  guides(fill = guide_legend(title = "Genes")) 

# plot three figures on top of one another
multiplot_temp <- gene_plot / coho_plot / chum_esrb_plot + 
  plot_layout(heights = c(0.2, 1, 1),
              guides = "collect") & 
  theme(legend.position = "right", legend.justification = 'top')

# combine FST label to other plots 
multiplot_esrb <- (y_lab - multiplot_temp) + # patchwork uses hyphen to allow for lefthand additions
  plot_layout(widths = c(1, 12))
multiplot_esrb


##### 2C) ESR1 FST ##############################################################

# which region of chr29 to plot (first and last position)
# panel spanning larger region (Fig 5)
xstart.er1 = 24.7
xend.er1 = 26.4

# where were boundaries for allele-based PCA
pca.start.er1 = 25414060
pca.end.er1 = 25501622

###### Sockeye er1 FST ##############################

sock_Fst <- read.delim2("./results/fst/allele/chumrun_NC_068449.1_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")
colnames(sock_Fst) <- c("region", "chrName", "midPos", "Nsites", "Fst")

sock_Fst <- sock_Fst %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst),
         chr = 29) %>%
  dplyr::select(chr, midPos, Fst)

sock_Fst$Fst[sock_Fst$Fst < 0] <- 0 # remove negative Fst values

# filter both pops to desired start and end point
sock_df <- sock_Fst %>%
  filter(midPos > xstart.er1, midPos < xend.er1)

###### Coho er1 FST ###########################

coho_Fst <- read.delim2("./results/fst/allele/coho-chum_NC_068449.1_EE-LL_minInd0.3_minDepthHalf.sfs.pbs.fst.txt",
                        row.names = NULL,sep = "\t")
colnames(coho_Fst) <- c("region", "chrName", "midPos", "Nsites", "Fst")

coho_Fst <- coho_Fst %>%
  mutate(midPos = as.numeric(midPos)/1e6,
         Fst = as.numeric(Fst),
         chr = 29) %>%
  dplyr::select(chr, midPos, Fst)

coho_Fst$Fst[coho_Fst$Fst < 0] <- 0 # remove negative Fst values

# filter both pops to desired start and end point
coho_df <- coho_Fst %>%
  filter(midPos > xstart.er1, midPos < xend.er1) 

###### Add Genes ##############################################################

# only the region of interest
gff_region <- gff_chr29 %>%
  mutate(start.pos = start.pos/1e6, fin.pos = fin.pos/1e6) %>% 
  filter(fin.pos > xstart.er1, start.pos < xend.er1)

# create new columns for genes and exons from ID
gff_region_exon <- gff_region %>%
  filter(exon == "exon") %>%
  mutate(gene = str_match(ID, gene_pattern)[,2],  # gene abbr.      
         exonID = str_match(ID,exon_pattern)[,2]) # mRNA name and exon number

exons_df <- gff_region_exon[,c(4,5,7:8)] # only retain columns of interest

# this file was manually edited for plotting purposes based on above description
# This filtered file removed all uncharacterized loci except those of great interest (er1, six genes)
genes_df <- read.csv("./data/R/chr29_chum_genes_exons_filtered.csv", 
                     header = T, row.names = NULL) %>%
  mutate(beg.pos = as.numeric(beg.pos), end.pos = as.numeric(end.pos),
         y.min = as.numeric(y.min), y.max = as.numeric(y.max))

genes_df$geneAbbr <- factor(genes_df$geneAbbr, levels = genes_df$geneAbbr) # set factors for plotting columns 
exons_to_plot <- inner_join(exons_df, genes_df, by = "gene")
exons_to_plot$geneAbbr <- factor(exons_to_plot$geneAbbr, levels = unique(exons_to_plot$geneAbbr)) # factor based on gene name

unique(exons_to_plot$gene)

######### edit exons to plot - some are so small they don't register
exons_to_plot <- exons_to_plot %>% 
  mutate(beg.pos = ifelse(beg.pos < xstart.er1, xstart.er1, beg.pos),
         end.pos = ifelse(end.pos > xend.er1, xend.er1, end.pos),
         exon_length = (fin.pos - start.pos)*1e6,
         fin.pos = ifelse(exon_length > 1200, fin.pos, fin.pos + 0.5*(1200 - exon_length)/1e6),
         beg.pos = ifelse(exon_length > 1200, beg.pos, beg.pos - 0.5*(1200 - exon_length)/1e6),
         # run again
         fin.pos = ifelse(exon_length > 1200, fin.pos, fin.pos + 0.5*(1200 - exon_length)/1e6),
         beg.pos = ifelse(exon_length > 1200, beg.pos, beg.pos - 0.5*(1200 - exon_length)/1e6))

# this exon in er1 is too small even with the above edit
# make slightly larger so it is visible in plot
exons_to_plot$fin.pos[which(exons_to_plot$geneAbbr == "er1")[length(which(exons_to_plot$geneAbbr == "er1"))]] <- 25.456 # changed from 27.994277

###### Make genes with FST < 0.5 gray w/o legend ####################
highfst_er1 <- rbind(coho_df, sock_df) %>%
  filter(Fst > 0.5) %>%
  distinct(midPos)

lowfst_er1_genes <- gff_region %>%
  filter(exon == 'gene') %>%
  mutate(gene = str_match(ID, gene_pattern)[,2]) %>%  # gene abbr.      
  rowwise() %>%
  mutate(highfst = any(highfst_er1$midPos >= start.pos & highfst_er1$midPos <= fin.pos)) %>%
  ungroup() %>%
  filter(highfst == F)

highfst_genes <- genes_df %>%
  filter(!(gene %in% lowfst_er1_genes$gene))

lowfst_genes <- genes_df %>%
  filter((gene %in% lowfst_er1_genes$gene))

# only keep exons from genes that have color codes
highfst_exons <- filter(exons_to_plot, gene %in% highfst_genes$gene)
lowfst_exons <- filter(exons_to_plot, gene %in% lowfst_genes$gene)

mypalette <- genes_df$color # color based on color column
names(mypalette) <- levels(genes_df$geneAbbr)

####### Plotting #################
# Set the general themes
theme_set(
  theme( 
    axis.text.x = element_text(angle = 0, size = 18, color = "black", vjust = 0.5),
    axis.title.x = element_text(angle = 0, size = 20, color = "black"),
    axis.text.y = element_text(angle = 0, size = 18, color = "black", vjust = 0.5),
    axis.title.y = element_text(size = 22, angle = 90,
                                margin = margin(t = 0, r = 8, b = 0, l = 0)),
    strip.text.y = element_text(angle = 0), panel.grid.major = element_line(color = "gray90"),
    axis.line = element_line(), panel.grid.minor.x = element_blank(),
    panel.border = element_rect(color = "black", fill = "NA"),
    legend.position = "none", panel.background = element_rect(fill = "white")
  )
)

coho_plot <- ggplot() +
  geom_point(data = coho_df, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") + 
  geom_vline(xintercept = pca.start.er1/1e6, linetype = "dashed", color = "navyblue") +
  geom_vline(xintercept = pca.end.er1/1e6, linetype = "dashed", color = "navyblue") +
  labs(y = "Coho") +
  scale_y_continuous(breaks = seq(0, 1, by = 0.5),
                     limits = c(-0.01, 1.02),
                     expand = expansion(mult = c(0.01, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.er1, xend.er1),
                     breaks = seq(0, 100, by = 0.4)) +
  theme(
    plot.margin = unit(c(0.1,0.15,0.05,0.05), "cm"),
    axis.title.x = element_blank(), axis.text.x = element_blank())

chum_er1_plot <- ggplot() +
  geom_point(data = sock_df, aes(x = midPos, y = Fst), 
             size = 2, alpha = 0.6, color = "gray10") + 
  geom_vline(xintercept = pca.start.er1/1e6, linetype = "dashed", color = "navyblue") +
  geom_vline(xintercept = pca.end.er1/1e6, linetype = "dashed", color = "navyblue") +
  labs(x="Chromosome Position (Mb)", y = "Chum") +
  scale_y_continuous(limits = c(-0.02, 1.02),
                     breaks = seq(0, 1, by = 0.5),
                     expand = expansion(mult = c(0.001, 0.001))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.er1, xend.er1),
                     breaks = seq(0, 100, by = 0.4)) +          
  theme(plot.margin = unit(c(0.15,0.05,0.1,0.05), "cm"))

####### plot genes/exons
gene_plot <- ggplot() +
  geom_rect(data = highfst_exons, aes(xmin = beg.pos, xmax = end.pos, ymin = y.min, ymax = y.max,
                                      fill = geneAbbr)) +
  geom_rect(data = highfst_exons, aes(xmin = start.pos, xmax = fin.pos, ymin = 0, ymax = 0.1,
                                      fill = geneAbbr)) +
  geom_rect(data = lowfst_exons, aes(xmin = beg.pos, xmax = end.pos, ymin = y.min, ymax = y.max),
            fill = "gray70", alpha = 0.7) +
  geom_rect(data = lowfst_exons, aes(xmin = start.pos, xmax = fin.pos, ymin = 0, ymax = 0.1),
            fill = "gray70", alpha = 0.7) +
  theme_minimal() +
  scale_y_continuous(limits = c(0, 0.1),
                     expand = expansion(mult = c(0, 0))) +
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.01)),
                     limits = c(xstart.er1, xend.er1)) +
  scale_fill_manual(values = mypalette) +
  theme(
    axis.text = element_blank(), axis.ticks = element_blank(),
    axis.title.y = element_blank(), axis.title.x = element_blank(),
    strip.text.y = element_blank(), axis.line = element_blank(),
    panel.background = element_blank(), panel.border = element_blank(),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.spacing = unit(0.1,"lines"),
    legend.text = element_text(size = 24), legend.title = element_text(size = 26),
    legend.position = "right", plot.margin = unit(c(0.05,0,0,0), "cm")) +
  guides(fill = guide_legend(title = "Genes")) 

# plot three figures on top of one another
multiplot_temp <- gene_plot / coho_plot / chum_er1_plot + 
  plot_layout(heights = c(0.2, 1, 1),
              guides = "collect") & 
  theme(legend.position = "right", legend.justification = 'top')

# combine FST label to other plots 
multiplot_er1 <- (y_lab - multiplot_temp) + # patchwork uses hyphen to allow for lefthand additions
  plot_layout(widths = c(1, 12))
multiplot_er1


##### C/D) Combine #####################
twofst_cowplot <- plot_grid(multiplot_lrrc9, multiplot_esrb,
                            rel_widths = c(5,5), ncol = 2, nrow = 1,
                            labels = c('C','D'), label_fontfamily = "helvetica",
                            label_size = 30, label_colour = "black")
twofst_cowplot

##### FIG 3) Combine A-D #####################

fig3 <- plot_grid(pca_five_legend, NULL, twofst_cowplot,
                  rel_heights = c(1,0.1, 1), 
                  ncol = 1, nrow = 3, align = 'v')
fig3

jpeg(paste0("./figures/figure3_panelsA-D_",format(Sys.Date(),"%Y%m%d"),".jpg"), 
     width = 18, height = 20, res = 200, units = "in")
print(fig3)
dev.off()
