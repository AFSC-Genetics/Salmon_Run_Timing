# COMBINE PCAs and FST for LRRC9, ESRB, and ESR1 
# FIGURE 3:
#  A) ESR1 COMBINE PCAS
#  B) ESRB COMBINE PCAS
#  C) LRRC9 COMBINE PCAS 

packages_needed <- c("ggplot2", "scales", "ggpubr", "ggrepel", "stringr", 
                     "data.table", "plyr","tools","gtools","reshape2", 
                     "patchwork", "cowplot", "tidyverse")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}
rm(i, packages_needed)

# metadata
allmeta_df <- read.csv("./data/raw/fourspecies_runtiming_metadata.csv", header = T)

#### 1C) LRRC9 ---------------------------------------------------------------

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
pink_pca.vectors$Runtime <- factor(pink_pca.vectors$Runtime, levels = c("Early", "Late"))

##### Sockeye - Euclide ####################

# read in the covariance matrix
sock_cov <- as.matrix(read.table("./results/pca/euclide_NC_068455.1_minInd0.3_lrrc9.cov"))
sock_e <- eigen(sock_cov) # calculate eigenvector values
sock_e_vectors <- as.data.frame(sock_e$vectors)
sock_varPC1 <- (sock_e$values[1]/sum(sock_e$values))*100 # PC1 variance
sock_varPC2 <- (sock_e$values[2]/sum(sock_e$values))*100 # PC2 variance

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

sock_pca.vectors = as_tibble(cbind(sock_meta, sock_e_vectors))
sock_pca.vectors$Runtime <- factor(sock_pca.vectors$Runtime, levels = c("Early", "Late"))

##### Chum LRRC9 ###

# read in the covariance matrix
chum_cov <- as.matrix(read.table("./results/pca/chumrun_NC_068455.1_lrrc9_minInd0.3_minDepthHalf.cov"))
chum_e <- eigen(chum_cov) # calculate eigenvector values
chum_e_vectors <- as.data.frame(chum_e$vectors)
chum_varPC1 <- (chum_e$values[1]/sum(chum_e$values))*100 # PC1 variance
chum_varPC2 <- (chum_e$values[2]/sum(chum_e$values))*100 # PC2 variance

chum_bam_df <- read.table("./data/bams/chumrun_bamslist.txt", header = F) %>%
  mutate(ABLG = as.numeric(gsub('[^0-9]','', V1))) %>% select(ABLG)

chum_meta <- allmeta_df %>%
  filter(Species == 'Chum') %>%
  mutate(ABLG = as.numeric(sub('ABLG','',sampleID))) %>%
  inner_join(chum_bam_df, ., by = "ABLG")

chum_pca.vectors = as_tibble(cbind(chum_meta, chum_e_vectors))
chum_pca.vectors$Runtime <- factor(chum_pca.vectors$Runtime, levels = c("Early", "Late"))

#### 1B) ESRB -----------------------------------------------------------------

###### COHO ###########

# read in the covariance matrix
coho_cov <- as.matrix(read.table("./results/pca/coho-chum_NC_068449.1_s25414060_e25501622_esrb_minInd0.3_minDepthHalf.cov"))
coho_e <- eigen(coho_cov)
coho_e_vectors <- as.data.frame(coho_e$vectors)
coho_varPC1 <- (coho_e$values[1]/sum(coho_e$values))*100 # PC1 variance
coho_varPC2 <- (coho_e$values[2]/sum(coho_e$values))*100 # PC2 variance

coho_bam_df <- read.table("./data/bams/coho-chum_filtered_bamslist.txt", header = F) %>%
  mutate(ABLG = as.numeric(gsub("[^0-9]", "", V1))) %>% select(-V1)

coho_meta <- allmeta_df %>%
  filter(Species == 'Coho') %>%
  mutate(ABLG = as.numeric(gsub('ABLG','',sampleID))) %>%
  inner_join(coho_bam_df, ., by = "ABLG")

coho_pca.vectors = as_tibble(cbind(coho_meta, coho_e_vectors))
coho_pca.vectors$Runtime <- factor(coho_pca.vectors$Runtime, levels = c("Early", "Late"))

##### CHUM ####################

# read in the covariance matrix
chumE_cov <- as.matrix(read.table("./results/pca/chumrun_NC_068449.1_s25414060_e25501622_esrb_minInd0.3_minDepthHalf.cov"))
chumE_e <- eigen(chumE_cov)
chumE_e_vectors <- as.data.frame(chumE_e$vectors)
chumE_varPC1 <- (chumE_e$values[1]/sum(chumE_e$values))*100 # PC1 variance
chumE_varPC2 <- (chumE_e$values[2]/sum(chumE_e$values))*100 # PC2 variance

chumE_pca.vectors = as_tibble(cbind(chum_meta, chumE_e_vectors))
chumE_pca.vectors$Runtime <- factor(chumE_pca.vectors$Runtime, levels = c("Early", "Late"))

#### 1A) ESR1 ----------------------------------------------------------------

###### COHO ###########

# read in the covariance matrix
cohoE1_cov <- as.matrix(read.table(paste0("./results/pca/coho-chum_NC_068428.1_s23264160_e23645630_minInd0.3_minDepthHalf_esr1.cov")))
cohoE1_e <- eigen(cohoE1_cov)
cohoE1_e_vectors <- as.data.frame(cohoE1_e$vectors)
cohoE1_varPC1 <- (cohoE1_e$values[1]/sum(cohoE1_e$values))*100 # PC1 variance
cohoE1_varPC2 <- (cohoE1_e$values[2]/sum(cohoE1_e$values))*100 # PC2 variance

cohoE1_pca.vectors = as_tibble(cbind(coho_meta, cohoE1_e_vectors))
cohoE1_pca.vectors$Runtime <- factor(cohoE1_pca.vectors$Runtime, levels = c("Early", "Late"))

##### SOCKEYE EXPAND ####################

# read in the covariance matrix
sockE1_cov <- as.matrix(read.table(paste0("./results/pca/euclide_NC_068428.1_s23264160_e23645630_minInd0.3_esr1.cov")))
sockE1_e <- eigen(sockE1_cov)
sockE1_e_vectors <- as.data.frame(sockE1_e$vectors)
sockE1_varPC1 <- (sockE1_e$values[1]/sum(sockE1_e$values))*100 # PC1 variance
sockE1_varPC2 <- (sockE1_e$values[2]/sum(sockE1_e$values))*100 # PC2 variance

sockE1_pca.vectors = as_tibble(cbind(sock_meta, sockE1_e_vectors))
sockE1_pca.vectors$Runtime <- factor(sockE1_pca.vectors$Runtime, levels = c("Early", "Late"))

#### Plotting --------------------------------------------------------- ########

# color blind palette colors
cb_palette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
                "#0072B2", "#D55E00", "#CC79A7", "#999999")

# Prep Plot Formatting
mypalette <- c('Early' =  "#F0E442",
               'Late'  =  "#0072B2")

# for science, consistent font sizing
width2col <- 7.25  # 2 col, in inches
width1col <- 3.6  # 1 col, in inches

largeFS <- 10 # Panel Labels
mediumFS <- 8  # Axis titles, legend titles
smallFS <- 6  # Axis text, annotations, captions

theme_set(
  theme(legend.text=element_text(size=mediumFS),
        legend.title = element_text(size = largeFS),
        axis.text = element_text(size = smallFS, angle = 0),
        axis.title.x = element_text(size = smallFS),
        axis.title.y = element_text(size = smallFS, angle=90, margin=margin(0,0,0,0)),
        plot.title = element_text(size = mediumFS, hjust = 0.5, margin=margin(0,0,2,0)),
        panel.background = element_rect(fill = "white"), 
        panel.spacing = unit(0,"lines"),
        strip.text.y = element_text(angle = 0),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()
  )
)

pointSize = 1.2

## --- LRRC9 ------------------------------------------------------------ ##

# PINK LRRC9
pink_lrrc9 <- ggplot(data = pink_pca.vectors, 
                     aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = pointSize, color = "gray20") +
  scale_fill_manual(name = expression('Pink'~italic(lrrc9)), values = mypalette) +
  scale_shape_manual(name = expression('Pink'~italic(lrrc9)), values = c(21,22)) +
  scale_x_continuous(breaks = c(-0.05, 0, 0.05, 0.1, 0.15)) +
  geom_vline(xintercept = 0.07, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle("Pink") + 
  # ggtitle(expression('Pink'~italic(lrrc9))) +
  labs(x = paste0("PC1 (",round(pink_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(pink_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none")

# SOCKEYE LRRC9
sock_lrrc9 <- ggplot(data = sock_pca.vectors, 
                     aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = pointSize, color = "gray20") +
  scale_fill_manual(name = expression('Sockeye'~italic(lrrc9)), values = mypalette) +
  scale_shape_manual(name = expression('Sockeye'~italic(lrrc9)), values = c(21,22)) +
  geom_vline(xintercept = 0.1, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle("Sockeye") + 
  # ggtitle(expression('Sockeye'~italic(lrrc9))) +
  scale_x_reverse(breaks = c(-0.1, -0.05, 0, 0.05, 0.1)) +
  labs(x = paste0("PC1 (",round(sock_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(sock_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none")

# CHUM LRRC9
chum_lrrc9 <- ggplot(data = chum_pca.vectors, 
                     aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = pointSize, color = "gray20") +
  scale_fill_manual(name = expression('Chum'~italic(lrrc9)), values = mypalette) +
  scale_shape_manual(name = expression('Chum'~italic(lrrc9)), values = c(21,22)) +
  geom_vline(xintercept = 0.1, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle('Chum') +
  # ggtitle(expression('Chum'~italic(lrrc9))) +
  scale_x_reverse(breaks = c(-0.05, 0, 0.05, 0.1, 0.15)) +
  labs(x = paste0("PC1 (",round(chum_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(chum_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none")

## --- ESR2B ---------------------------------------------------------------------##

# COHO ESR2B
coho_esrb <- ggplot(data = coho_pca.vectors, 
                    aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = pointSize, color = "gray20") +
  scale_fill_manual(name = expression('Coho'~italic(esrb)), values = mypalette) +
  scale_shape_manual(name = expression('Coho'~italic(esrb)), values = c(21,22)) +
  scale_x_continuous(breaks = c(-.1,-0.05, 0, 0.05, 0.1, 0.15)) +
  geom_vline(xintercept = 0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle('Coho') +
  # ggtitle(expression('Coho'~italic(esr2b))) +
  labs(x = paste0("PC1 (",round(coho_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(coho_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none")

# CHUM ESR2B
chum_esrb <- ggplot(data = chumE_pca.vectors, 
                    aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = pointSize, color = "gray20") +
  scale_fill_manual(name = 'Run Timing', values = mypalette) +
  scale_shape_manual(name = 'Run Timing', values = c(21,22)) +
  scale_x_continuous(breaks = c(-0.1,-0.05, 0, 0.05, 0.1, 0.15)) +
  geom_vline(xintercept = 0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle('Chum') +
  # ggtitle(expression('Chum'~italic(esr2b))) +
  labs(x = paste0("PC1 (",round(chumE_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(chumE_varPC2, digits = 1),"%)")) +
  theme(legend.position = "right") +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5,
                             override.aes = list(size=smallFS)))

## --- ESR1 --------------------------------------------------------------- ##

# COHO ESR1
coho_er1 <- ggplot(data = cohoE1_pca.vectors, 
                   aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = pointSize, color = "gray20") +
  scale_fill_manual(name = expression('Coho'~italic(er1)), values = mypalette) +
  scale_shape_manual(name = expression('Coho'~italic(er1)), values = c(21,22)) +
  scale_x_continuous(breaks = c(-.1,-0.05, 0, 0.05, 0.1)) +
  geom_vline(xintercept = 0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle('Coho') +
  # ggtitle(expression('Coho'~italic(esr1))) +
  labs(x = paste0("PC1 (",round(cohoE1_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(cohoE1_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none")

# SOCKEYE
sock_er1 <- ggplot(data = sockE1_pca.vectors, 
                   aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = pointSize, color = "gray20") +
  scale_fill_manual(name = 'Run Timing', values = mypalette) +
  scale_shape_manual(name = 'Run Timing', values = c(21,22)) +
  scale_x_reverse(breaks = c(-0.1, 0,  0.1, 0.2,0.3)) +
  geom_vline(xintercept = 0.15, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.1, color = "gray30", alpha = 0.5, linetype = "dashed") +
  ggtitle('Sockeye') +
  # ggtitle(expression('Sockeye'~italic(esr1))) +
  labs(x = paste0("PC1 (",round(sockE1_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(sockE1_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none") +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5,
                             override.aes = list(size=5)))

## -- COMBINE ------------------------------------------------------------ ##
legend <- get_legend(chum_esrb)
chum_esrb <- chum_esrb + theme(legend.position = "none")

esr1_title <- ggdraw() + 
  draw_label(expression(italic(esr1)),
             x = 0, hjust = 0, vjust=0.75,
             fontfamily = "sans",
             size = largeFS) + 
  theme(plot.margin = margin(0, 0, 0, 7))

esrb_title <- ggdraw() + 
  draw_label(expression(italic(esr2b)),
             x = 0, hjust = 0, vjust=0.75,
             fontfamily = "sans",
             size = largeFS) + 
  theme(plot.margin = margin(0, 0, 0, 7))

lrrc9_title <- ggdraw() + 
  draw_label(expression(italic(lrrc9)),
             x = 0, hjust = 0, vjust=0.75,
             fontfamily = "sans",
             size = largeFS) + 
  theme(plot.margin = margin(0, 0, 0, 7))

# --- Version 1: Longer - Row by Gene --------
# 2 cols
pca_long <- plot_grid(esr1_title,   NULL,
                      sock_er1,     coho_er1,
                      esrb_title,   NULL,
                      chum_esrb,    coho_esrb, 
                      lrrc9_title,  NULL,
                      chum_lrrc9,   sock_lrrc9,
                      pink_lrrc9,   legend,
                      rel_heights = c(0.1,1,0.1,1,0.1,1,1), rel_widths = c(1,1),
                      nrow = 7, align = 'v', axis= "lr",
                      labels = c('A','',
                                 '','',
                                 'B','',
                                 '','',
                                 'C','',
                                 '',''),
                      label_fontfamily = "sans",
                      label_size = largeFS, label_colour = "black")

# for import into MS
png(str_c("./figures/pca/Figure3_Allele_PCAs_Long_PanelGeneName_",format(Sys.Date(),"%Y%m%d"),".png"),
    width= width1col, height= 1.9*width1col, res = 500, units = "in")
pca_long
dev.off()

# svg for science 
svglite(str_c("./figures/pca/Figure3_Allele_PCAs_Long_PanelGeneName_",format(Sys.Date(),"%Y%m%d"),".svg"),
        width= width1col, height= 1.9*width1col)
print(pca_long)
dev.off()

# --- Version 2: Wider - Estrogen genes on top, LRRC9 on bottom --------

# 2 rows
pca_flip <- plot_grid(esr1_title,NULL,esrb_title,
                      sock_er1, coho_er1, coho_esrb,
                      lrrc9_title,NULL,NULL,
                      chum_lrrc9, pink_lrrc9,sock_lrrc9,
                      rel_heights = c(0.1,1,0.1,1), rel_widths = c(1,1,1),
                      nrow = 4, align = 'v',
                      labels = c('A','','B',
                                 '','','',
                                 'C','',''),
                      label_fontfamily = "sans",
                      label_size = largeFS, label_colour = "black")

# rightmost plot separately (difficult otherwise)
extras_flip <- plot_grid(NULL, chum_esrb, NULL,legend, rel_heights = c(0.1,1,0.1,1), ncol=1)

# combine left and right
pca_three_genes_flip <- plot_grid(pca_flip, extras_flip,
                                  rel_widths=c(3, 1),
                                  align = 'hv')

# for import into MS
png(str_c("./figures/pca/Figure3_Allele_PCAs_Wide_PanelGeneName_",format(Sys.Date(),"%Y%m%d"),".png"),
    width= width2col, height= 0.5*width2col, res = 500, units = "in")
pca_three_genes_flip
dev.off()

# svg for science 
svglite(str_c("./figures/pca/Figure3_Allele_PCAs_Wide_PanelGeneName_",format(Sys.Date(),"%Y%m%d"),".svg"),
        width= width2col, height= 0.5*width2col)
print(pca_three_genes_flip)
dev.off()


