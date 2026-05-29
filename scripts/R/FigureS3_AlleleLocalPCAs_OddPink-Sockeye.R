# Supplemental local PCAs for Sockeye (Whitefish) and Pink (Odd)

library(stringr)
library(tools)
library(tidyverse)
library(patchwork)
library(cowplot)

# setwd("../../Salmon_Run_Timing")

### Import Metadata ###

allmeta_df <- read.csv("./data/raw/fourspecies_runtiming_metadata.csv", header = T)
pick_df <- read.csv("./data/raw/pick_creek_metadata.csv", header = T)
pink_odd_meta <- read.csv("./data/raw/pinkOdd_collection_07172024_cleaned.csv", header = T)

##### Whitefish Addition ############
# read in the covariance matrix
sockall_cov <- as.matrix(read.table("./results/pca/sock-all_NC_068455.1_minInd0.3_lrrc9.cov"))
  sockall_e <- eigen(sockall_cov) # calculate eigenvector values
  sockall_e_vectors <- as.data.frame(sockall_e$vectors)
  sockall_varPC1 <- (sockall_e$values[1]/sum(sockall_e$values))*100 #Variance explained by PC1
  sockall_varPC2 <- (sockall_e$values[2]/sum(sockall_e$values))*100 #Variance explained by PC2

sockall_bam_df <- read.table("./data/bams/all_sockeye_bamslist.txt", header = F) %>%
  dplyr::mutate(temp = basename(file_path_sans_ext(V1)),
                temp = gsub("^[^_]*_", "", temp), # remove everything after 1st underscore
                sampleID = str_extract(temp, "[^_]+")) %>%
  select(-c(V1, temp))

sockall_meta <- allmeta_df %>% filter(Species=="Sockeye")

# join those two dataframes
sockall_popFID <- inner_join(sockall_bam_df, sockall_meta, by = "sampleID") %>%
  mutate(Runtime = Runtime %>% 
           gsub("Early Stream",'Early Crk. (Teal)',.) %>% 
           gsub("Late Stream",'Late Crk. (Whitefish)',.) %>%
           gsub("Late Beach",'Late Beach (Anvil)',.))

##combine row names (population info) with the covariance matrix
sockall_pca.vectors <- as_tibble(cbind(sockall_popFID, sockall_e_vectors))

sockall_pca.vectors$Runtime <- factor(sockall_pca.vectors$Runtime, levels = c(unique(sockall_pca.vectors$Runtime)[-1],unique(sockall_pca.vectors$Runtime)[1]))
  sockall_Palette <-c("khaki", "skyblue1", "orchid")
  names(sockall_Palette) <- c(unique(sockall_pca.vectors$Runtime)[-1],unique(sockall_pca.vectors$Runtime)[1])

theme_set(
  theme( 
    legend.text=element_text(size=16), 
    legend.title = element_text(size = 18,v.just=1),
    legend.position = "bottom",
    plot.title = element_text(size = 20,hjust=0.5),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text = element_text(angle = 0, size = 14), axis.title = element_text(size = 16),
    panel.background = element_rect(fill = "white"), 
    panel.spacing = unit(0,"lines"), strip.text.y = element_text(angle = 0)
  )
)

sockall_lrrc9 <- ggplot(data = sockall_pca.vectors, 
                        aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = 3, color = "gray20") +
  scale_fill_manual(name = "Run Timing", values = sockall_Palette) +
  scale_shape_manual(name = "Run Timing", values = c(21,23,22)) +
  ggtitle(expression('Sockeye'~italic(lrrc9)))+
  geom_vline(xintercept = 0.1, color = "gray20", alpha = 0.5, linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = -0.05, color = "gray20", alpha = 0.5, linetype = "dashed", linewidth = 1) +
  scale_x_reverse(breaks = c(-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15)) +
  labs(x = paste0("PC1 (",round(sockall_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(sockall_varPC2, digits = 1),"%)")) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5,
                             title.theme = element_text(vjust = 1),
                             nrow = 3, byrow = TRUE,
                             override.aes = list(size =3)))
sockall_lrrc9

##### Odd-Pink Addition ############

# read in the covariance matrix
pinkO_cov <- as.matrix(read.table("./results/pca/pink-odd_NC_068455.1_lrrc9_minInd0.3.cov"))
  pinkO_e <- eigen(pinkO_cov)
  pinkO_e_vectors <- as.data.frame(pinkO_e$vectors)
  pinkO_varPC1 <- (pinkO_e$values[1]/sum(pinkO_e$values))*100 #Variance explained by PC1
  pinkO_varPC2 <- (pinkO_e$values[2]/sum(pinkO_e$values))*100 #Variance explained by PC2

# call in bams
pinkO_bam_df <- read.table("./data/bams/pink-odd_filtered_bamslist.txt", header = F) %>%
  mutate(V1 = basename(file_path_sans_ext(V1)),
         temp = gsub("^[^_]*_", "", V1), # remove everything after 1st underscore
         sampleID = str_extract(temp, "[^_]+")) %>%
  select(-c(V1, temp))

# call in some metadata
pinkO_meta <- pink_odd_meta %>%
  mutate(sampleID = paste0("ABLG",ABLG),
         Runtime = ifelse(Runtime == "late","Late","Early"))

# join those two dataframes
pinkO_pop_df <- inner_join(pinkO_bam_df, pinkO_meta, by = "sampleID")

##combine row names (population info) with the covariance matrix
pinkO_pca.vectors = as_tibble(cbind(pinkO_pop_df, pinkO_e_vectors))

pinkO_pca.vectors$Runtime <- factor(pinkO_pca.vectors$Runtime, levels = c("Early", "Late"))
  mypalette <- c("goldenrod1", "royalblue3")
  names(mypalette) <- levels(pinkO_pca.vectors$Runtime)

### PINK ODD
pinkOdd_lrrc9 <- ggplot(data = pinkO_pca.vectors, 
                        aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = 3, color = "gray20") +
  scale_fill_manual(name = "Run Timing", values = mypalette) +
  scale_shape_manual(name = "Run Timing", values = c(21,22)) +
  ggtitle(expression('Pink-Odd'~italic(lrrc9)))+
  geom_vline(xintercept = 0, color = "gray20", alpha = 0.5, linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = 0.15, color = "gray20", alpha = 0.5, linetype = "dashed", linewidth = 1) +
  scale_y_continuous(breaks = c(-.3, -.2,-.1, 0, 0.1)) +
  scale_x_reverse(breaks = c(-0.05, 0, 0.05, 0.1, 0.15, 0.2)) +
  labs(x = paste0("PC1 (",round(pinkO_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(pinkO_varPC2, digits = 1),"%)")) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5,
                             title.theme = element_text(vjust = 1),
                             nrow = 2, byrow = TRUE,
                             override.aes = list(size =3)))
pinkOdd_lrrc9

#### Sockeye ESR1 #########################################

er1_bam_df <- read.table("./data/bams/sock-chum_plusPick_bamslist.txt", header = F) %>%
  dplyr::mutate(temp = basename(file_path_sans_ext(V1)),
                temp = gsub("^[^_]*_", "", temp), # remove everything after 1st underscore
                sampleID = str_extract(temp, "[^_]+")) %>%
  select(sampleID)

er1_sock_meta <- allmeta_df %>% 
  filter(Species=="Sockeye") %>%
  select(-Species) %>%
  rbind(., pick_df %>% 
          dplyr::rename(sampleID = sample) %>%
          mutate(Runtime = paste(Runtime,"- Pick Creek")) %>%
          select(sampleID, Runtime)) %>%
  inner_join(er1_bam_df, ., by = "sampleID")

# read in the covariance matrix
er1.cov <- as.matrix(read.table("./results/pca/pick-all_NC_068428.1_s23264160_e23645630_minInd0.3_esr1.cov"))
er1.e <- eigen(er1.cov)
  er1_e_vectors <- as.data.frame(er1.e$vectors)
  er1.varPC1 <- (er1.e$values[1]/sum(er1.e$values))*100 # PC1 variance
  er1.varPC2 <- (er1.e$values[2]/sum(er1.e$values))*100 # PC2 variance

# combine row names (population info) with the covariance matrix
er1.pca.vectors <- er1_sock_meta %>%
  mutate(Runtime = Runtime %>% 
           gsub("Early Stream",'Early Crk. (Teal)',.) %>% 
           gsub("Late Stream",'Late Crk. (Whitefish)',.) %>%
           gsub("- ","Crk. (",.) %>% gsub("Pick Creek","Pick)",.) %>%
           gsub("Late Beach",'Late Beach (Anvil)',.)) %>% 
  bind_cols(er1_e_vectors) %>%
  as_tibble()

er1.pca.vectors$Runtime <- factor(er1.pca.vectors$Runtime, 
                                  levels = c(unique(er1.pca.vectors$Runtime)[-3],unique(er1.pca.vectors$Runtime)[3]))

  er1.palette <- c("goldenrod1", "royalblue3", "khaki", "skyblue1", "orchid")
  er1.shapes <- c(24,24,21,23,22)
  names(er1.palette) <- levels(er1.pca.vectors$Runtime)
  names(er1.shapes) <- levels(er1.pca.vectors$Runtime)

# plot pca with cutoffs
sockall_er1 <- ggplot(data = er1.pca.vectors, 
       aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = 3, color = "gray20") +
  scale_fill_manual(name = "Run Timing", values = er1.palette) +
  scale_color_manual(name = "Run Timing", values = c("black", "black", "firebrick4", "gold4", "skyblue1")) +
  scale_shape_manual(name = "Run Timing", values = er1.shapes) +
  ggtitle(expression('Sockeye'~italic(esr1)))+
  geom_vline(xintercept = 0.075, color = "gray20", alpha = 0.5, linetype = "dashed", linewidth = 1) +
  geom_vline(xintercept = -0.04, color = "gray20", alpha = 0.5, linetype = "dashed", linewidth = 1) +
  scale_x_continuous(breaks = c(-0.05, 0, 0.05, 0.1, 0.15, 0.2)) +
  labs(x = paste0("PC1 (",round(er1.varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(er1.varPC2, digits = 1),"%)")) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5,
                             title.theme = element_text(vjust = 1),
                             nrow = 3, byrow = TRUE,
                             override.aes = list(size =3)))

sockall_er1
  
###### COMBINE SURPLUS #################################
surplus_pcas <- plot_grid(sockall_er1, NULL, sockall_lrrc9, NULL, pinkOdd_lrrc9,
                          rel_widths = c(2,0.1,2,0.1,2), 
                          nrow = 1, align = 'h',
                          labels = c('A)','','B)','','C)'), 
                          label_fontfamily = "ArialMT",
                          label_size = 30, label_colour = "black")

#jpeg(paste0("../Salmon_runtiming/2024_fourspecies/figures/pca/combine/supplfig_genes_pca_allele_cutoff_twoPanels_legend_",format(Sys.Date(),"%Y%m%d"),".jpg"),
jpeg(paste0("./figures/pca/combine/supplfig_genes_pca_allele_cutoff_twoPanels_legend_",format(Sys.Date(),"%Y%m%d"),".jpg"),
     width = 20, height = 7, res = 300, units = "in")
surplus_pcas
dev.off()
