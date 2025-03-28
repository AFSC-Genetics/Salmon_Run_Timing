# Supplemental local PCAs for Sockeye (Whitefish) and Pink (Odd)

library(here)
library(stringr)
library(tools)
library(tidyverse)
library(patchwork)
library(cowplot)

setwd("C:/Users/Natasha.Howe/Work/Salmon-Genetics/Salmon_runtiming/2024_fourspecies")

##### Whitefish Addition ############
# read in the covariance matrix
sockall_cov <- as.matrix(read.table("../2024_sockeye/results/pca/sock-all_NC_068455.1_minInd0.3_lrrc9.cov"))
  sockall_e <- eigen(sockall_cov) # calculate eigenvector values
  sockall_e_vectors <- as.data.frame(sockall_e$vectors)
  sockall_e_per <- sockall_e$values/sum(sockall_e$values) # percent explained by each component

sockall_bam_df <- read.table("../2024_sockeye/sedna_files/bams/all_sockeye_bamslist.txt", header = F)

# convert bam to FID with ABLG
sockall_FID <- sockall_bam_df %>%
  dplyr::mutate(temp = basename(file_path_sans_ext(V1)),
                temp = gsub("^[^_]*_", "", temp), # remove everything after 1st underscore
                sampleID = str_extract(temp, "[^_]+")) %>%
  select(-c(V1, temp))

# call in some metadata
sockall_meta <- read.csv("../2024_sockeye/data/raw/wood_runtiming_metadata_anvil_added.csv", 
                         header = T) 

# join those two dataframes
sockall_popFID <- inner_join(sockall_FID, sockall_meta, by = "sampleID") %>%
  mutate(Runtime = sub('Creek','Stream',Runtime))

##combine row names (population info) with the covariance matrix
sockall_pca.vectors = as_tibble(cbind(sockall_popFID, sockall_e_vectors))

# determine the variance explained as a percent
sockall_pca.eigenval.sum = sum(sockall_e$values) #sum of eigenvalues
  sockall_varPC1 <- (sockall_e$values[1]/sockall_pca.eigenval.sum)*100 #Variance explained by PC1
  sockall_varPC2 <- (sockall_e$values[2]/sockall_pca.eigenval.sum)*100 #Variance explained by PC2

# create dataframe for color designation
#sockall_df_popMetaData <- data.frame(Runtime = )

sockall_Palette <- c("goldenrod1", "royalblue3", "orchid")
sockall_pca.vectors$Runtime <- factor(sockall_pca.vectors$Runtime, levels = c("Early Stream", "Late Stream", "Late Beach"))
names(sockall_Palette) <- c("Early Stream", "Late Stream", "Late Beach")

theme_set(
  theme( 
    legend.text=element_text(size=16),
    legend.title = element_text(size = 18),
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(angle = 0, size = 14),
    axis.title = element_text(size = 16),
    panel.background = element_rect(fill = "white"), 
    panel.spacing = unit(0,"lines"),
    strip.text.y = element_text(angle = 0)
  )
)

sockall_lrrc9 <- ggplot(data = sockall_pca.vectors, 
                        aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = 3, color = "gray20") +
  scale_fill_manual(name = expression('Sockeye'~italic(lrrc9)), values = sockall_Palette) +
  scale_shape_manual(name = expression('Sockeye'~italic(lrrc9)), values = c(21,22,23)) +
  geom_vline(xintercept = 0.1, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = -0.05, color = "gray30", alpha = 0.5, linetype = "dashed") +
  scale_x_reverse(breaks = c(-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15)) +
  labs(x = paste0("PC1 (",round(sockall_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(sockall_varPC2, digits = 1),"%)")) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5,
                             override.aes = list(size =3)))
sockall_lrrc9

##### Odd-Pink Addition ############

# read in the covariance matrix
pinkO_cov <- as.matrix(read.table("../2024_pinkOdd/results/pca/pink-odd_NC_068455.1_lrrc9_minInd0.3.cov"))
  pinkO_e <- eigen(pinkO_cov)
  pinkO_e_vectors <- as.data.frame(pinkO_e$vectors)
  pinkO_e_per <- pinkO_e$values/sum(pinkO_e$values) # percent explained by each component

# call in bams
pinkO_bam_df <- read.table("../2024_pinkOdd/sedna_files/bams/pink-odd_filtered_bamslist.txt", header = F)

# convert bam to FID with ABLG
pinkO_FID <- pinkO_bam_df %>%
  mutate(V1 = basename(file_path_sans_ext(V1)),
         temp = gsub("^[^_]*_", "", V1), # remove everything after 1st underscore
         sampleID = str_extract(temp, "[^_]+")) %>%
  select(-c(V1, temp))

# call in some metadata
pinkO_meta <- read.csv("../2024_pinkOdd/data/raw/pinkOdd_collection_07172024_cleaned.csv", 
                       header = T) %>%
  mutate(sampleID = paste0("ABLG",ABLG),
         Runtime = ifelse(Runtime == "late","Late","Early"))

# join those two dataframes
pinkO_popFID <- inner_join(pinkO_FID, pinkO_meta, by = "sampleID")

##combine row names (population info) with the covariance matrix
pinkO_pca.vectors = as_tibble(cbind(pinkO_popFID, pinkO_e_vectors))

# determine the variance explained as a percent
pinkO_pca.eigenval.sum = sum(pinkO_e$values) #sum of eigenvalues
  pinkO_varPC1 <- (pinkO_e$values[1]/pinkO_pca.eigenval.sum)*100 #Variance explained by PC1
  pinkO_varPC2 <- (pinkO_e$values[2]/pinkO_pca.eigenval.sum)*100 #Variance explained by PC2

df_popMetaData <- data.frame(Runtime = c("Early", "Late"),color = c("goldenrod1", "royalblue3"))

mypalette <- c("goldenrod1", "royalblue3")
pinkO_pca.vectors$Runtime <- factor(pinkO_pca.vectors$Runtime, levels = c("Early", "Late"))
names(mypalette) <- levels(pinkO_pca.vectors$Runtime)

### PINK ODD
pinkOdd_lrrc9 <- ggplot(data = pinkO_pca.vectors, 
                        aes(x=V1, y=V2, fill = Runtime, shape = Runtime)) + 
  geom_point(alpha = 0.7, size = 3, color = "gray20") +
  scale_fill_manual(name = expression('Pink-Odd'~italic(lrrc9)), values = mypalette) +
  scale_shape_manual(name = expression('Pink-Odd'~italic(lrrc9)), values = c(21,22)) +
  geom_vline(xintercept = 0, color = "gray30", alpha = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.15, color = "gray30", alpha = 0.5, linetype = "dashed") +
  scale_x_reverse(breaks = c(-0.05, 0, 0.05, 0.1, 0.15, 0.2)) +
  labs(x = paste0("PC1 (",round(pinkO_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(pinkO_varPC2, digits = 1),"%)")) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5,
                             override.aes = list(size =3)))
pinkOdd_lrrc9

###### COMBINE SURPLUS #################################
surplus_pcas <- plot_grid(pinkOdd_lrrc9, NULL, sockall_lrrc9,
                          rel_widths = c(1,0.1,1), 
                          ncol = 3, align = 'h',
                          labels = c('A','','B'), 
                          label_fontfamily = "ArialMT",
                          label_size = 30, label_colour = "black")

jpeg(paste0("./figures/pca/combine/supplfig_twospp_genes_pca_allele_cutoff_twoPanels_legend_",format(Sys.Date(),"%Y%m%d"),".jpg"),
     width = 14, height = 7, res = 300, units = "in")
surplus_pcas
dev.off()
