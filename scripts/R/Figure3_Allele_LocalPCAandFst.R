# COMBINE PCA, LRRC9 LOCAL, AND ESRB LOCAL PCAS AND FST
# FIGURE 3:
#  A) LRRC9 COMBINE PCAS 
#  B) ESRB COMBINE PCAS
#  C) LRRC9 ALLELE FST
#  D) ESRB ALLELE FST

packages_needed <- c("ggplot2", "scales", "ggpubr", "ggrepel", "stringr", 
                     "data.table", "plyr","tools","gtools","reshape2", 
                     "patchwork", "cowplot", "tidyverse")

for(i in 1:length(packages_needed)){
  if(!(packages_needed[i] %in% installed.packages())){install.packages(packages_needed[i])}
  library(packages_needed[i], character.only = TRUE)
}

rm(i, packages_needed)

#setwd()

METADATAFILE <- "./data/raw/fourspecies_runtiming_metadata.csv"

#### A) LRRC9 #################################################################

##### PINK ######################

# read in the covariance matrix
pink_cov <- as.matrix(read.table("./results/pca/pink-chum_NC_068455.1_lrrc9_minInd0.3.cov"))
  pink_e <- eigen(pink_cov) # calculate eigenvector values
  pink_e_vectors <- as.data.frame(pink_e$vectors)
  pink_e_per <- pink_e$values/sum(pink_e$values) # percent explained by each component

# call in bams
pink_bam_df <- read.table("./data/bams/pink-chum_filtered_bamslist.txt", header = F)

# convert bam to FID with ABLG
pink_FID <- pink_bam_df %>%
  mutate(ABLG = as.numeric(gsub("[^0-9]","", V1))) %>% select(ABLG)

# call in some metadata
pink_meta <- read.csv(METADATAFILE, header = T) %>%
  filter(Species == 'Pink') %>%
  mutate(ABLG = as.numeric(sub('ABLG','',sampleID))) 

# join those two dataframes
pink_popFID <- inner_join(pink_FID, pink_meta, by = "ABLG")

##combine row names (population info) with the covariance matrix
pink_pca.vectors = as_tibble(cbind(pink_popFID, pink_e_vectors))

# determine the variance explained as a percent
pink_pca.eigenval.sum = sum(pink_e$values) #sum of eigenvalues
  pink_varPC1 <- (pink_e$values[1]/pink_pca.eigenval.sum)*100 # PC1 variance
  pink_varPC2 <- (pink_e$values[2]/pink_pca.eigenval.sum)*100 # PC2 variance

# Plotting
theme_set(
  theme( 
    legend.text=element_text(size=16), legend.title = element_text(size = 18),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    axis.text = element_text(angle = 0, size = 14), axis.title = element_text(size = 16),
    legend.position = "top", panel.background = element_rect(fill = "white"), 
    panel.spacing = unit(0,"lines"), strip.text.y = element_text(angle = 0)
  )
)

# create dataframe for color designation
mypalette <- c("goldenrod1", "royalblue3")
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
  ggtitle(expression('Pink'~italic(lrrc9))) +
  labs(x = paste0("PC1 (",round(pink_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(pink_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5, margin=margin(0,0,10,0)))

##### EUCLIDE SOCKEYE ####################

# read in the covariance matrix
sock_cov <- as.matrix(read.table("./results/pca/euclide_NC_068455.1_minInd0.3_lrrc9.cov"))
  sock_e <- eigen(sock_cov) # calculate eigenvector values
  sock_e_vectors <- as.data.frame(sock_e$vectors)
  sock_e_per <- sock_e$values/sum(sock_e$values) # percent explained by each component

sock_bam_df <- read.table("./data/bams/all_sockeye_bamslist.txt", header = F)

# convert bam to FID with ABLG
sock_FID <- sock_bam_df %>%
  dplyr::mutate(temp = basename(file_path_sans_ext(V1)),
                temp = gsub("^[^_]*_", "", temp), # remove everything after 1st underscore
                sampleID = str_extract(temp, "[^_]+")) %>%
  select(sampleID)

# call in some metadata
sock_meta <- read.csv(METADATAFILE, header = T) %>%
  filter(Species == "Sockeye",
         Runtime != "Late Stream") %>% # remove whitefish
  mutate(Runtime = if_else(Runtime == 'Late Beach', 'Late', 'Early'))

# join those two dataframes
sock_popFID <- inner_join(sock_FID, sock_meta, by = "sampleID")

##combine row names (population info) with the covariance matrix
sock_pca.vectors = as_tibble(cbind(sock_popFID, sock_e_vectors))

# determine the variance explained as a percent
sock_pca.eigenval.sum = sum(sock_e$values) #sum of eigenvalues
  sock_varPC1 <- (sock_e$values[1]/sock_pca.eigenval.sum)*100 # PC1 variance
  sock_varPC2 <- (sock_e$values[2]/sock_pca.eigenval.sum)*100 # PC2 variance

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
  ggtitle(expression('Sockeye'~italic(lrrc9))) +
  scale_x_reverse(breaks = c(-0.15, -0.1, -0.05, 0, 0.05, 0.1, 0.15)) +
  labs(x = paste0("PC1 (",round(sock_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(sock_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5, margin=margin(0,0,10,0)))

##### CHUM ###################

# read in the covariance matrix
chum_cov <- as.matrix(read.table("./results/pca/chumrun_NC_068455.1_lrrc9_minInd0.3_minDepthHalf.cov"))
  chum_e <- eigen(chum_cov) # calculate eigenvector values
  chum_e_vectors <- as.data.frame(chum_e$vectors)
  chum_e_per <- chum_e$values/sum(chum_e$values) # percent explained by each component

# call in bams
chum_bam_df <- read.table("./data/bams/chumrun_bamslist.txt", header = F)

# convert bam to FID with ABLG
chum_FID <- chum_bam_df %>%
  mutate(ABLG = as.numeric(gsub('[^0-9]','', V1))) %>% select(ABLG)

# call in metadata
chum_meta <- read.csv(METADATAFILE, header = T) %>%
  filter(Species == 'Chum') %>%
  mutate(ABLG = as.numeric(sub('ABLG','',sampleID))) 

# join those two dataframes
chum_popFID <- inner_join(chum_FID, chum_meta, by = "ABLG")

##combine row names (population info) with the covariance matrix
chum_pca.vectors = as_tibble(cbind(chum_popFID, chum_e_vectors))

# determine the variance explained as a percent
chum_pca.eigenval.sum = sum(chum_e$values) #sum of eigenvalues
  chum_varPC1 <- (chum_e$values[1]/chum_pca.eigenval.sum)*100 # PC1 variance
  chum_varPC2 <- (chum_e$values[2]/chum_pca.eigenval.sum)*100 # PC2 variance

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
  ggtitle(expression('Chum'~italic(lrrc9))) +
  scale_x_reverse(breaks = c(-0.05, 0, 0.05, 0.1, 0.15)) +
  labs(x = paste0("PC1 (",round(chum_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(chum_varPC2, digits = 1),"%)")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 20, hjust = 0.5, margin=margin(0,0,10,0)))


#### B) ESRB ###################################################################

###### COHO ###########

# read in the covariance matrix
coho_cov <- as.matrix(read.table("./results/pca/coho-chum_NC_068449.1_s25414060_e25501622_esrb_minInd0.3_minDepthHalf.cov"))
  coho_e <- eigen(coho_cov)
  coho_e_vectors <- as.data.frame(coho_e$vectors)
  coho_e_per <- coho_e$values/sum(coho_e$values)

coho_bam_df <- read.table("./data/bams/coho-chum_filtered_bamslist.txt", header = F)

# convert bam to FID with ABLG
coho_FID <- coho_bam_df %>%
  mutate(ABLG = as.numeric(gsub("[^0-9]", "", V1))) %>% select(-V1)

# call in some metadata
coho_meta <- read.csv(METADATAFILE, header = T) %>%
  filter(Species == 'Coho') %>%
  mutate(ABLG = as.numeric(gsub('ABLG','',sampleID)))

# join those two dataframes
coho_popFID <- inner_join(coho_FID, coho_meta, by = "ABLG")

##combine row names (population info) with the covariance matrix
coho_pca.vectors = as_tibble(cbind(coho_popFID, coho_e_vectors))

# determine the variance explained as a percent
coho_pca.eigenval.sum = sum(coho_e$values) #sum of eigenvalues
  coho_varPC1 <- (coho_e$values[1]/coho_pca.eigenval.sum)*100 # PC1 variance
  coho_varPC2 <- (coho_e$values[2]/coho_pca.eigenval.sum)*100 # PC2 variance

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
  ggtitle(expression('Coho'~italic(esrb))) +
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
  chumE_e_per <- chumE_e$values/sum(chumE_e$values) # percent explained by each component

# call in bams
chum_bam_df <- read.table("./data/bams/chumrun_bamslist.txt", header = F)

# convert bam to FID with ABLG
chum_FID <- chum_bam_df %>%
  mutate(ABLG = as.numeric(gsub("[^0-9]", "", V1))) %>% select(-V1)

# call in metadata
chum_meta <- read.csv(METADATAFILE, header = T) %>%
  filter(Species == 'Chum') %>%
  mutate(ABLG = as.numeric(sub('ABLG','',sampleID))) 

# join those two dataframes
chum_popFID <- inner_join(chum_FID, chum_meta, by = "ABLG")

##combine row names (population info) with the covariance matrix
chumE_pca.vectors = as_tibble(cbind(chum_popFID, chumE_e_vectors))

# determine the variance explained as a percent
chumE_pca.eigenval.sum = sum(chumE_e$values) #sum of eigenvalues
  chumE_varPC1 <- (chumE_e$values[1]/chumE_pca.eigenval.sum)*100 # PC1 variance
  chumE_varPC2 <- (chumE_e$values[2]/chumE_pca.eigenval.sum)*100 # PC2 variance

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
  ggtitle(expression('Chum'~italic(esrb))) +
  labs(x = paste0("PC1 (",round(chumE_varPC1, digits = 1),"%)"), 
       y= paste0("PC2 (",round(chumE_varPC2, digits = 1),"%)")) +
  theme(legend.position = "right",
        legend.text=element_text(size=32), legend.title = element_text(size=35),
        plot.title = element_text(size=20, hjust=0.5, margin=margin(0,0,10,0))) +
  guides(fill = guide_legend(title.position="top", title.hjust = 0.5,
                             override.aes = list(size=8)))

#### A/B) COMBINE ##########################################################

legend <- get_legend(chum_esrb)
chum_esrb <- chum_esrb + theme(legend.position = "none")

# cowplot of leftmost 4 plot gric
pca_four <- plot_grid(NULL,chum_lrrc9, pink_lrrc9,
                      NULL,NULL,NULL,
                      NULL,chum_esrb, coho_esrb,
                      rel_heights = c(2,0.1,2), rel_widths = c(0.09,1,1), 
                      nrow = 3, align = 'v',
                      labels = c('A','','','','','','B','',''), 
                      label_fontfamily = "helvetica",
                      label_size = 30, label_colour = "black")

# add sockeye lrrc9 and legend to the right
pca_five_legend <- plot_grid(pca_four,
                             plot_grid(sock_lrrc9, legend, ncol=1),
                             rel_widths=c(2, 1),
                             align = 'hv')
pca_five_legend

jpeg("./figures/pca/combine/fourspp_genes_pca_allele_cutoff_fivePanels_legend_20250205.jpg",
     width = 14, height = 9, res = 300, units = "in")
pca_five_legend
dev.off()

#### C) LRRC9 FST #############################################################

# Define Boundaries
  xstart.lrrc9 = 27.86
  xend.lrrc9 = 28.24
  pca.start.lrrc9 = 28128954
  pca.end.lrrc9 = 28169980

###### GFF NCBI Data ###################################

# find exons from gff file for genes of interest (from NCBI chum reference genome)
gff_df <- read.delim('../../Ref_genomes/chumV2/genomic.gff', header = F, comment.char = "#")
  gff_df <- gff_df[,c(1:5,9)] # remove excess columns
  colnames(gff_df) <- c("chrName", "RefSeq","exon","start.pos","fin.pos", "ID")

# only keep chr29 & 35
gff_chr35 <- gff_df %>%
  filter(chrName == "NC_068455.1")
gff_chr29 <- gff_df %>%
  filter(chrName == "NC_068449.1")
rm(gff_df)

########## Pink lrrc9 FST ###################
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

##### Chum lrrc9 FST ######################
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

############ Add Genes ########################################################

# only the region of interest
gff_region <- gff_chr35 %>%
  mutate(start.pos = start.pos/1e6, fin.pos = fin.pos/1e6) %>% 
  filter(fin.pos > xstart.lrrc9, start.pos < xend.lrrc9)

# prep pattern for str_match below
gene_pattern <- "gene=\\s*(.*?)\\s*;"           # keep string btwn "gene=" & ":product" 
exon_pattern <- "ID=exon-\\s*(.*?)\\s*;Parent"  # keep string btwn "exon=" & ";Parent" 
descr_pattern <- ";description=\\s*(.*?)\\s*;"  # keep description

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

##### Make genes with FST < 0.5 gray w/o legend ####################
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

################# Plotting ###########################
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


#### D) ESRB FST ##############################################################

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

##### Coho esrb FST ###########################

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

##### Make genes with FST < 0.5 gray w/o legend ####################
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


#### C/D) Combine #####################
twofst_cowplot <- plot_grid(multiplot_lrrc9, multiplot_esrb,
                            rel_widths = c(5,5), ncol = 2, nrow = 1,
                            labels = c('C','D'), label_fontfamily = "helvetica",
                            label_size = 30, label_colour = "black")
twofst_cowplot

#### FIG 3) Combine A-D #####################

fig3 <- plot_grid(pca_five_legend, NULL, twofst_cowplot,
                  rel_heights = c(1,0.1, 1), 
                  ncol = 1, nrow = 3, align = 'v')
fig3

jpeg(paste0("./figures/figure3_panelsA-D_",format(Sys.Date(),"%Y%m%d"),".jpg"), 
     width = 18, height = 20, res = 200, units = "in")
print(fig3)
dev.off()
